import glob
from collections import OrderedDict
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from Bio.PDB.Polypeptide import index_to_one, index_to_three, one_to_index, three_to_one
from scipy.stats import pearsonr
from torch.nn.functional import softmax
from torch.utils.data import DataLoader, Dataset

from cavity_model import (
    CavityModel,
    DDGDataset,
    DDGToTensor,
    DownstreamModel,
    DownstreamModelSimple,
    ResidueEnvironmentsDataset,
    ToTensor,
)


def train_val_split(
    parsed_pdb_filenames: List[str],
    TRAIN_VAL_SPLIT: float,
    DEVICE: str,
    BATCH_SIZE: int,
):
    """
    Helper function to perform training and validation split of ResidueEnvironments. Note that
    we do the split on PDB level not on ResidueEnvironment level due to possible leakage.
    """
    n_train_pdbs = int(len(parsed_pdb_filenames) * TRAIN_VAL_SPLIT)
    filenames_train = parsed_pdb_filenames[:n_train_pdbs]
    filenames_val = parsed_pdb_filenames[n_train_pdbs:]

    to_tensor_transformer = ToTensor(DEVICE)

    dataset_train = ResidueEnvironmentsDataset(filenames_train, transformer=to_tensor_transformer)
    dataset_val = ResidueEnvironmentsDataset(filenames_val, transformer=to_tensor_transformer)

    dataloader_train = DataLoader(
        dataset_train,
        batch_size=BATCH_SIZE,
        shuffle=True,
        collate_fn=to_tensor_transformer.collate_cat,
        drop_last=True,
    )
    # TODO: Fix it so drop_last doesn't have to be True when calculating validation accuracy.
    dataloader_val = DataLoader(
        dataset_val,
        batch_size=BATCH_SIZE,
        shuffle=False,
        collate_fn=to_tensor_transformer.collate_cat,
        drop_last=True,
    )

    print(
        f"Training data set includes {len(filenames_train)} pdbs with "
        f"{len(dataset_train)} environments."
    )
    print(
        f"Validation data set includes {len(filenames_val)} pdbs with "
        f"{len(dataset_val)} environments."
    )

    return dataloader_train, dataset_train, dataloader_val, dataset_val


def _train_step(
    batch_x: torch.Tensor,
    batch_y: torch.Tensor,
    cavity_model_net: CavityModel,
    optimizer: torch.optim.Adam,
    loss_function: torch.nn.CrossEntropyLoss,
) -> Tuple[torch.Tensor, float]:
    """
    Helper function to take a training step
    """
    cavity_model_net.train()
    optimizer.zero_grad()
    batch_y_pred = cavity_model_net(batch_x)
    loss_batch = loss_function(batch_y_pred, torch.argmax(batch_y, dim=-1))
    loss_batch.backward()
    optimizer.step()
    return (batch_y_pred, loss_batch.detach().cpu().item())


def _eval_loop(
    cavity_model_net: CavityModel,
    dataloader_val: DataLoader,
    loss_function: torch.nn.CrossEntropyLoss,
) -> Tuple[float, float]:
    """
    Helper function to perform an eval loop
    """
    # Eval loop. Due to memory, we don't pass the whole eval set to the model
    labels_true_val = []
    labels_pred_val = []
    loss_batch_list_val = []
    for batch_x_val, batch_y_val in dataloader_val:
        cavity_model_net.eval()
        batch_y_pred_val = cavity_model_net(batch_x_val)

        loss_batch_val = loss_function(batch_y_pred_val, torch.argmax(batch_y_val, dim=-1))
        loss_batch_list_val.append(loss_batch_val.detach().cpu().item())

        labels_true_val.append(torch.argmax(batch_y_val, dim=-1).detach().cpu().numpy())
        labels_pred_val.append(torch.argmax(batch_y_pred_val, dim=-1).detach().cpu().numpy())
    acc_val = np.mean((np.reshape(labels_true_val, -1) == np.reshape(labels_pred_val, -1)))
    loss_val = np.mean(loss_batch_list_val)
    return acc_val, loss_val


def train_loop(
    dataloader_train: DataLoader,
    dataloader_val: DataLoader,
    cavity_model_net: CavityModel,
    loss_function: torch.nn.CrossEntropyLoss,
    optimizer: torch.optim.Adam,
    prefix: str,
    EPOCHS: int,
    PATIENCE_CUTOFF: int,
):
    """
    Helper function to perform training loop for the Cavity Model.
    """
    current_best_epoch_idx = -1
    current_best_loss_val = 1e4
    patience = 0
    epoch_idx_to_model_path = {}
    for epoch in range(EPOCHS):
        labels_true = []
        labels_pred = []
        loss_batch_list = []
        for batch_x, batch_y in dataloader_train:
            # Take train step
            batch_y_pred, loss_batch = _train_step(
                batch_x, batch_y, cavity_model_net, optimizer, loss_function
            )
            loss_batch_list.append(loss_batch)

            labels_true.append(torch.argmax(batch_y, dim=-1).detach().cpu().numpy())
            labels_pred.append(torch.argmax(batch_y_pred, dim=-1).detach().cpu().numpy())

        # Train epoch metrics
        acc_train = np.mean((np.reshape(labels_true, -1) == np.reshape(labels_pred, -1)))
        loss_train = np.mean(loss_batch_list)

        # Validation epoch metrics
        acc_val, loss_val = _eval_loop(cavity_model_net, dataloader_val, loss_function)

        print(
            f"Epoch {epoch:2d}. Train loss: {loss_train:5.3f}. "
            f"Train Acc: {acc_train:4.2f}. Val loss: {loss_val:5.3f}. "
            f"Val Acc {acc_val:4.2f}"
        )

        # Save model
        model_path = f"cavity_models/{prefix}_model_epoch_{epoch:02d}.pt"
        epoch_idx_to_model_path[epoch] = model_path
        torch.save(cavity_model_net.state_dict(), model_path)

        # Early stopping
        if loss_val < current_best_loss_val:
            current_best_loss_val = loss_val
            current_best_epoch_idx = epoch
            patience = 0
        else:
            patience += 1
        if patience > PATIENCE_CUTOFF:
            print("Early stopping activated.")
            break

    best_model_path = epoch_idx_to_model_path[current_best_epoch_idx]
    print(
        f"Best epoch idx: {current_best_epoch_idx} with validation loss: "
        f"{current_best_loss_val:5.3f} and model_path: "
        f"{best_model_path}"
    )
    return best_model_path


def populate_dfs_with_resenvs(
    ddg_data_dict: Dict[str, pd.DataFrame], resenv_datasets_look_up: Dict[str, Dataset]
):
    """
    Helper function populate ddG dfs with the WT ResidueEnvironment objects.
    """
    print(
        "Dropping data points where residue is not defined in structure "
        "or due to missing parsed pdb file"
    )
    # Add wt residue environments to standard ddg data dataframes
    for ddg_data_key in ddg_data_dict.keys():
        resenvs_ddg_data = []
        for idx, row in ddg_data_dict[ddg_data_key].iterrows():
            resenv_key = (
                f"{row['pdbid']}{row['chainid']}_" f"{row['variant'][1:-1]}{row['variant'][0]}"
            )
            try:
                if "symmetric" in ddg_data_key:
                    ddg_data_key_adhoc_fix = "symmetric"
                else:
                    ddg_data_key_adhoc_fix = ddg_data_key
                resenv = resenv_datasets_look_up[ddg_data_key_adhoc_fix][resenv_key]
                resenvs_ddg_data.append(resenv)
            except KeyError:
                resenvs_ddg_data.append(np.nan)
        ddg_data_dict[ddg_data_key]["resenv"] = resenvs_ddg_data
        n_datapoints_before = ddg_data_dict[ddg_data_key].shape[0]
        ddg_data_dict[ddg_data_key].dropna(inplace=True)
        n_datapoints_after = ddg_data_dict[ddg_data_key].shape[0]
        print(
            f"dropped {n_datapoints_before - n_datapoints_after:4d} / "
            f"{n_datapoints_before:4d} data points from dataset {ddg_data_key}"
        )

        # Add wt and mt idxs to df
        ddg_data_dict[ddg_data_key]["wt_idx"] = ddg_data_dict[ddg_data_key].apply(
            lambda row: one_to_index(row["variant"][0]), axis=1
        )
        ddg_data_dict[ddg_data_key]["mt_idx"] = ddg_data_dict[ddg_data_key].apply(
            lambda row: one_to_index(row["variant"][-1]), axis=1
        )


def populate_dfs_with_nlls_and_nlfs(
    ddg_data_dict: Dict[str, pd.DataFrame],
    cavity_model_infer_nets: List[CavityModel],
    DEVICE: str,
    BATCH_SIZE: int,
    EPS: float,
):
    """
    Helper function to populate ddG dfs with predicted negative-log-likelihoods and negative-log-frequencies
    """

    # Load PDB amino acid frequencies used to approximate unfolded states
    pdb_nlfs = -np.log(np.load("data/pdb_frequencies.npz")["frequencies"])

    # Load IDP amino acid frequences to approximate unfolded states to dict
    idp_nlfs = {
        "A": -np.log(0.06856891151135473),
        "C": -np.log(0.03244909945184025),
        "D": -np.log(0.058633516053249804),
        "E": -np.log(0.09025058731401722),
        "F": -np.log(0.023296789350039156),
        "G": -np.log(0.06528974158183241),
        "H": -np.log(0.009641738449490995),
        "I": -np.log(0.04113645262333594),
        "K": -np.log(0.07845536413469067),
        "L": -np.log(0.07852877838684416),
        "M": -np.log(0.018794048551292093),
        "N": -np.log(0.04023101018010963),
        "P": -np.log(0.07703602192638997),
        "Q": -np.log(0.047841620986687546),
        "R": -np.log(0.0668559122944401),
        "S": -np.log(0.07938527799530148),
        "T": -np.log(0.05716523101018011),
        "V": -np.log(0.04835552075176194),
        "W": -np.log(0.003205755677368833),
        "Y": -np.log(0.014878621769772905),
    }

    # Add predicted Nlls and NLFs to ddG dataframes
    for ddg_data_key in ddg_data_dict.keys():
        df = ddg_data_dict[ddg_data_key]

        # Perform predictions on matched residue environments
        ddg_resenvs = list(df["resenv"].values)
        ddg_resenv_dataset = ResidueEnvironmentsDataset(ddg_resenvs, transformer=ToTensor(DEVICE))

        # Define dataloader for resenvs matched to ddG data
        ddg_resenv_dataloader = DataLoader(
            ddg_resenv_dataset,
            batch_size=BATCH_SIZE,
            shuffle=False,
            collate_fn=ToTensor(DEVICE).collate_cat,
            drop_last=False,
        )

        for i, cavity_model_i in enumerate(cavity_model_infer_nets):
            # Infer NLLs
            pred_nlls = []
            for batch_x, _ in ddg_resenv_dataloader:
                batch_pred_nlls = (
                    -torch.log(softmax(cavity_model_i(batch_x), dim=-1) + EPS)
                    .detach()
                    .cpu()
                    .numpy()
                )
                pred_nlls.append(batch_pred_nlls)
            pred_nlls_list = [row for row in np.vstack(pred_nlls)]

            # Add NLLs to dataframe
            df[f"nlls_{i}"] = pred_nlls_list

            # Isolate WT and MT NLLs and add to datafra
            df[f"wt_nll_{i}"] = df.apply(lambda row: row[f"nlls_{i}"][row["wt_idx"]], axis=1)
            df[f"mt_nll_{i}"] = df.apply(lambda row: row[f"nlls_{i}"][row["mt_idx"]], axis=1)

            # Add PDB database statistics negative log frequencies to df
            df[f"wt_nlf_{i}"] = df.apply(lambda row: pdb_nlfs[row["wt_idx"]], axis=1)
            df[f"mt_nlf_{i}"] = df.apply(lambda row: pdb_nlfs[row["mt_idx"]], axis=1)

            # Add vanilla predictions without taking into acount unfolded state
            df[f"ddg_pred_ultra_vanilla_{i}"] = df.apply(
                lambda row: row[f"mt_nll_{i}"] - row[f"wt_nll_{i}"],
                axis=1,
            )

            # Add ddG prediction (without downstream model)
            df[f"ddg_pred_no_ds_{i}"] = df.apply(
                lambda row: row[f"mt_nll_{i}"]
                - row[f"mt_nlf_{i}"]
                - row[f"wt_nll_{i}"]
                + row[f"wt_nlf_{i}"],
                axis=1,
            )

            # Add IDP database statistics negative log frequencies to df
            df[f"wt_idp_nlf_{i}"] = df.apply(
                lambda row: idp_nlfs[index_to_one(row["wt_idx"])], axis=1
            )
            df[f"mt_idp_nlf_{i}"] = df.apply(
                lambda row: idp_nlfs[index_to_one(row["mt_idx"])], axis=1
            )

            # Add ddG prediction with IDP statistics (without downstream model)
            df[f"ddg_pred_idp_no_ds_{i}"] = df.apply(
                lambda row: row[f"mt_nll_{i}"]
                - row[f"mt_idp_nlf_{i}"]
                - row[f"wt_nll_{i}"]
                + row[f"wt_idp_nlf_{i}"],
                axis=1,
            )


def augment_with_reverse_mutation(ddg_data_dict: Dict[str, pd.DataFrame]):
    """
    Helper function that augments the ddg dfs with the reverse mutations.
    The dict contains deep copies of the original dataframes before copying.
    """

    ddg_data_dict_augmented = OrderedDict()
    for ddg_key in ["dms", "protein_g", "guerois"]:
        ddg_data_df_augmented = ddg_data_dict[ddg_key].copy(deep=True).drop(columns="resenv")
        rows_augmented = []
        for row in ddg_data_df_augmented.iterrows():
            row_cp = row[1].copy(deep=True)
            row_cp.name = str(row_cp.name) + "_augmented"

            # Augment
            row_cp["variant"] = (
                row_cp["variant"][-1] + row_cp["variant"][1:-1] + row_cp["variant"][0]
            )
            row_cp["ddg"] = -1.0 * row_cp["ddg"]
            row_cp["ddg_pred_no_ds_0"] = -1.0 * row_cp["ddg_pred_no_ds_0"]

            (
                wt_idx,
                mt_idx,
                wt_nll,
                mt_nll,
                wt_nlf,
                mt_nlf,
                wt_nll_md,
                fragment_nll_wt_given_mt,
                fragment_nll_wt_given_wt,
                fragment_nll_mt_given_mt,
                fragment_nll_mt_given_wt,
            ) = (
                row_cp["mt_idx"],
                row_cp["wt_idx"],
                row_cp["mt_nll_0"],
                row_cp["wt_nll_0"],
                row_cp["mt_nlf_0"],
                row_cp["wt_nlf_0"],
                row_cp["mt_nll_md_0"],
                row_cp["fragment_nll_mt_given_wt_0"],
                row_cp["fragment_nll_mt_given_mt_0"],
                row_cp["fragment_nll_wt_given_wt_0"],
                row_cp["fragment_nll_wt_given_mt_0"],
            )

            (
                row_cp["wt_idx"],
                row_cp["mt_idx"],
                row_cp["wt_nll_0"],
                row_cp["mt_nll_0"],
                row_cp["wt_nlf_0"],
                row_cp["mt_nlf_0"],
                row_cp["wt_nll_md_0"],
                row_cp["fragment_nll_wt_given_mt_0"],
                row_cp["fragment_nll_wt_given_wt_0"],
                row_cp["fragment_nll_mt_given_mt_0"],
                row_cp["fragment_nll_mt_given_wt_0"],
            ) = (
                wt_idx,
                mt_idx,
                wt_nll,
                mt_nll,
                wt_nlf,
                mt_nlf,
                wt_nll_md,
                fragment_nll_wt_given_mt,
                fragment_nll_wt_given_wt,
                fragment_nll_mt_given_mt,
                fragment_nll_mt_given_wt,
            )

            rows_augmented.append(row_cp)
        ddg_data_df_augmented = ddg_data_df_augmented.append(rows_augmented)
        ddg_data_dict_augmented[ddg_key] = ddg_data_df_augmented

    return ddg_data_dict_augmented


def get_ddg_training_dataloaders(ddg_data_dict_augmented, BATCH_SIZE_DDG, SHUFFLE_DDG, transformer):
    ddg_dataloaders_train_dict = {}
    for key in ddg_data_dict_augmented.keys():
        ddg_dataset_aug = DDGDataset(ddg_data_dict_augmented[key], transformer=transformer())
        ddg_dataloader_aug = DataLoader(
            ddg_dataset_aug,
            batch_size=BATCH_SIZE_DDG,
            shuffle=SHUFFLE_DDG,
            drop_last=True,
        )
        ddg_dataloaders_train_dict[key] = ddg_dataloader_aug

    return ddg_dataloaders_train_dict


def get_ddg_validation_dataloaders(
    ddg_data_dict, transformer, keys=["dms", "protein_g", "guerois"]
):
    """
    Helper function that return validation set dataloaders for ddg data.
    """
    ddg_dataloaders_val_dict = {}
    for key in keys:
        ddg_dataset = DDGDataset(ddg_data_dict[key], transformer=transformer())
        ddg_dataloader = DataLoader(
            ddg_dataset,
            batch_size=len(ddg_dataset),  # Full dataset as batch size
            shuffle=False,
            drop_last=False,
        )
        ddg_dataloaders_val_dict[key] = ddg_dataloader

    return ddg_dataloaders_val_dict

import sklearn.metrics

def train_downstream_and_evaluate(
    ddg_dataloaders_train_dict,
    ddg_dataloaders_val_dict,
    DEVICE,
    LEARNING_RATE_DDG,
    EPOCHS_DDG,
):
    pearsons_r_results_dict = {}
    msa_results_dict = {}
    for train_key in ["dms", "protein_g", "guerois"]:
        # Define model
        downstream_model_net = DownstreamModel().to(DEVICE)
        loss_ddg = torch.nn.MSELoss()
        optimizer_ddg = torch.optim.Adam(downstream_model_net.parameters(), lr=LEARNING_RATE_DDG)
        for epoch in range(EPOCHS_DDG):
            
            # First epoch evaluation
            if epoch == 0:
                for val_key in ["dms", "protein_g", "guerois"]:
#                     if val_key == train_key:
#                         continue

                    val_batch = next(iter(ddg_dataloaders_val_dict[val_key]))
                    val_batch_x, val_batch_y = (
                        val_batch["x_"].to(DEVICE),
                        val_batch["y_"],
                    )
                    val_batch_y_pred = (
                        downstream_model_net(val_batch_x).reshape(-1).detach().cpu().numpy()
                    )
                    pearson_r = pearsonr(val_batch_y_pred, val_batch_y)[0]
                    msa_all = sklearn.metrics.mean_absolute_error(val_batch_y_pred, val_batch_y)
                    if train_key not in pearsons_r_results_dict:
                        pearsons_r_results_dict[train_key] = {}
                        msa_results_dict[train_key] = {}
                    if val_key not in pearsons_r_results_dict[train_key]:
                        pearsons_r_results_dict[train_key][val_key] = []
                        msa_results_dict[train_key][val_key]= []
 
                    pearsons_r_results_dict[train_key][val_key].append(pearson_r)
                    msa_results_dict[train_key][val_key].append(msa_all)
                
            
            for batch in ddg_dataloaders_train_dict[train_key]:
                x_batch, y_batch = batch["x_"].to(DEVICE), batch["y_"].to(DEVICE)

                downstream_model_net.train()
                optimizer_ddg.zero_grad()
                y_batch_pred = downstream_model_net(x_batch).squeeze()
                loss_ddg_batch = loss_ddg(y_batch_pred, y_batch)
                loss_ddg_batch.backward()
                optimizer_ddg.step()

            # Evaluation
            for val_key in ["dms", "protein_g", "guerois"]:
#                     if val_key == train_key:
#                         continue

                val_batch = next(iter(ddg_dataloaders_val_dict[val_key]))
                val_batch_x, val_batch_y = (
                    val_batch["x_"].to(DEVICE),
                    val_batch["y_"],
                )
                val_batch_y_pred = (
                    downstream_model_net(val_batch_x).reshape(-1).detach().cpu().numpy()
                )
                pearson_r = pearsonr(val_batch_y_pred, val_batch_y)[0]
                msa_all = sklearn.metrics.mean_absolute_error(val_batch_y_pred, val_batch_y)
                if train_key not in pearsons_r_results_dict:
                    pearsons_r_results_dict[train_key] = {}
                    msa_results_dict[train_key] = {}
                if val_key not in pearsons_r_results_dict[train_key]:
                    pearsons_r_results_dict[train_key][val_key] = []
                    msa_results_dict[train_key][val_key]= []

                pearsons_r_results_dict[train_key][val_key].append(pearson_r)
                msa_results_dict[train_key][val_key].append(msa_all)

            if epoch % 10 == 0:
                print(train_key, f"epoch {epoch:3d}")
    return pearsons_r_results_dict, msa_results_dict


def train_downstream_and_evaluate_simple(
    ddg_dataloaders_train_dict,
    ddg_dataloaders_val_dict,
    DEVICE,
    LEARNING_RATE_DDG,
    EPOCHS_DDG,
):
    pearsons_r_results_dict = {}
    msa_results_dict = {}
    for train_key in ["dms", "protein_g", "guerois"]:
        # Define model
        downstream_model_net = DownstreamModelSimple().to(DEVICE)
        loss_ddg = torch.nn.MSELoss()
        optimizer_ddg = torch.optim.Adam(downstream_model_net.parameters(), lr=LEARNING_RATE_DDG)
        for epoch in range(EPOCHS_DDG):
            
            # First epoch evaluation
            if epoch == 0:
                for val_key in ["dms", "protein_g", "guerois"]:
#                     if val_key == train_key:
#                         continue

                    val_batch = next(iter(ddg_dataloaders_val_dict[val_key]))
                    val_batch_x, val_batch_y = (
                        val_batch["x_"].to(DEVICE),
                        val_batch["y_"],
                    )
                    val_batch_y_pred = (
                        downstream_model_net(val_batch_x).reshape(-1).detach().cpu().numpy()
                    )
                    pearson_r = pearsonr(val_batch_y_pred, val_batch_y)[0]
                    msa_all = sklearn.metrics.mean_absolute_error(val_batch_y_pred, val_batch_y)
                    if train_key not in pearsons_r_results_dict:
                        pearsons_r_results_dict[train_key] = {}
                        msa_results_dict[train_key] = {}
                    if val_key not in pearsons_r_results_dict[train_key]:
                        pearsons_r_results_dict[train_key][val_key] = []
                        msa_results_dict[train_key][val_key]= []
 
                    pearsons_r_results_dict[train_key][val_key].append(pearson_r)
                    msa_results_dict[train_key][val_key].append(msa_all)
                
            
            for batch in ddg_dataloaders_train_dict[train_key]:
                x_batch, y_batch = batch["x_"].to(DEVICE), batch["y_"].to(DEVICE)

                downstream_model_net.train()
                optimizer_ddg.zero_grad()
                y_batch_pred = downstream_model_net(x_batch).squeeze()
                loss_ddg_batch = loss_ddg(y_batch_pred, y_batch)
                loss_ddg_batch.backward()
                optimizer_ddg.step()

            # Evaluation
            for val_key in ["dms", "protein_g", "guerois"]:
#                     if val_key == train_key:
#                         continue

                val_batch = next(iter(ddg_dataloaders_val_dict[val_key]))
                val_batch_x, val_batch_y = (
                    val_batch["x_"].to(DEVICE),
                    val_batch["y_"],
                )
                val_batch_y_pred = (
                    downstream_model_net(val_batch_x).reshape(-1).detach().cpu().numpy()
                )
                pearson_r = pearsonr(val_batch_y_pred, val_batch_y)[0]
                msa_all = sklearn.metrics.mean_absolute_error(val_batch_y_pred, val_batch_y)
                if train_key not in pearsons_r_results_dict:
                    pearsons_r_results_dict[train_key] = {}
                    msa_results_dict[train_key] = {}
                if val_key not in pearsons_r_results_dict[train_key]:
                    pearsons_r_results_dict[train_key][val_key] = []
                    msa_results_dict[train_key][val_key]= []

                pearsons_r_results_dict[train_key][val_key].append(pearson_r)
                msa_results_dict[train_key][val_key].append(msa_all)

            if epoch % 10 == 0:
                print(train_key, f"epoch {epoch:3d}")
    return pearsons_r_results_dict, msa_results_dict



def _trim_left_flank(left_flank: str):
    if len(left_flank) <= 5:
        padding = "-" * (5 - len(left_flank))
        return padding + left_flank
    else:
        return left_flank[-5:]


def _trim_right_flank(right_flank: str):
    if len(right_flank) <= 5:
        padding = "-" * (5 - len(right_flank))
        return right_flank + padding
    else:
        return right_flank[:5]


def add_flanking_seq_fragments(ddg_data_dict: Dict, dataset: str, pdb_filename: str):

    if "left_flank" not in ddg_data_dict[dataset].columns:
        ddg_data_dict[dataset]["left_flank"] = np.nan
    if "wt_restype" not in ddg_data_dict[dataset].columns:
        ddg_data_dict[dataset]["wt_restype"] = np.nan
    if "mt_restype" not in ddg_data_dict[dataset].columns:
        ddg_data_dict[dataset]["mt_restype"] = np.nan
    if "right_flank" not in ddg_data_dict[dataset].columns:
        ddg_data_dict[dataset]["right_flank"] = np.nan

    pdbid = pdb_filename.split(r"/")[-1][0:4].upper()

    # # Load SEQRES
    # chain_id_to_seq_res = {}
    # for record in SeqIO.parse(pdb_filename, "pdb-seqres"):
    #     seq_res = str(record.seq)
    #     chain_id = record.id[-1]
    #     chain_id_to_seq_res[chain_id] = seq_res
    #     print(record.annotations)

    # # Load PDBSEQ
    # from Bio.SeqIO.PdbIO import PdbAtomIterator
    # chain_id_to_pdb_seq = {}
    # with open(pdb_filename) as handle:
    #     for record in PdbAtomIterator(handle):
    #         pdb_seq = str(record.seq)
    #         chain_id = record.id[-1]
    #         chain_id_to_pdb_seq[chain_id] = pdb_seq

    from Bio.PDB.PDBParser import PDBParser

    p = PDBParser()
    model_first = p.get_structure(pdbid, pdb_filename)[0]
    chain_id_to_pdb_seq = {}
    chain_id_to_pdb_residue_numbers = {}
    for chain in model_first:
        pdb_seq = []
        pdb_residue_numbers = []
        for residue in chain.get_residues():
            if residue.resname.strip() in [index_to_three(i) for i in range(20)]:
                pdb_residue_numbers.append(residue.id[1])
                pdb_seq.append(three_to_one(residue.resname.strip()))
        chain_id_to_pdb_seq[chain.id] = "".join(pdb_seq)
        chain_id_to_pdb_residue_numbers[chain.id] = pdb_residue_numbers

    for idx, row in ddg_data_dict[dataset].iterrows():
        if row["pdbid"] == pdbid:
            residue_number = int(row["variant"][1:-1])
            chain_id = row["chainid"]

            pdb_sequence = chain_id_to_pdb_seq[chain_id]
            resid = chain_id_to_pdb_residue_numbers[chain_id].index(residue_number)

            if row["variant"][0] == pdb_sequence[resid]:
                ddg_data_dict[dataset].loc[idx, "left_flank"] = _trim_left_flank(
                    pdb_sequence[:resid]
                )
                ddg_data_dict[dataset].loc[idx, "right_flank"] = _trim_right_flank(
                    pdb_sequence[resid + 1 :]
                )
                ddg_data_dict[dataset].loc[idx, "wt_restype"] = row["variant"][0]
                ddg_data_dict[dataset].loc[idx, "mt_restype"] = row["variant"][-1]
            else:
                print("WRONG", row[["pdbid", "variant"]])


# TODO Rename probabilities to NLL
def infer_probabilities_for_center_residues(
    ddg_data_dict: Dict,
    dataset: str,
    cavity_model_infer_nets: CavityModel,
    DEVICE: str,
    EPS: float,
    *,
    is_wt: bool = True,
    stride: int = 50,
):
    """
    Add the inferred wt and mt probabilities for the center residues of simulated sequence
    fragments. This is used for modelling the unfolded state.
    """

    for j, cavity_model_infer_net in enumerate(cavity_model_infer_nets):
        # Lets get the probabilities for the center residues
        fragments_nlls_wt = []
        fragments_nlls_mt = []
        for i, row in ddg_data_dict[dataset].iterrows():
            pdb_id = row["pdbid"]
            print(i, pdb_id)

            if is_wt:
                npz_wc = (
                    f"data/data_{dataset}/simulate_seq_fragments_{dataset}/samples_{pdb_id}_"
                    + row["variant"]
                    + "/*.npz"
                )
            else:
                npz_wc = (
                    f"data/data_{dataset}/simulate_seq_fragments_{dataset}/samples_{pdb_id}_"
                    + row["variant"]
                    + "_mt"
                    + "/*.npz"
                )

            # Create a dataset of all residue environments for this fragment
            try:
                fragment_dataset_wt = ResidueEnvironmentsDataset(
                    sorted(glob.glob(npz_wc))[::stride], transformer=None
                )
            except ValueError as e:
                print(f"Variant {row['variant']} failed with error: {e}.")
                print(npz_wc)
                fragments_nlls_wt.append(np.array([np.nan, np.nan, np.nan]))
                fragments_nlls_mt.append(np.array([np.nan, np.nan, np.nan]))
                continue

            # Extract the ones with pdb_residue_number 0 (center)
            res_envs_center_wt = [
                res_env for res_env in fragment_dataset_wt if res_env.pdb_residue_number == 0
            ]
            # Create dataset and dataloder
            dataset_center = ResidueEnvironmentsDataset(
                res_envs_center_wt, transformer=ToTensor(DEVICE)
            )
            dataloader_center = DataLoader(
                dataset_center,
                batch_size=len(dataset_center),
                shuffle=False,
                drop_last=False,
                collate_fn=ToTensor(DEVICE).collate_cat,
            )
            # Get X tensor
            x_batch, _ = next(iter(dataloader_center))

            # Predict
            y_pred = (
                -torch.log(softmax(cavity_model_infer_net(x_batch), dim=-1) + EPS)
                .detach()
                .cpu()
                .numpy()
            )
            fragments_nlls_wt.append(y_pred[:, row["wt_idx"]])
            fragments_nlls_mt.append(y_pred[:, row["mt_idx"]])

        if is_wt:
            ddg_data_dict[dataset][f"fragment_nll_wt_given_wt_{j}"] = fragments_nlls_wt
            ddg_data_dict[dataset][f"fragment_nll_mt_given_wt_{j}"] = fragments_nlls_mt
        else:
            ddg_data_dict[dataset][f"fragment_nll_wt_given_mt_{j}"] = fragments_nlls_wt
            ddg_data_dict[dataset][f"fragment_nll_mt_given_mt_{j}"] = fragments_nlls_mt


def add_ddg_preds_with_unfolded_state(ddg_data_dict: Dict, dataset: str, n_models: int):
    """
    Calculates ddg prediction using the first only the WT fragments, then the MT fragments,
    and then both. Note that is still without any down stream model
    """
    for i in range(n_models):
        ddg_data_dict[dataset][f"ddg_pred_wt_phaistos_no_ds_{i}"] = ddg_data_dict[dataset].apply(
            lambda row: row[f"mt_nll_{i}"]
            - row[f"fragment_nll_mt_given_wt_{i}"].mean()
            - row[f"wt_nll_{i}"]
            + row[f"fragment_nll_wt_given_wt_{i}"].mean(),
            axis=1,
        )

        ddg_data_dict[dataset][f"ddg_pred_mt_phaistos_no_ds_{i}"] = ddg_data_dict[dataset].apply(
            lambda row: row[f"mt_nll_{i}"]
            - row[f"fragment_nll_mt_given_mt_{i}"].mean()
            - row[f"wt_nll_{i}"]
            + row[f"fragment_nll_wt_given_mt_{i}"].mean(),
            axis=1,
        )

        ddg_data_dict[dataset][f"ddg_pred_wt_and_mt_phaistos_no_ds_{i}"] = ddg_data_dict[
            dataset
        ].apply(
            lambda row: row[f"mt_nll_{i}"]
            - (
                row[f"fragment_nll_mt_given_wt_{i}"].mean()
                + row[f"fragment_nll_mt_given_mt_{i}"].mean()
            )
            / 2
            - row[f"wt_nll_{i}"]
            + (
                row[f"fragment_nll_wt_given_wt_{i}"].mean()
                + row[f"fragment_nll_wt_given_mt_{i}"].mean()
            )
            / 2,
            axis=1,
        )


def _get_residue_map(dataset: str, pdb_id: str, chain_id: str):
    conversion_filepath = (
        f"data/data_{dataset}/molecular_dynamics/"
        f"residue_number_mapping/{pdb_id}_{chain_id}_mapping.txt"
    )
    with open(conversion_filepath, "r") as f:
        residue_map = {line.strip().split()[0]: line.strip().split()[1] for line in f}

    return residue_map


def infer_molecular_dynamics_nlls(
    ddg_data_dict: Dict,
    dataset: str,
    DEVICE: str,
    EPS: float,
    cavity_model_infer_nets: List[CavityModel],
    *,
    stride: int = 400,
    is_npt_ensemble=False,
    use_residue_number_map=False,
):
    """
    Infer negative log likelihoods for MD simulations of folded state, i.e. only wild type.
    """
    md_datasets = {}
    for j, cavity_model_infer_net in enumerate(cavity_model_infer_nets):
        md_nlls_wt_given_wt = []
        md_nlls_mt_given_wt = []

        residue_numbering_maps = {}
        for i, row in ddg_data_dict[dataset].iterrows():
            print(i)
            if not is_npt_ensemble:
                md_pdbs = glob.glob(
                    f"data/data_{dataset}/molecular_dynamics/pdbs_parsed/{row['pdbid']}_*npz"
                )[::stride]
            else:
                md_pdbs = glob.glob(
                    f"data/data_{dataset}/molecular_dynamics_npt/pdbs_parsed/{row['pdbid'].lower()}_*npz"
                )[::stride]

            wt_restype = row["variant"][0]
            chain_id = row["chainid"]
            pdb_id = row["pdbid"]

            # Store whole parsed pdbs for quicker parsing
            if pdb_id in md_datasets:
                md_pdbs_dataset = md_datasets[pdb_id]
            else:
                md_pdbs_dataset = ResidueEnvironmentsDataset(sorted(md_pdbs), transformer=None)
                md_datasets[pdb_id] = md_pdbs_dataset

            # Store residue numbering conversion maps. The maps convert from MD residue numbering
            # to the orignal residue numbering
            if (pdb_id, chain_id) in residue_numbering_maps:
                residue_number_map = residue_numbering_maps[(pdb_id, chain_id)]
            else:
                residue_number_map = _get_residue_map(dataset, pdb_id, chain_id)
                residue_numbering_maps[(pdb_id, chain_id)] = residue_number_map

            # print(residue_number_map)

            # The MD pdb numbers need to be renumbered back to the original residue numbers
            md_pdb_residue_number = int(row["variant"][1:-1])
            if use_residue_number_map:
                pdb_residue_number = int(residue_number_map[str(md_pdb_residue_number)])
            else:
                pdb_residue_number = int(md_pdb_residue_number)

            # Extract resenvs that match restype and residue number
            extracted_resenvs_md = []
            for resenv in md_pdbs_dataset:
                if (
                    wt_restype == index_to_one(resenv.restype_index)
                    and pdb_residue_number == resenv.pdb_residue_number
                ):
                    extracted_resenvs_md.append(resenv)

            if len(extracted_resenvs_md) == 0:
                print(f"Did not find matchin resenvs for {row['variant']}")

            dataset_variant_md = ResidueEnvironmentsDataset(
                extracted_resenvs_md, transformer=ToTensor(DEVICE)
            )
            dataloader_variant_md = DataLoader(
                dataset_variant_md,
                batch_size=len(dataset_variant_md),
                shuffle=False,
                drop_last=False,
                collate_fn=ToTensor(DEVICE).collate_cat,
            )
            # Get X tensor
            x_batch, _ = next(iter(dataloader_variant_md))

            # Predict
            y_pred = (
                -torch.log(softmax(cavity_model_infer_net(x_batch), dim=-1) + EPS)
                .detach()
                .cpu()
                .numpy()
            )

            md_nlls_wt_given_wt.append(y_pred[:, row["wt_idx"]])
            md_nlls_mt_given_wt.append(y_pred[:, row["mt_idx"]])
        ddg_data_dict[dataset][f"wt_nll_md_{j}"] = md_nlls_wt_given_wt
        ddg_data_dict[dataset][f"mt_nll_md_{j}"] = md_nlls_mt_given_wt


def add_ddg_preds_with_md_simulations(ddg_data_dict: Dict, dataset: str, n_models: int):

    for i in range(n_models):
        # MD only
        ddg_data_dict[dataset][f"ddg_pred_md_no_ds_all_{i}"] = ddg_data_dict[dataset].apply(
            lambda row: row[f"mt_nll_md_{i}"].mean() - row[f"wt_nll_md_{i}"].mean(),
            axis=1,
        )
        # MD + PDB
        ddg_data_dict[dataset][f"ddg_pred_md_pdb_statistics_no_ds_all_{i}"] = ddg_data_dict[
            dataset
        ].apply(
            lambda row: row[f"mt_nll_md_{i}"].mean()
            - row[f"mt_nlf_{i}"]
            - row[f"wt_nll_md_{i}"].mean()
            + row[f"wt_nlf_{i}"],
            axis=1,
        )

        # MD + IDP statistics
        ddg_data_dict[dataset][f"ddg_pred_md_idp_statistics_no_ds_all_{i}"] = ddg_data_dict[
            dataset
        ].apply(
            lambda row: row[f"mt_nll_md_{i}"].mean()
            - row[f"mt_idp_nlf_{i}"]
            - row[f"wt_nll_md_{i}"].mean()
            + row[f"wt_idp_nlf_{i}"],
            axis=1,
        )

        # MD + phaistos WT and MT statistics
        ddg_data_dict[dataset][
            f"ddg_pred_md_phaistos_mt_and_wt_statistics_no_ds_{i}"
        ] = ddg_data_dict[dataset].apply(
            lambda row: row[f"mt_nll_md_{i}"].mean()
            - (
                row[f"fragment_nll_mt_given_wt_{i}"].mean()
                + row[f"fragment_nll_mt_given_mt_{i}"].mean()
            )
            / 2
            - row[f"wt_nll_md_{i}"].mean()
            + (
                row[f"fragment_nll_wt_given_wt_{i}"].mean()
                + row[f"fragment_nll_wt_given_mt_{i}"].mean()
            )
            / 2,
            axis=1,
        )

    # ddg_data_dict[dataset]["ddg_pred_md_pdb_statistics_no_ds_1"] = ddg_data_dict[dataset].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["mt_nlf"]
    #     - row["wt_nll_md"][0:1].mean()
    #     + row["wt_nlf"],
    #     axis=1,
    # )
    # ddg_data_dict[dataset]["ddg_pred_md_pdb_statistics_no_ds_2"] = ddg_data_dict[dataset].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["mt_nlf"]
    #     - row["wt_nll_md"][0:2].mean()
    #     + row["wt_nlf"],
    #     axis=1,
    # )
    # ddg_data_dict[dataset]["ddg_pred_md_pdb_statistics_no_ds_5"] = ddg_data_dict[dataset].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["mt_nlf"]
    #     - row["wt_nll_md"][0:5].mean()
    #     + row["wt_nlf"],
    #     axis=1,
    # )
    # ddg_data_dict[dataset]["ddg_pred_md_pdb_statistics_no_ds_20"] = ddg_data_dict[dataset].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["mt_nlf"]
    #     - row["wt_nll_md"][0:20].mean()
    #     + row["wt_nlf"],
    #     axis=1,
    # )
    # ddg_data_dict[dataset]["ddg_pred_md_pdb_statistics_no_ds_all"] = ddg_data_dict[dataset].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["mt_nlf"]
    #     - row["wt_nll_md"].mean()
    #     + row["wt_nlf"],
    #     axis=1,
    # )
    # ddg_data_dict[dataset]["ddg_pred_md_pdb_statistics_no_ds_200"] = ddg_data_dict[dataset].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["mt_nlf"]
    #     - row["wt_nll_md"][0:200].mean()
    #     + row["wt_nlf"],
    #     axis=1,
    # )

    # # MD + idp statistics
    # ddg_data_dict[dataset]["ddg_pred_md_idp_statistics_no_ds"] = ddg_data_dict[dataset].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["mt_idp_nlf"]
    #     - row["wt_nll_md"].mean()
    #     + row["wt_idp_nlf"],
    #     axis=1,
    # )

    # # MD + phaistos WT statistics
    # ddg_data_dict[dataset]["ddg_pred_md_phaistos_wt_statistics_no_ds"] = ddg_data_dict[
    #     dataset
    # ].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["fragment_nll_mt_given_wt"].mean()
    #     - row["wt_nll_md"].mean()
    #     + row["fragment_nll_wt_given_wt"].mean(),
    #     axis=1,
    # )

    # # MD + phaistos MT statistics
    # ddg_data_dict[dataset]["ddg_pred_md_phaistos_mt_statistics_no_ds"] = ddg_data_dict[
    #     dataset
    # ].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - row["fragment_nll_mt_given_mt"].mean()
    #     - row["wt_nll_md"].mean()
    #     + row["fragment_nll_wt_given_mt"].mean(),
    #     axis=1,
    # )

    # # MD + phaistos WT and MT statistics
    # ddg_data_dict[dataset]["ddg_pred_md_phaistos_mt_and_wt_statistics_no_ds"] = ddg_data_dict[
    #     dataset
    # ].apply(
    #     lambda row: row["mt_nll_md"].mean()
    #     - (row["fragment_nll_mt_given_wt"].mean() + row["fragment_nll_mt_given_mt"].mean()) / 2
    #     - row["wt_nll_md"].mean()
    #     + (row["fragment_nll_wt_given_wt"].mean() + row["fragment_nll_wt_given_mt"].mean()) / 2,
    #     axis=1,
    # )


def get_predictions_both_structures(ddg_data_dict: Dict, n_models: int):
    # Rename columns so they specify if it is the direct or inverse direction
    symmetric_direct_df = ddg_data_dict["symmetric_direct"].copy(deep=True)
    symmetric_direct_df.columns = [
        name + "_dir" if "_dir" not in name else name for name in symmetric_direct_df.columns
    ]
    symmetric_inverse_df = ddg_data_dict["symmetric_inverse"].copy(deep=True)
    symmetric_inverse_df.columns = [
        name + "_inv" if "_inv" not in name else name for name in symmetric_inverse_df.columns
    ]

    # Inner merge both dataframes
    ddg_data_dict["symmetric_both"] = pd.merge(
        symmetric_direct_df,
        symmetric_inverse_df,
        how="inner",
        left_on="merge_column_dir",
        right_on="merge_column_inv",
    )

    for i in range(n_models):
        ddg_data_dict["symmetric_both"][f"ddg_pred_no_ds_both_dir_{i}"] = ddg_data_dict[
            "symmetric_both"
        ].apply(lambda row: row[f"wt_nll_{i}_inv"] - row[f"wt_nll_{i}_dir"], axis=1)

        ddg_data_dict["symmetric_both"][f"ddg_pred_no_ds_both_inv_{i}"] = ddg_data_dict[
            "symmetric_both"
        ].apply(lambda row: row[f"wt_nll_{i}_dir"] - row[f"wt_nll_{i}_inv"], axis=1)

    # ddg_data_dict["symmetric_both"]["ddg_pred_no_ds_both_dir"] = ddg_data_dict[
    #     "symmetric_both"
    # ].apply(lambda row: 0.5 * (row["ddg_pred_no_ds_dir"] - row["ddg_pred_no_ds_inv"]), axis=1)

    # ddg_data_dict["symmetric_both"]["ddg_pred_no_ds_both_inv"] = ddg_data_dict[
    #     "symmetric_both"
    # ].apply(lambda row: 0.5 * (row["ddg_pred_no_ds_inv"] - row["ddg_pred_no_ds_dir"]), axis=1)


def output_sequence_fragments_to_csv(ddg_data_dict: Dict):
    # Add flanking sequence fragments for protein g
    raw_pdbs = glob.glob("data/data_protein_g/pdbs_raw/*.pdb")
    for raw_pdb in raw_pdbs:
        add_flanking_seq_fragments(
            ddg_data_dict,
            "protein_g",
            raw_pdb,
        )

    # Add flanking sequence fragments for guerois
    raw_pdbs = glob.glob("data/data_guerois/pdbs_raw/*.pdb")
    for raw_pdb in raw_pdbs:
        add_flanking_seq_fragments(
            ddg_data_dict,
            "guerois",
            raw_pdb,
        )

    # Add flanking sequence fragments for dms
    raw_pdbs = glob.glob("data/data_dms/pdbs_raw/*.pdb")
    for raw_pdb in raw_pdbs:
        add_flanking_seq_fragments(
            ddg_data_dict,
            "dms",
            raw_pdb,
        )

    # Output CSVs (So Wouter can simulate them)
    ddg_data_dict["protein_g"][
        ["pdbid", "variant", "left_flank", "right_flank", "wt_restype", "mt_restype"]
    ].to_csv("data/data_protein_g/sequence_flanks_protein_g.csv")
    ddg_data_dict["dms"][
        ["pdbid", "variant", "left_flank", "right_flank", "wt_restype", "mt_restype"]
    ].to_csv("data/data_dms/sequence_flanks_dms.csv")
    ddg_data_dict["guerois"][
        ["pdbid", "variant", "left_flank", "right_flank", "wt_restype", "mt_restype"]
    ].to_csv("data/data_guerois/sequence_flanks_guerois.csv")
