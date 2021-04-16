import pandas as pd
import sys
from Bio import SeqIO
import glob
import re


def trim_left_flank(left_flank: str):
    if len(left_flank) <= 5:
        padding = "-" * (5 - len(left_flank))
        return padding + left_flank
    else:
        return left_flank[-5:]


def trim_right_flank(right_flank: str):
    if len(right_flank) <= 5:
        padding = "-" * (5 - len(right_flank))
        return right_flank + padding
    else:
        return right_flank[:5]


if __name__ == "__main__":
    csv_filename = sys.argv[1]
    mut_df = pd.read_csv(csv_filename)
    dataset = re.compile(r".+?/data_(.+?)/.+\.csv").match(csv_filename).groups()[0]
    pdb_info_dict = {}
    for pdb_filename in glob.glob(f"data/data_{dataset}/pdbs_raw/*.pdb"):
        for record in SeqIO.parse(pdb_filename, "pdb-seqres"):
            pdb_info_dict[record.id] = str(record.seq)

    left_flanks = []
    wt_residues = []
    mt_residues = []
    right_flanks = []
    for idx, row in mut_df.iterrows():
        record_id = f"{row['pdbid']}:{row['chainid']}"
        residue_number = int(row["variant"][1:-1])
        wt_residue = row["variant"][0]
        mt_residue = row["variant"][-1]

        seq = pdb_info_dict[record_id]
        resid = residue_number - 1
        left_flank = trim_left_flank(seq[:resid])
        right_flank = trim_right_flank(seq[resid + 1:])

        try:
            if wt_residue == seq[resid]:
                left_flanks.append(left_flank)
                wt_residues.append(wt_residue)
                right_flanks.append(right_flank)
                mt_residues.append(mt_residue)
            else:
                print(f"Mismatch problem with {record_id} {row['variant']}")
                left_flanks.append("")
                wt_residues.append(wt_residue)
                right_flanks.append("")
                mt_residues.append(mt_residue)
        except IndexError:
            print(f"Got index error for {record_id} {row['variant']}")

    mut_df["left_sequence_flank"] = left_flanks
    mut_df["wt_residue"] = wt_residues
    mut_df["right_sequence_flank"] = right_flanks
    mut_df["mt_residue"] = mt_residues

    mut_df.to_csv("temp.csv")
