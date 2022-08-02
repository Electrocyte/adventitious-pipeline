#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 10:15:27 2021

@author: mangi
"""

import pandas as pd
import os


def process(lines=None):
    ks = ["name", "sequence", "optional", "quality"]
    return {k: v for k, v in zip(ks, lines)}


def generate_fastq_df(fastq_file: str) -> pd.DataFrame:
    rows = []
    n = 4
    with open(fastq_file, "r") as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process(lines)
                rows.append(record)
                lines = []
    df = pd.DataFrame.from_records(rows)
    df = df.rename(
        columns={
            "name": "readID",
            "sequence": "Sequence",
            "optional": "Plus",
            "quality": "Quality",
        }
    )

    df.Sequence = df.Sequence.str.replace("U", "T")
    return df


def fastq_formatter(directory: str, df: pd.DataFrame, sample: str) -> str:
    lines = []
    for _index, row in df.iterrows():
        lines.append(row["readID"])
        lines.append(row["Sequence"])
        lines.append(row["Plus"])
        lines.append(row["Quality"])
    os.makedirs(f"{directory}", exist_ok=True)
    out_file = f"{directory}/U2T_{sample}.fastq"
    if not os.path.isfile(out_file):
        print(f"Saving file as: {out_file}")
        with open(out_file, "w") as f:
            f.write("\n".join(lines))
    return out_file


def main(filename: str, in_dir: str, sample: str) -> None:
    fq_df = generate_fastq_df(filename)
    out_file = fastq_formatter(in_dir, fq_df, sample)


if __name__ == "__main__":

    sample = "20210323_RNA_PAO1unfTCPA_38700000CFU_1st_2h"
    in_dir = (
        f"/home/james/SequencingData/Direct_RNA/analysis/sample_data/{sample}/trimmed/"
    )
    filename = f"{in_dir}/trimmed_{sample}.fastq"
    fq_df = generate_fastq_df(filename)
    fastq_formatter(in_dir, fq_df, sample)
