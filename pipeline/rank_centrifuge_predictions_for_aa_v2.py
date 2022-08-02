#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 13:58:24 2021

@author: mangi
"""

import glob
import pandas as pd


# called in mp_metagenomic_assessment_v4.py
def generate_centrifuge_ranks(
    adventitious_agent: str, centrifuge_file: str, directory: str, kingdom: str
) -> pd.DataFrame:
    title_aa = adventitious_agent.title()
    lower_aa = adventitious_agent.lower()

    input_file = glob.glob(f"{centrifuge_file}")
    print(input_file, centrifuge_file)
    df = pd.read_csv(input_file[0])
    df[["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]] = df[
        "sample"
    ].str.split("_", 5, expand=True)
    df["concentration_CFU"] = df["concentration_CFU"].astype(str) + "CFU"
    df["date"] = df["date"].astype(str)
    df["batch"] = df["batch"].astype(str)

    df_sorted = df.sort_values(
        [
            "date",
            "strain",
            "numUniqueReads",
        ],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    aa_vals = 0
    ranks = []
    for sample in df_sorted["sample"].unique():
        df_sorted_intermediate = df_sorted.loc[df_sorted["sample"] == sample]
        date, na, strain, cfu, batch, time = sample.split("_")
        cfu = cfu.replace("CFU", "")
        isolate_df = df_sorted_intermediate.loc[
            (df_sorted_intermediate["date"] == date)
            & (df_sorted_intermediate["batch"] == batch)
            & (df_sorted_intermediate["strain"] == strain)
        ]
        isolate_df.reset_index(inplace=True, drop=True)
        isolate_df = isolate_df.drop_duplicates(
            subset=[
                "name",
                "date",
                "NA",
                "strain",
                "concentration_CFU",
                "batch",
                "duration_h",
            ],
            keep="first",
        )
        isolate_df.reset_index(inplace=True, drop=True)

        aa_local_indices = isolate_df[
            isolate_df["name"].str.contains(
                f"{title_aa}|{lower_aa}|{adventitious_agent}"
            )
        ].index
        total_indices = len(isolate_df.index)

        aa_vals = isolate_df[
            isolate_df["name"].str.contains(
                f"{title_aa}|{lower_aa}|{adventitious_agent}"
            )
        ].reset_index(drop=True)
        if len(aa_vals["numUniqueReads"]):
            aa_vals = aa_vals["numUniqueReads"][0]
        else:
            aa_vals = 0

        if len(aa_local_indices) > 0:
            rank = f"{str(int(aa_local_indices.values[0])+1)}/{str(total_indices)}"
            AA_rank = f"{str(int(aa_local_indices.values[0])+1)}"
        else:
            rank = "No target reads found"
            AA_rank = "No target reads found"
        rank_info = [
            date,
            na,
            strain,
            f"{cfu}CFU",
            batch,
            time,
            rank,
            AA_rank,
            total_indices,
            aa_vals,
        ]
        ranks.append(rank_info)

    df_ranks = pd.DataFrame(
        ranks,
        columns=[
            "date",
            "NA",
            "strain",
            "concentration_CFU",
            "batch",
            "duration_h",
            "rank",
            "AA",
            "no_ranks",
            "AA_reads",
        ],
    )
    df_ranks_sorted = df_ranks.sort_values(
        ["date", "rank"], ascending=[True, True]
    ).reset_index(drop=True)

    save_name = f"{directory}/describe_rank_analysis_for_{kingdom}_{adventitious_agent}_from_centrifuge_classification.csv"
    save_name = save_name.replace(" ", "_")
    print(save_name)
    df_ranks_sorted.to_csv(save_name, index=False)
    return df_ranks_sorted


if __name__ == "__main__":
    directory = "E:/SequencingData/Viral_CHO/analysis/"
    centrifuge_file = "centrifuge_virus_describe_meta.csv"
    centrifuge_file = "describe_OCS_centrifuge_meta_virus_no_host_all"
    filename = f"{directory}/{centrifuge_file}"
    kingdom = "virus"
    adventitious_agent = "Minute virus"
    generate_centrifuge_ranks(adventitious_agent, filename, directory)
