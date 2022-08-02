#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 13:44:54 2021

@author: mangi
"""

import glob
import pandas as pd
import re
import numpy as np
from pathlib import Path

# from pipeline.ML_tool import hs_ML_check_v5, hs_use_ML_model_v5
from pipeline import rank_BLAST_predictions_for_aa_v2


def count_reads(row):
    counter = 0
    if "numReads" in row and "numUniqueReads" in row:
        if row["numReads"] > 0:
            counter += 1
        if row["numUniqueReads"] > 0:
            counter += 1
        if "OCS" in row.index:
            counter += 2

    if "length_count" in row:
        if row["length_count"] > 0:
            counter += 2

    if "max_bitscore" in row:
        if row["max_bitscore"] > 200:
            counter += 1
    if "bitscore_max" in row:
        if row["bitscore_max"] > 200:
            counter += 1
    if "OCS" in row.index:
        counter += 2
    return counter


def cent_get(directory: str) -> pd.DataFrame:
    cent_glob = glob.glob(directory, recursive=True)
    df = pd.read_csv(cent_glob[0], sep=",")
    df["analysis"] = "centrifuge"
    df["c_score"] = df.apply(count_reads, axis=1)
    df[["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]] = df[
        "sample"
    ].str.split("_", 5, expand=True)
    df["concentration_CFU"] = df["concentration_CFU"].astype(np.int64)
    df["batch"] = df["batch"].astype(np.int64)
    df["duration_h"] = df["duration_h"].astype(np.int64)
    df["date"] = df["date"].astype(np.int64)
    return df


def BLASTn_get_one(BLASTn_directory: str) -> pd.DataFrame:
    number = re.compile(r"(^\d+)")
    all_data = re.compile(
        r"(\d+)_(\w+)_(\w.+)_(\d+[PC]FU)_(\w+)_(\d+)(h{0,1})(fast){0,1}\/"
    )
    BLASTn_globs = glob.glob(BLASTn_directory)
    dfs = []
    for BLASTn_glob in BLASTn_globs:
        df = pd.read_csv(BLASTn_glob, sep=",")

        if "\\" in BLASTn_glob:
            BLASTn_glob = BLASTn_glob.replace("\\", "/")

        if not df.empty:
            is_data = all_data.findall(BLASTn_glob)
            date = is_data[0][0]
            NA = is_data[0][1]
            strain = is_data[0][2]
            CFU = is_data[0][3]
            CFU = number.findall(CFU)[0]
            batch = is_data[0][4]
            batch = number.findall(batch)[0]
            duration = is_data[0][5]
            duration = number.findall(duration)[0]

            df["date"] = int(date)
            df["NA"] = NA
            df["strain"] = strain
            df["concentration_CFU"] = int(CFU)
            df["batch"] = int(batch)
            df["duration_h"] = int(duration)
            dfs.append(df)

    BLASTn_dfs = pd.concat(dfs)
    BLASTn_dfs["species"] = BLASTn_dfs["species"].astype(str)
    BLASTn_dfs["name"] = BLASTn_dfs["genus"] + " " + BLASTn_dfs["species"]
    BLASTn_dfs = BLASTn_dfs.drop(["sample", "genus", "species"], axis=1)
    BLASTn_dfs["analysis"] = "hsblastN"
    BLASTn_dfs["b_score"] = BLASTn_dfs.apply(count_reads, axis=1)
    return BLASTn_dfs


def clean_up_out(combined: pd.DataFrame) -> pd.DataFrame:
    combined["c_score"] = combined["c_score"].fillna(0)
    combined["b_score"] = combined["b_score"].fillna(0)
    combined["numUniqueReads"] = combined["numUniqueReads"].fillna(0)
    combined["totalUniqReads"] = combined["totalUniqReads"].fillna(0)
    combined["prediction_score"] = combined["b_score"] + combined["c_score"]
    combined = combined.drop(columns=["b_score", "c_score"])
    return combined


def run_merge_bn_cent(
    directory: str,
    no_host: str,
    BLASTn_name: str,
    cent_df: pd.DataFrame,
    hs_species: dict,
) -> (str, str, pd.DataFrame, pd.DataFrame, pd.DataFrame, list):
    BLASTn_directory = (
        f"{directory}/analysis/sample_data/2*/blastN/{BLASTn_name}/describe_top_ten.csv"
    )

    # find all blastn outputs
    BLASTn_df = BLASTn_get_one(BLASTn_directory)
    BLASTn_df["hs_database"] = BLASTn_name
    groupby = cent_df.merge(
        BLASTn_df,
        how="outer",
        on=[
            "date",
            "NA",
            "strain",
            "concentration_CFU",
            "batch",
            "duration_h",
            "name",
            "analysis",
        ],
    )
    combined_shorts = clean_up_out(groupby)
    priority_list = combined_shorts.loc[combined_shorts["prediction_score"] > 2]

    # original version
    save = (
        f"{directory}/analysis/describe_{BLASTn_name}{no_host}_combined_cent_blastn.csv"
    )
    priority_save = f"{directory}/analysis/describe_{BLASTn_name}{no_host}_priority_combined_cent_blastn.csv"
    BLASTn_agg = f"{directory}/analysis/describe_{BLASTn_name}{no_host}_agg.csv"

    # # handle a single species of adventitious agent
    BLAST_spp_ranked = []
    combined_shorts.to_csv(save, index=False)
    priority_list.to_csv(priority_save, index=False)
    BLASTn_df.to_csv(BLASTn_agg, index=False)

    ########### FOR ML ###########
    analysis_directory = f"{directory}/analysis/"
    HSBLASTn_all_dir = f"{directory}/analysis/sample_data/2*/blastN/{BLASTn_name}/describe_all_predictions.csv"
    HSBLASTn_df_all = BLASTn_get_one(HSBLASTn_all_dir)
    filename = f"describe_{BLASTn_name}{no_host}_all_agg.csv"
    HSBLASTn_agg = f"{analysis_directory}{filename}"
    print(f"Saving BLAST aggregate file to: {HSBLASTn_agg}")
    HSBLASTn_df_all.to_csv(HSBLASTn_agg, index=False)
    ########### FOR ML ###########

    # handle multiple different adventitious agent search species
    if isinstance(hs_species, dict):
        for AA, acronym in hs_species.items():
            rank_BLAST_predictions_for_aa_v2.generate_BLAST_ranks(
                AA, HSBLASTn_agg, analysis_directory, BLASTn_name
            )
            BLAST_spp_ranked.append(AA)

    return (
        HSBLASTn_agg,
        filename,
        HSBLASTn_df_all,
        priority_list,
        combined_shorts,
        BLAST_spp_ranked,
    )


def main(
    directory: str,
    barcoded: bool,
    kingdom: str,
    BLASTn_name: str,
    meta_data_name: str,
    no_host: str,
    cent_df_str: str,
    hs_species: dict,
) -> (str, str, list):
    print(f"Centrifuge summary file: {cent_df_str}")
    cent_df = cent_get(f"{cent_df_str}")

    (
        HSBLASTn_agg,
        filename,
        HSBLASTn_df_all,
        priority_list,
        combined_shorts,
        BLAST_spp_ranked,
    ) = run_merge_bn_cent(directory, no_host, BLASTn_name, cent_df, hs_species)

    if "priority" in meta_data_name:
        return (
            f"describe_{BLASTn_name}{no_host}_priority_combined_cent_blastn.csv",
            "priority",
            BLAST_spp_ranked,
        )
    if "priority" not in meta_data_name:
        return (
            f"describe_{BLASTn_name}{no_host}_combined_cent_blastn.csv",
            "nofilter",
            BLAST_spp_ranked,
        )


if __name__ == "__main__":
    # directory = "E:/SequencingData/Viral_CHO/"
    directory = "/mnt/usersData/DNA/"
    barcoded = False
    meta_data_name = "priority"
    kingdom = "virus"
    BLASTn_name = "cviral"
    no_host = "_no_host"
    # cent_df_str = f"{directory}/analysis/centrifuge_bacteria_describe_meta.csv"
    # cent_df_str = f"{directory}/analysis/centrifuge_fungus_describe_meta.csv"
    cent_df_str = f"{directory}/analysis/centrifuge_virus_describe_meta.csv"
    build_new_model = False
    # hs_species = {"Minute virus":"MVM"}
    hs_species = {}

    main(
        directory,
        barcoded,
        kingdom,
        BLASTn_name,
        meta_data_name,
        no_host,
        cent_df_str,
        hs_species,

    )
