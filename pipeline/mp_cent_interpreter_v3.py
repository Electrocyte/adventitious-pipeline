#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 13:11:58 2021

@author: mangi
"""


import pandas as pd
import glob
import re



def extract_experimental_data(filename, all_data, number):
    is_data = all_data.findall(filename)

    date = is_data[0][0]
    NA = is_data[0][1]
    strain = is_data[0][2]
    CFU = is_data[0][3]
    CFU = number.findall(CFU)[0]
    batch = is_data[0][4]
    batch = number.findall(batch)[0]
    duration = is_data[0][5]
    duration = number.findall(duration)[0]

    df = pd.read_csv(filename, sep="\t")
    df["totalUniqReads"] = sum(df["numUniqueReads"])
    df["date"] = int(date)
    df["NA"] = NA
    df["strain"] = strain
    df["concentration_CFU"] = int(CFU)
    df["batch"] = int(batch)
    df["duration_h"] = int(duration)
    return df


def find_agg_describe(directory: str, kingdom: str):
    files = f"{directory}/analysis/sample_data/**/centrifuge/describe_{kingdom}_out.csv"
    print(files)
    glob_files = glob.glob(files)
    describe_dfs = []
    for file in glob_files:
        df = pd.read_csv(file)
        describe_dfs.append(df)
    cat_df = pd.concat(describe_dfs)
    return cat_df


def run_aggregate_centrifuge_files(
    directory: str, kingdom: str, no_host: str
) -> (str, str, pd.DataFrame, pd.DataFrame):
    number = re.compile(r"(^\d+)")
    all_data = re.compile(r"(\d+)_(\w+)_(\w.+)_(\w+)_(\w+)_(\w+)\/")
    files_to_find = f"{directory}/analysis/sample_data/**/centrifuge/2*{kingdom}{no_host}_centrifuge_report.tsv"
    print(f"Files to parse: {files_to_find}")
    centrifuge_files = glob.glob(files_to_find, recursive=True)

    cat_df = find_agg_describe(directory, kingdom)
    cat_df.to_csv(f"{directory}/analysis/describe_summary_{kingdom}.csv", index=False)

    experimental_dfs = []
    for file in centrifuge_files:
        if "\\" in file:
            file = file.replace("\\", "/")
        experimental_df = extract_experimental_data(file, all_data, number)
        if not experimental_df.empty:
            experimental_df = experimental_df[experimental_df.numUniqueReads > 0]
            experimental_dfs.append(experimental_df)

    meta_experimental_dfs = pd.concat(experimental_dfs)
    filename = f"centrifuge_analysis_for_{kingdom}{no_host}_all.csv"
    save1 = f"{directory}/analysis/{filename}"
    meta_experimental_dfs.to_csv(save1, index=False)

    meta_experimental_dfs["sample"] = meta_experimental_dfs[
        ["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]
    ].apply(lambda row: "_".join(row.values.astype(str)), axis=1)
    to_drop = [
        "taxRank",
        "genomeSize",
        "date",
        "NA",
        "strain",
        "concentration_CFU",
        "batch",
        "duration_h",
    ]
    small_meta_experimental_dfs = meta_experimental_dfs.drop(to_drop, axis=1)

    merge_save = f"{directory}/analysis/centrifuge_{kingdom}_describe_meta.csv"
    print(f"Saving to: {merge_save}")
    merge_describe_meta = pd.merge(
        small_meta_experimental_dfs,
        cat_df,
        how="left",
        left_on=["sample", "taxID"],
        right_on=["sample", "taxID"],
    )
    merge_describe_meta.to_csv(merge_save, index=False)

    return merge_save, meta_experimental_dfs


def main(
    directory: str,
    kingdom: str,
    no_host: str,
    expt_species: dict,
) -> str:

    merge_save, meta_experimental_dfs = run_aggregate_centrifuge_files(
        directory, kingdom, no_host
    )

    return merge_save


if __name__ == "__main__":
    directory = "/mnt/usersData/DNA/"
    kingdom = "bacteria"
    no_host = "_no_host"
    expt_species = ""
    build_new_model = False
    override_path = ""
    main(directory, kingdom, no_host, expt_species, build_new_model, override_path)

# time ./pipeline/mp_cent_interpreter_v3.py
