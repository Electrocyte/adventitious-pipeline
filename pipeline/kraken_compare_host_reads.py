#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 11:41:06 2021

@author: mangi
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import os


def clean_df(df: pd.DataFrame, status: str, species: str) -> pd.DataFrame:
    df = df.loc[df["taxName"] == species]
    df = df.drop(
        labels=[
            "reads",
            "kmers",
            "dup",
            "taxID",
            "rank",
            "date",
            "NA",
            "strain",
            "concentration_CFU",
            "batch",
            "duration_h",
        ],
        axis=1,
    )
    new_names = [
        (i, i + f"_{status}") for i in df.columns if i not in ["sample", "taxName"]
    ]
    df.rename(columns=dict(new_names), inplace=True)
    return df


def draw_host_reads_as_percent(
    df: pd.DataFrame, title: str, directory: str, save: bool
) -> None:
    f, ax = plt.subplots(figsize=(15, 15))

    ax = sns.barplot(
        x="Host reads as % of total reads", y="Samples", data=df, palette="bright"
    )

    # ax.set_xscale('log')
    ax.set_xlabel("Host reads as % of total reads", size=40)
    ax.set_ylabel("Sample", size=40)

    plt.title(title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    # plt.legend(loc=4, prop={'size': 30}, markerscale=4)

    for y, value in enumerate(df["Host reads as % of total reads"]):
        str_value = f'{"{:.0f}".format(value)}'
        ax.text(value, y, str_value, size=20)
    print(f"Saving to: {directory}/{title}.png")
    if save:
        plt.savefig(f"{directory}/{title}.png", dpi=300, bbox_inches="tight")
    plt.show()


def draw_reads_from_kraken(
    df: pd.DataFrame, title: str, directory: str, save: bool
) -> None:
    f, ax = plt.subplots(figsize=(15, 15))

    ax = sns.barplot(
        x="Additional host reads found by krakenuniq (%)",
        y="Samples",
        data=df,
        palette="bright",
    )

    # ax.set_xscale('log')
    ax.set_xlabel("Additional host reads found by krakenuniq (%)", size=40)
    ax.set_ylabel("Sample", size=40)

    plt.title(title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    # plt.legend(loc=4, prop={'size': 30}, markerscale=4)

    for y, value in enumerate(df["Additional host reads found by krakenuniq (%)"]):
        str_value = f'{"{:.2f}".format(value)}'
        ax.text(value, y, str_value, size=20)
    print(f"Saving to: {directory}/{title}.png")
    if save:
        plt.savefig(f"{directory}/{title}.png", dpi=300, bbox_inches="tight")
    plt.show()


def draw_reads_v_host(df: pd.DataFrame, title: str, directory: str, save: bool) -> None:
    f, ax = plt.subplots(figsize=(15, 15))

    ax = sns.scatterplot(
        x="Host reads as % of total reads",
        y="total_reads",
        data=df,
        palette="bright",
        hue="batch",
        s=200,
    )

    ax.set_xlabel("Host reads as % of total reads", size=40)
    ax.set_ylabel("Total Reads", size=40)

    plt.title(title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    plt.legend(loc=2, prop={"size": 20}, markerscale=4)

    for y, value in enumerate(df["Host reads as % of total reads"]):
        sample_name = df.strain[y]
        yy = df.total_reads[y]
        # str_value = f'{"{:.0f}".format(value)}%, {sample_name}'
        str_value = f"{sample_name}"
        ax.text(value, yy, str_value, size=20)
    print(f"Saving to: {directory}/{title}.png")
    if save:
        plt.savefig(f"{directory}/{title}.png", dpi=300, bbox_inches="tight")
    plt.show()


def draw_reads_v_host_reads(
    df: pd.DataFrame, title: str, directory: str, save: bool
) -> None:
    f, ax = plt.subplots(figsize=(15, 15))

    ax = sns.scatterplot(
        x="kraken_host_reads",
        y="total_reads",
        data=df,
        palette="bright",
        hue="batch",
        s=200,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Total host reads", size=40)
    ax.set_ylabel("Total Reads", size=40)

    ax.set_xlim([100, max(df["kraken_host_reads"])])
    ax.set_ylim([100, max(df["total_reads"])])

    plt.title(title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    plt.legend(loc=4, prop={"size": 20}, markerscale=4)
    print(f"Saving to: {directory}/{title}.png")
    if save:
        plt.savefig(f"{directory}/{title}.png", dpi=300, bbox_inches="tight")
    plt.show()


def draw_reads_v_host_reads_log(
    df: pd.DataFrame, title: str, directory: str, save: bool
) -> None:
    f, ax = plt.subplots(figsize=(15, 15))

    ax = sns.scatterplot(
        x="logcentrifuge_host_reads",
        y="logkraken_host_reads",
        data=df,
        palette="bright",
        hue="batch",
        s=200,
    )

    ax.set_xlabel("Centrifuge host reads (log(x))", size=40)
    ax.set_ylabel("Kraken host reads (log(y))", size=40)

    ax.set_xlim([-1, max(df["logcentrifuge_host_reads"])])
    ax.set_ylim([-1, max(df["logkraken_host_reads"])])

    plt.title(title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    plt.legend(loc=4, prop={"size": 20}, markerscale=4)
    print(f"Saving to: {directory}/{title}.png")
    if save:
        plt.savefig(f"{directory}/{title}.png", dpi=300, bbox_inches="tight")
    plt.show()


def draw_coverage_from_kraken(
    df: pd.DataFrame, title: str, directory: str, save: bool
) -> None:
    f, ax = plt.subplots(figsize=(15, 15))

    ax = sns.barplot(
        x="coverage per kraken (%)", y="Samples", data=df, palette="bright"
    )

    # ax.set_xscale('log')
    ax.set_xlabel("Coverage for host genome as per kraken (%)", size=40)
    ax.set_ylabel("Sample", size=40)

    plt.title(title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    # plt.legend(loc=4, prop={'size': 30}, markerscale=4)

    for y, value in enumerate(df["coverage per kraken (%)"]):
        str_value = f'{"{:.2f}".format(value)}'
        ax.text(value, y, str_value, size=20)
    print(f"Saving to: {directory}/{title}.png")
    if save:
        plt.savefig(f"{directory}/{title}.png", dpi=300, bbox_inches="tight")
    plt.show()


def draw_coverage_from_kraken_minimap2(
    df: pd.DataFrame, title: str, directory: str, save: bool
) -> None:
    f, ax = plt.subplots(figsize=(15, 15))

    ax = sns.scatterplot(
        x="coverage per kraken (%)",
        y="percent_covered_mean",
        data=df,
        palette="bright",
        hue="batch",
        s=200,
    )

    # ax.set_xscale('log')
    ax.set_xlabel("Kraken coverage (%)", size=40)
    ax.set_ylabel("Minimap2 coverage (%)", size=40)

    ax.set_xlim([-1, 105])
    ax.set_ylim([-1, 105])

    plt.title(title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    plt.legend(loc=1, prop={"size": 20}, markerscale=4)

    print(f"Saving to: {directory}/{title}.png")
    if save:
        plt.savefig(f"{directory}/{title}.png", dpi=300, bbox_inches="tight")
    plt.show()


def log_not_zero_hs(row):
    if row["hsbn_host_reads"] > 0:
        return np.log10(row["hsbn_host_reads"])
    if row["hsbn_host_reads"] == 0:
        return row["hsbn_host_reads"]


def log_not_zero_ce(row):
    if row["centrifuge_host_reads"] > 0:
        return np.log10(row["centrifuge_host_reads"])
    if row["centrifuge_host_reads"] == 0:
        return row["centrifuge_host_reads"]


def log_not_zero_kr(row):
    if row["kraken_host_reads"] > 0:
        return np.log10(row["kraken_host_reads"])
    if row["kraken_host_reads"] == 0:
        return row["kraken_host_reads"]


def run_host_kraken_charting(
    directory: str,
    kraken_cov_folder: str,
    AA_analysis: str,
    species: str,
    save: bool,
    minimap_host: str,
    host_name_analysis: str,
    samples: list,
):

    analysis_directory = f"{directory}/analysis"
    trim = f"{analysis_directory}/{kraken_cov_folder}/aggregate_kraken_coverage.csv"

    no_host = f"{analysis_directory}/{kraken_cov_folder}/aggregate_host_removed_kraken_coverage.csv"

    # will likely cause issues with new samples, might need to only use the large file
    metagenome_reads = f"{analysis_directory}/{AA_analysis}/metagenome_reads.csv"
    print(analysis_directory, trim, no_host, metagenome_reads)
    trim_df = pd.read_csv(trim, index_col=False)
    if "Unnamed: 0" in trim_df.columns:
        trim_df = trim_df.drop(labels=["Unnamed: 0"], axis=1)
    t_status = "host"

    nh_df = pd.read_csv(no_host)
    nh_status = "no_host"

    trim_df = clean_df(trim_df, t_status, species)
    nh_df = clean_df(nh_df, nh_status, species)

    # duplicated names causes error:
    # raise ValueError("Columns must be same length as key")
    # ValueError: Columns must be same length as key
    try:
        trim_df[
            ["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]
        ] = trim_df["sample"].str.split("_", 5, expand=True)
        nh_df[
            ["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]
        ] = nh_df["sample"].str.split("_", 5, expand=True)
    except ValueError:
        print("Columns exist! Skipping.")

    df = pd.merge(trim_df, nh_df, on=["sample", "taxName"], how="outer")
    df = df.fillna(0)

    metagenome_reads_glob = glob.glob(metagenome_reads)[0]

    mg_df = pd.read_csv(metagenome_reads_glob)
    mg_df = mg_df.loc[mg_df["sample"].isin(samples)]
    mg_df.reset_index(inplace=True, drop=True)
    mg_df = mg_df.drop(
        labels=["hs_classified", "cent_classified", "hs_unc", "ce_unc"], axis=1
    )
    all_df = pd.merge(df, mg_df, on=["sample"], how="outer")

    all_df["Total_non-host_reads"] = all_df["hs_all"] - all_df[f"taxReads_{t_status}"]

    all_df["Total host reads removed by krakenuniq (%)"] = (
        all_df[f"taxReads_{t_status}"] / all_df["hs_all"] * 100
    )
    all_df["Additional host reads found by krakenuniq (%)"] = (
        all_df[f"taxReads_{nh_status}"] / all_df["hs_all"] * 100
    )
    all_df["Host reads as % of total reads"] = (
        all_df[f"taxReads_{t_status}"] / all_df["hs_all"] * 100
    )
    all_df["cov_host"] = all_df["cov_host"] * 100

    all_df.rename(
        columns={
            "taxReads_no_host": "reads unique to krakenuniq",
            "taxReads_host": "kraken_host_reads",
            "hs_all": "total_reads",
            "cov_host": "coverage per kraken (%)",
        },
        inplace=True,
    )

    columns2save = [
        "sample",
        "taxName",
        "coverage per kraken (%)",
        "hsbn_host_reads",
        "centrifuge_host_reads",
        "kraken_host_reads",
        "Total_non-host_reads",
        "total_reads",
        "Host reads as % of total reads",
        "reads unique to krakenuniq",
        "Additional host reads found by krakenuniq (%)",
        "Total host reads removed by krakenuniq (%)",
    ]

    summary_df = all_df[columns2save]
    summary_df.to_csv(
        f"{analysis_directory}/{kraken_cov_folder}/Removing_host_reads.csv", index=False
    )
    print(
        f"Saving summary df to: {analysis_directory}/{kraken_cov_folder}/Removing_host_reads.csv"
    )

    out_dir = f"{analysis_directory}/{kraken_cov_folder}/"
    out_dir_batch = f"{analysis_directory}/{kraken_cov_folder}/charts"
    os.makedirs(out_dir_batch, exist_ok=True)

    summary_df[
        ["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]
    ] = summary_df["sample"].str.split("_", 5, expand=True)
    summary_df["Samples"] = summary_df["strain"] + "-" + summary_df["batch"].astype(str)

    for batch in list(trim_df["batch"].unique()):
        batch_df = summary_df.loc[summary_df["batch"] == batch]
        if not batch_df.empty:
            try:
                draw_host_reads_as_percent(
                    batch_df,
                    f"All assessed host reads as a percentage of total reads{batch}",
                    out_dir_batch,
                    save,
                )
                draw_reads_from_kraken(
                    batch_df,
                    f"Host reads not found by centrifuge or HS-BLASTn{batch}",
                    out_dir_batch,
                    save,
                )

                batch_df["loghsbn_host_reads"] = batch_df.apply(log_not_zero_hs, axis=1)
                batch_df["logcentrifuge_host_reads"] = batch_df.apply(
                    log_not_zero_ce, axis=1
                )
                batch_df["logkraken_host_reads"] = batch_df.apply(
                    log_not_zero_kr, axis=1
                )

                draw_reads_v_host(
                    batch_df,
                    f"Reads generated compared to percent of identified host reads{batch}",
                    out_dir_batch,
                    save,
                )
                draw_reads_v_host_reads(
                    batch_df,
                    f"Reads generated compared to identified host reads{batch}",
                    out_dir_batch,
                    save,
                )
                draw_reads_v_host_reads_log(
                    batch_df,
                    f"Krakenuniq host reads vs centrifuge host reads{batch}",
                    out_dir_batch,
                    save,
                )
                draw_coverage_from_kraken(
                    batch_df,
                    f"Krakenuniq coverage for host reads{batch}",
                    out_dir_batch,
                    save,
                )
            except KeyError:
                print("df does not contain useful data")

    coverage_minimap2 = f"{analysis_directory}/{host_name_analysis}/{minimap_host}_host_centrifuge_to_coverage_agg.csv"
    cms_df = pd.read_csv(coverage_minimap2)
    cms_df["date"] = cms_df["date"].astype(str)
    cms_df["concentration_CFU"] = cms_df["concentration_CFU"].astype(str) + "CFU"
    cms_df["batch"] = cms_df["batch"].astype(str)
    cms_df["duration_h"] = cms_df["duration_h"].astype(str)

    all_cov = pd.merge(
        summary_df,
        cms_df,
        how="left",
        left_on=["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"],
        right_on=["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"],
    )

    draw_coverage_from_kraken_minimap2(
        all_cov,
        "Krakenuniq coverage vs minimap2 alignment & samtools for host reads",
        out_dir,
        save,
    )


if __name__ == "__main__":
    directory = "E:/SequencingData/Viral_CHO/"
    AA_analysis = "coverage_summary_cviral_virus_virus_no_host_OCS"
    species = "Cricetulus griseus"
    kraken_cov_folder = f"kraken_coverage_host_{species.replace(' ', '_')}"
    save = True
    minimap_host = "chinese_hamster_chromosomes"
    host_name_analysis = f"coverage_summary_{minimap_host}_host_original"
    samples = [
        "20211018_ssDNA_MVM1-100CHOK1-Concentrate_50000CFU_20_36",
        "20211018_ssDNA_MVM1-100CHOK1-Supernatant_50000CFU_20_36",
        "20211018_ssDNA_MVM1-10CHOK1-Concentrate_500000CFU_20_36",
        "20211018_ssDNA_MVM1-10CHOK1-Filtrate_500000CFU_20_36",
        "20211018_ssDNA_MVM1-10CHOK1-Supernatant_500000CFU_20_36",
    ]
    run_host_kraken_charting(
        directory, kraken_cov_folder, AA_analysis, species, save, samples
    )
