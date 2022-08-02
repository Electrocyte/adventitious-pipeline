#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 16:18:22 2021

@author: mangi
"""

import pandas as pd
from pathlib import Path
import os
import glob
import re


def import_data(directory):
    # check if file is empty or not
    if os.stat(directory).st_size > 0:
        if "centrifuge" in directory:
            df = pd.read_csv(directory, sep="\t")
        if "blastN" in directory:
            df = pd.read_csv(
                directory,
                sep="\t",
                names=[
                    "qseqid",
                    "sseqid",
                    "pident",
                    "length",
                    "mismatches",
                    "gap_opens",
                    "q_start",
                    "q_end",
                    "s_start",
                    "s_end",
                    "evalue",
                    "bitscore",
                ],
            )
        return df


def RVDB_ref_handler(row):
    compiled = re.compile(r"(\w+)\|(\w+)\|([\w\.\d]+)\|(.*[\.\-\w]+)")
    return compiled.findall(row)[0][2]


def clean_qseqid(row):
    compile_clean = re.compile(r"(\w+\-\w+\-\w+\-\w+\-\w+$)")
    if len(row) < 40:
        return row
    if len(row) > 40:
        cleaned = compile_clean.findall(row)
        if len(cleaned) > 0:
            len_row = len(cleaned[0])
            remove = len_row - 36
            return cleaned[0][remove:]


def extract_sample(sample: str, filename: str) -> pd.DataFrame:
    df = pd.read_csv(f"{filename}")
    date, NA, strain, concentration_CFU, batch, duration_h = sample.split("_")
    batch = re.sub("[^0-9]", "", batch)
    duration_h = re.sub("[^0-9]", "", duration_h)
    CFU_free = concentration_CFU.replace("CFU", "")

    sub_df = df.loc[
        (df["date"] == int(date))
        & (df["NA"] == NA)
        & (df["strain"] == strain)
        & (df["concentration_CFU"] == int(CFU_free))
        & (df["batch"] == int(batch))
        & (df["duration_h"] == int(duration_h))
    ]

    return sub_df


def process(lines=None):
    ks = ["name", "sequence", "optional", "quality"]
    return {k: v for k, v in zip(ks, lines)}


# changed loading algorithm to process fastq file up to 60X faster than pandas!
def generate_fastq_df(fastq_dir: str, read_ids: list) -> pd.DataFrame:
    rows = []
    n = 4
    with open(fastq_dir, "r") as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process(lines)
                rows.append(record)
                lines = []
    clean_df = pd.DataFrame.from_records(rows)
    clean_df = clean_df.rename(
        columns={
            "name": "readID",
            "sequence": "Sequence",
            "optional": "Plus",
            "quality": "Quality",
        }
    )

    clean_df[["readID_A", "readID_B"]] = clean_df["readID"].str.split(
        " ", n=1, expand=True
    )
    clean_df[["readID_C", "readID_D"]] = clean_df["readID_A"].str.split(
        "@", n=1, expand=True
    )
    clean_df = clean_df.drop(["readID_A", "readID_B", "readID_C"], axis=1)
    fq_df = clean_df[clean_df["readID_D"].isin(read_ids)]
    fq_df = fq_df.drop(["readID_D"], axis=1)
    return fq_df


def fastq_formatter(new_fq: str, df: pd.DataFrame) -> None:
    lines = []
    for _index, row in df.iterrows():
        lines.append(row["readID"])
        lines.append(row["Sequence"])
        lines.append(row["Plus"])
        lines.append(row["Quality"])
    if not os.path.isfile(new_fq):
        if len(lines) > 0:
            print(f"Saving to: {new_fq}")
            with open(new_fq, "w") as f:
                f.write("\n".join(lines))


def choose_classifier(
    BLAST_spp_ranked: list, meta_df: pd.DataFrame, unique_ID: str, classifier: str
) -> list:
    flat_list = []
    if len(BLAST_spp_ranked) > 0:
        genera = [x.replace("_", " ") for x in BLAST_spp_ranked]
        TPs = []
        for genus in genera:
            TP = meta_df.loc[meta_df["name"].str.contains(genus)]
            if len(TP) > 0:
                if classifier == "centrifuge":
                    TP = list(TP[unique_ID].unique().astype(int))

                elif classifier == "blastn":
                    TP = list(TP[unique_ID].unique())

                TPs.append(TP)
        flat_list = list(set([item for sublist in TPs for item in sublist]))
    return flat_list


def process_n_save(
    input_df: pd.DataFrame,
    unique_ID: str,
    flat_list: list,
    all_reads: pd.DataFrame,
    binned_reads: str,
    classifier: str,
) -> None:
    TP_reads = list(input_df.loc[input_df[unique_ID].isin(flat_list)].readID)
    FP_reads = list(input_df.loc[~input_df[unique_ID].isin(flat_list)].readID)
    fq_TP_df = generate_fastq_df(all_reads, TP_reads)
    fq_FP_df = generate_fastq_df(all_reads, FP_reads)
    fq_TP_out_name = f"{binned_reads}/TP_reads_{classifier}.fastq"
    fq_FP_out_name = f"{binned_reads}/FP_reads_{classifier}.fastq"
    fastq_formatter(fq_TP_out_name, fq_TP_df)
    fastq_formatter(fq_FP_out_name, fq_FP_df)


def extract_non_host_reads(
    analysis_directory: str,
    sample: str,
    classifier: str,
    meta_df: pd.DataFrame,
    input_df: pd.DataFrame,
    all_reads: pd.DataFrame,
    unique_ID: str,
    BLAST_spp_ranked: list,
) -> None:
    TPs_BLAST_flat = []
    TPs_centrifuge_flat = []

    binned_reads = f"{analysis_directory}/sample_data/{sample}/binned_reads/"
    os.makedirs(binned_reads, exist_ok=True)

    if "centrifuge" in classifier:
        TPs_centrifuge_flat = choose_classifier(
            BLAST_spp_ranked, meta_df, unique_ID, classifier
        )
        process_n_save(
            input_df,
            unique_ID,
            TPs_centrifuge_flat,
            all_reads,
            binned_reads,
            classifier,
        )

    elif "blastn" in classifier:
        TPs_BLAST_flat = choose_classifier(
            BLAST_spp_ranked, meta_df, unique_ID, classifier
        )
        process_n_save(
            input_df, unique_ID, TPs_BLAST_flat, all_reads, binned_reads, classifier
        )


def main(
    samples: list,
    analysis_directory: str,
    meta_data_name: str,
    coverage_list: str,
    kingdom: str,
    bN_name: str,
    no_hosts: str,
    Direct_RNA: bool,
    host_bool: bool,
    BLAST_spp_ranked: list,
) -> None:

    for sample in samples:
        print(
            f"\nCurrent sample: {sample}, \ncoverage data: {coverage_list}, \nmeta data: {meta_data_name}"
        )

        # short list from metagenomic classification only
        meta_df = extract_sample(sample, meta_data_name)

        # short list from coverage analysis
        cov_df = extract_sample(sample, coverage_list)
        cov_df.rename(columns={"species": "name"}, inplace=True)

        troubleshooting = f"{analysis_directory}/sample_data/{sample}/"

        blastn_out = (
            f"{troubleshooting}blastN/trimmed_{sample}{no_hosts}_{bN_name}.hs_bn.tsv"
        )
        centrifuge_out = f"{troubleshooting}centrifuge/{sample}_{kingdom}{no_hosts}_centrifuge_troubleshooting_report.tsv"
        fastq_file = f"{analysis_directory}/sample_data/{sample}/trimmed/tr*q"

        if host_bool:
            fastq_file = f"{analysis_directory}/sample_data/{sample}/trimmed/no_host*q"

        if Direct_RNA:
            blastn_out = f"{troubleshooting}blastN/U2T_trimmed_{sample}{no_hosts}_{bN_name}.hs_bn.tsv"
            fastq_file = f"{analysis_directory}/sample_data/{sample}/trimmed/U*q"

        # load in fastq file and extract desired reads
        fastq_file_loc = glob.glob(fastq_file)
        if len(fastq_file_loc) == 0:
            fastq_file = f"{analysis_directory}/sample_data/{sample}/trimmed/tr*q"
            fastq_file_loc = glob.glob(fastq_file)

        if len(fastq_file_loc) > 0:
            all_reads = fastq_file_loc[0]
            if Path(all_reads).is_file():
                if os.stat(Path(all_reads)).st_size != 0:
                    if Path(centrifuge_out).is_file():
                        centrifuge_filesize = os.stat(centrifuge_out).st_size
                        if centrifuge_filesize > 0:
                            print("Binning reads from centrifuge")
                            # load in centrifuge data
                            cent_df = import_data(centrifuge_out)
                            cent_df = cent_df.loc[cent_df["seqID"] != "unclassified"]
                            cent_df.reset_index(inplace=True, drop=True)
                            # extract metagenomic classified reads
                            extract_non_host_reads(
                                analysis_directory,
                                sample,
                                "centrifuge",
                                meta_df,
                                cent_df,
                                all_reads,
                                "taxID",
                                BLAST_spp_ranked,
                            )

                    if Path(blastn_out).is_file():
                        blastn_filesize = os.stat(blastn_out).st_size
                        if blastn_filesize > 0:
                            print("Binning reads from HS-BLASTn")
                            # load in BLASTn data
                            blast_df = import_data(blastn_out)
                            blast_df.qseqid = blast_df.qseqid.apply(clean_qseqid)
                            blast_df = blast_df[blast_df["sseqid"].notna()]
                            if len(blast_df[blast_df["sseqid"].str.contains("acc")]):
                                blast_df["sseqid"] = blast_df["sseqid"].apply(
                                    RVDB_ref_handler
                                )
                            blast_df = blast_df.rename(columns={"qseqid": "readID"})
                            # extract metagenomic classified reads
                            extract_non_host_reads(
                                analysis_directory,
                                sample,
                                "blastn",
                                meta_df,
                                blast_df,
                                all_reads,
                                "sseqid",
                                BLAST_spp_ranked,
                            )

        else:
            print("No reads to bin")


if __name__ == "__main__":
    # github = "D:/GitHub/SMART-CAMP/"
    # top_dir = "/mnt/usersData/test_DNA/"
    directory = "D:/SequencingData/Harmonisation/DNA/"

    samples = [
        "20211020_DNA_Klebpneu-no-zymo_10000CFU_84_18",
        "20211020_DNA_Klebpneu-no-zymo_1000CFU_84_18",
        "20211020_DNA_Klebpneu-no-zymo_100CFU_84_18",
        "20211020_DNA_Klebpneu-no-zymo_10CFU_84_18",
        "20211020_DNA_Klebpneu-zymo_10000CFU_84_18",
        "20211020_DNA_Klebpneu-zymo_1000CFU_84_18",
        "20211020_DNA_Klebpneu-zymo_100CFU_84_18",
        "20211020_DNA_Klebpneu-zymo_10CFU_84_18",
        "20211020_DNA_TCKlebpneu-no-zymo_10000CFU_86_18",
        "20211020_DNA_TCKlebpneu-no-zymo_1000CFU_86_18",
        "20211020_DNA_TCKlebpneu-no-zymo_100CFU_86_18",
        "20211020_DNA_TCKlebpneu-no-zymo_10CFU_86_18",
        "20211020_DNA_TCKlebpneu-zymo_10000CFU_86_18",
        "20211020_DNA_TCKlebpneu-zymo_1000CFU_86_18",
        "20211020_DNA_TCKlebpneu-zymo_100CFU_86_18",
        "20211020_DNA_TCKlebpneu-zymo_10CFU_86_18",
        "20211021_aDNA_Klebpneu-no-zymo-16S_10000CFU_810_18",
        "20211021_aDNA_Klebpneu-no-zymo-16S_1000CFU_810_18",
        "20211021_aDNA_Klebpneu-no-zymo-16S_100CFU_810_18",
        "20211021_aDNA_Klebpneu-no-zymo-16S_10CFU_810_18",
        "20211021_aDNA_Klebpneu-zymo-16S_10000CFU_89_18",
        "20211021_aDNA_Klebpneu-zymo-16S_1000CFU_89_18",
        "20211021_aDNA_Klebpneu-zymo-16S_100CFU_89_18",
        "20211021_aDNA_Klebpneu-zymo-16S_10CFU_89_18",
        "20211021_aDNA_TCKlebpneu-no-zymo-16S_10000CFU_810_18",
        "20211021_aDNA_TCKlebpneu-no-zymo-16S_1000CFU_810_18",
        "20211021_aDNA_TCKlebpneu-no-zymo-16S_100CFU_810_18",
        "20211021_aDNA_TCKlebpneu-no-zymo-16S_10CFU_810_18",
        "20211021_aDNA_TCKlebpneu-zymo-16S_10000CFU_89_18",
        "20211021_aDNA_TCKlebpneu-zymo-16S_1000CFU_89_18",
        "20211021_aDNA_TCKlebpneu-zymo-16S_100CFU_89_18",
        "20211021_aDNA_TCKlebpneu-zymo-16S_10CFU_89_18",
    ]

    host_bool = False
    Direct_RNA = False
    # barcoded_sample = True
    # no_hosts = ""
    no_hosts = "_no_host"
    analysis_directory = f"{directory}/analysis/"
    # data_output = f"{analysis_directory}/coverage_summary_filter_bacteria_bacteria_filter_bacteria_priority"
    data_output = f"{analysis_directory}/coverage_summary_filter_bacteria_bacteria_16S_23S_no_host_OCS"
    # meta_data_name = "hs_filter_bacteria_priority_combined_cent_blastn.csv"
    meta_data_name = "hs_bacteria_OCS_filter_bacteria_OCS_no_host.csv"
    coverage_list = f"{data_output}/mp_centrifuge_to_coverage_agg*.csv"
    coverage_list = glob.glob(coverage_list)[0]
    meta_data_name_dir = f"{analysis_directory}/{meta_data_name}"
    kingdom = "bacteria"
    bN_name = "filter_bacteria"
    BLAST_spp_ranked = ["Minute_virus"]
    print(
        f"\n{analysis_directory}\n{data_output}\n{meta_data_name}\n{coverage_list}\n{meta_data_name_dir}"
    )

    main(
        samples,
        analysis_directory,
        meta_data_name_dir,
        coverage_list,
        kingdom,
        bN_name,
        no_hosts,
        Direct_RNA,
        host_bool,
        BLAST_spp_ranked,
    )


# time ./pipeline/bin_FP_TP.py
