#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:20:57 2021

@author: mangi
"""

# import argparse
import subprocess
from pathlib import Path
import os
import pandas as pd
from multiprocessing import Pool
from functools import partial

pd.options.mode.chained_assignment = None  # default='warn'

# time /home/james/krakenuniq/krakenuniq/krakenuniq --preload --db /mnt/usersData/krakenDB/vertebrate \
#     --threads 10 --fastq-input "/mnt/usersData/Viral_CHO/analysis/20210420_DNA_CHOK1_1000000CFU_3_42/trimmed/trimmed_20210420_DNA_CHOK1_1000000CFU_3_42.fastq" \
#         --report-file 20210420_DNA_CHOK1_1000000CFU_3_42_REPORTFILE.tsv \
#             --output "20210420_DNA_CHOK1_1000000CFU_3_42_READCLASSIFICATION.tsv"

# High number of kmers is good
# kmers: number of unique k-mers
# dup: average number of times each unique k-mer has been seen
# cov: coverage of the k-mers of the clade in the database


######################################
# handle threads and processes
def threads_n_processes(threads: str) -> (int, str):
    threads = int(threads)

    # control thread count
    if threads > 20:
        threads_per_process = str(int(threads / threads))

    if threads <= 20:
        if threads > 5:
            threads_per_process = str(int(threads / 5))

    processes_to_start = threads
    if processes_to_start > 20:
        processes_to_start = 20

    print(processes_to_start, threads_per_process)
    return processes_to_start, threads_per_process


######################################


def run_krakenuniq(
    krakenuniq_command: str,
    database: str,
    threads: str,
    full_fastq: str,
    report_file: str,
    output: str,
) -> None:
    krakenuniq = subprocess.run(
        [
            krakenuniq_command,
            "--preload",
            "--db",
            database,
            "--threads",
            threads,
            "--fastq-input",
            full_fastq,
            "--report-file",
            report_file,
            "--output",
            output,
        ],
        capture_output=True,
    )
    returncode = krakenuniq.returncode
    if returncode != 0:
        print("KrakenUniq failed to run")


######################################
def process_kraken_reports(sample: str, directory: str, no_host: bool) -> None:
    if not no_host:
        report_file_name = (
            f"{directory}/analysis/sample_data/{sample}/kraken/{sample}_REPORTFILE.tsv"
        )
        report_file_name_out = f"{directory}/analysis/sample_data/{sample}/kraken/{sample}_REPORTFILE_clean.tsv"
        df_out = f"{directory}/analysis/sample_data/{sample}/kraken/{sample}_processed_kreport.csv"

    if no_host:
        report_file_name = f"{directory}/analysis/sample_data/{sample}/kraken/{sample}_host_removed_REPORTFILE.tsv"
        report_file_name_out = f"{directory}/analysis/sample_data/{sample}/kraken/{sample}_host_removed_REPORTFILE_clean.tsv"
        df_out = f"{directory}/analysis/sample_data/{sample}/kraken/{sample}_host_removed_processed_kreport.csv"

    # clean up save file
    with open(report_file_name, "r") as f:
        lines = f.readlines()
        occurrences = [i for i, n in enumerate(lines) if "taxReads" in n]
        if len(occurrences) > 1:
            second_occurrence = occurrences[1]
            cleaned = lines[0:second_occurrence]
            with open(report_file_name_out, "w") as f:
                for line in cleaned:
                    f.write(line)
        else:
            report_file_name_out = report_file_name

    report_df = pd.read_csv(report_file_name, delimiter="\t", skiprows=3)
    report_df["taxName"] = report_df["taxName"].str.strip()

    tax_reads_df = report_df.loc[report_df["taxReads"] > 0]
    print(tax_reads_df)
    tax_reads_df["sample"] = sample
    tax_reads_df[
        ["date", "NA", "strain", "concentration_CFU", "batch", "duration_h"]
    ] = sample.split("_")
    tax_reads_df["concentration_CFU"] = tax_reads_df["concentration_CFU"].replace(
        "CFU", "", regex=True
    )

    tax_reads_df.to_csv(df_out, index=False)


######################################


def multiprocess_sample(
    directory: str,
    krakenuniq_command: str,
    no_host: bool,
    database: str,
    threads: str,
    sample: str,
) -> None:
    print(f"\nRunning sample: {sample}")
    kraken_out = f"{directory}/analysis/sample_data/{sample}/kraken/"  # <-----
    os.makedirs(kraken_out, exist_ok=True)

    # # host reads removed
    if no_host:
        full_fastq = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/no_host_{sample}.fastq"
        )
        report_file = f"{kraken_out}{sample}_host_removed_REPORTFILE.tsv"
        output = f"{kraken_out}{sample}_host_removed_READCLASSIFICATION.tsv"

    # # trimmed reads
    if not no_host:
        full_fastq = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
        )
        report_file = f"{kraken_out}{sample}_REPORTFILE.tsv"
        output = f"{kraken_out}{sample}_READCLASSIFICATION.tsv"

    print(
        f"{krakenuniq_command}, \n{database}, \n{full_fastq}, \n{report_file}, \n{output}\n"
    )

    if Path(full_fastq).is_file():
        if not Path(output).is_file():
            run_krakenuniq(
                krakenuniq_command, database, threads, full_fastq, report_file, output
            )
        process_kraken_reports(sample, directory, no_host)


def run_kraken_host_analysis(
    top_dir: str, threads: str, samples: list, kraken_cov_folder: str
) -> None:
    processes_to_start, threads_per_process = threads_n_processes(threads)
    kraken_coverage = f"{top_dir}/analysis/{kraken_cov_folder}/"
    os.makedirs(kraken_coverage, exist_ok=True)

    ######################################
    # MAY REQUIRE EDITING IF NEW SETUP
    krakenuniq_command = "/home/james/krakenuniq/krakenuniq/krakenuniq"  # <-----

    database = "/mnt/usersData/krakenDB/vertebrate/"  # <-----
    ######################################

    # generate kraken folder, find filenames, save files
    func_host = partial(
        multiprocess_sample,
        top_dir,
        krakenuniq_command,
        False,
        database,
        threads_per_process,
    )

    with Pool(processes=processes_to_start) as p:
        p.map(func_host, samples)  # process data_inputs iterable with pool

    func_host_removed = partial(
        multiprocess_sample,
        top_dir,
        krakenuniq_command,
        True,
        database,
        threads_per_process,
    )

    with Pool(processes=processes_to_start) as p:
        p.map(func_host_removed, samples)  # process data_inputs iterable with pool

    print("\nMultiprocessing and kraken coverage analysis complete")
    kraken_dfs = []
    kraken_dfs_nh = []
    kraken_df_nh_path = ""
    for sample in samples:
        kraken_df_path = f"{top_dir}/analysis/sample_data/{sample}/kraken/{sample}_processed_kreport.csv"
        kraken_df = pd.read_csv(kraken_df_path)
        kraken_dfs.append(kraken_df)
        if "TC" in sample or "CHOK1" in sample or "Jurkat" in sample:
            kraken_df_nh_path = f"{top_dir}/analysis/sample_data/{sample}/kraken/{sample}_host_removed_processed_kreport.csv"
            kraken_df_nh = pd.read_csv(kraken_df_nh_path)
            kraken_dfs_nh.append(kraken_df_nh)
    concat_df = pd.concat(kraken_dfs)
    concat_df.to_csv(f"{kraken_coverage}/aggregate_kraken_coverage.csv", index=False)
    print(f"\n{kraken_coverage}/aggregate_kraken_coverage.csv\n")

    if len(kraken_dfs_nh) > 0:
        concat_df_nh = pd.concat(kraken_dfs_nh)
        concat_df_nh.to_csv(
            f"{kraken_coverage}/aggregate_host_removed_kraken_coverage.csv", index=False
        )
        print(f"\n{kraken_coverage}/aggregate_host_removed_kraken_coverage.csv\n")


if __name__ == "__main__":
    top_dir = "/mnt/usersData/Viral_CHO/"
    threads = "20"
    samples = [
        "20211018_ssDNA_MVM1-100CHOK1-Concentrate_50000CFU_20_36",
        "20211018_ssDNA_MVM1-100CHOK1-Supernatant_50000CFU_20_36",
        "20211018_ssDNA_MVM1-10CHOK1-Concentrate_500000CFU_20_36",
        "20211018_ssDNA_MVM1-10CHOK1-Filtrate_500000CFU_20_36",
        "20211018_ssDNA_MVM1-10CHOK1-Supernatant_500000CFU_20_36",
    ]
    kraken_cov_folder = "kraken_coverage_host_Cricetulus_griseus"
    run_kraken_host_analysis(top_dir, threads, samples, kraken_cov_folder)
