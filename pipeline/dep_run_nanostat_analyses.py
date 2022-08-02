#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 09:21:08 2022

@author: mangi
"""

from pipeline import cat_nanostats, find_barcode_reads_for_QC
from multiprocessing import Pool
from functools import partial
import subprocess
import pandas as pd
import glob
import os
import argparse
from pathlib import Path


# handle threads and processes
def threads_n_processes(threads: str) -> (int, str):
    threads = int(threads)

    # control thread count
    if threads > 25:
        threads_per_process = str(int(threads / threads))
        # print(threads, threads_per_process, threads*int(threads_per_process))

    if threads <= 25 and threads > 2:
        threads_per_process = str(int((threads**-1) * 60))
        # print(threads, threads_per_process, threads*int(threads_per_process))

    if threads <= 2:
        threads_per_process = 20
        # print(threads, threads_per_process, threads*int(threads_per_process))

    processes_to_start = threads
    if processes_to_start > 25:
        processes_to_start = 25

    total_threads = threads * int(threads_per_process)
    print(
        f"Processes to start: {processes_to_start}, threads per process: {threads_per_process}, total threads: {total_threads}"
    )
    return processes_to_start, threads_per_process


# NanoPlot -t 20 --summary ~/SequencingData/Direct_RNA/{sample}/**/sequencing_summary_*.txt --loglength -o ~/SequencingData/Direct_RNA/analysis/sample_data/{sample}/summary-plots-log-transformed
def nanoplot(directory: str, threads: str, sample: str) -> None:
    seq_summ = glob.glob(f"{directory}/{sample}/**/sequencing_summary*.txt")
    # print(f'{directory}/{sample}/**/sequencing_summary*.txt')
    if len(seq_summ) > 0:
        seq_summ = seq_summ[0]
        plot_out = (
            f"{directory}/analysis/sample_data/{sample}/summary-plots-log-transformed/"
        )
        os.makedirs(plot_out, exist_ok=True)
        print(plot_out)
        if not Path(f"{plot_out}NanoStats.txt").is_file():
            subprocess.run(
                [
                    "NanoPlot",
                    "-t",
                    str(threads),
                    "--summary",
                    seq_summ,
                    "--loglength",
                    "-o",
                    plot_out,
                ],
                capture_output=True,
            )


# takes in pd.df with at least 'date', 'NA', 'strain', 'concentration_CFU', 'batch', 'duration_h'
def get_barcodes(df: pd.DataFrame) -> (list, list):
    print(df)
    dates = df["date"].to_list()
    NAs = df["NA"].to_list()
    strains = df["strain"].to_list()
    CFUs = df["concentration_CFU"].to_list()
    batches = df["batch"].to_list()
    time = df["duration_h"].to_list()
    sample_id = [list(e) for e in zip(dates, NAs, strains, CFUs, batches, time)]
    identities = df["Identifier"].to_list()
    barcodes = df["Barcode"].to_list()

    final_samples = []
    for sample in sample_id:
        intermed_sample = []
        for item in sample:
            intermed_sample.append(str(item))
        intermed_sample = "_".join(intermed_sample)
        final_samples.append(intermed_sample)

    folder_builder = [list(e) for e in zip(final_samples, identities, barcodes)]
    return final_samples, folder_builder


def run_nanostat(directory: str, threads: int, experiment_for_processing: str) -> None:
    non_barcoded = False
    barcoded = False

    processes_to_start, threads_per_process = threads_n_processes(threads)

    if ".txt" in experiment_for_processing:
        non_barcoded = True
    elif ".csv" in experiment_for_processing:
        barcoded = True

    print(f"Barcoded: {barcoded}; non-barcoded: {non_barcoded}")

    # for non-barcoded

    if non_barcoded:
        samples = [
            line.rstrip("\n").split(",") for line in open(experiment_for_processing)
        ]
        samples = [item for sublist in samples for item in sublist]

        mp_func = partial(nanoplot, directory, threads_per_process)
        with Pool(processes=processes_to_start) as p:
            p.map(mp_func, samples)

    # for barcoded

    if barcoded:
        barcode_df = pd.read_csv(experiment_for_processing)
        samples, builder = get_barcodes(barcode_df)

        builder_dict = {}
        for i in builder:
            if len(builder_dict) == 0:
                builder_dict[i[1]] = [i[0]]
                j = i[1]
            if i[1] == j:
                builder_dict[i[1]].append(i[0])
            else:
                if i[1] != j:
                    builder_dict[i[1]] = [i[0]]
                    j = i[1]

        identifiers = []
        dict_identifiers = {}
        for index, sample in enumerate(samples):
            if builder[index][0] == sample:
                identifier = builder[index][1]
                identifiers.append(identifier)
                value_id = identifier.replace("B", "").replace("S", "")
                dict_identifiers[identifier] = f"_{value_id}_"

        # check if batch exists
        unfounds = []
        for ex_key, ex_val in dict_identifiers.items():
            fastq_globs = find_barcode_reads_for_QC.find_fastqs(ex_val, directory)
            print(fastq_globs)
            if len(fastq_globs) > 0:
                pass

            elif len(fastq_globs) == 0:
                unfounds.append(ex_val)

        not_found = list(set(unfounds))

        to_use_sample_name = []
        for item in not_found:
            for k, v in dict_identifiers.items():
                if item == v:
                    to_use_sample_name.append(k)

        for k, v in builder_dict.items():
            for items in to_use_sample_name:
                if items == k:
                    dict_identifiers[k] = v

        print("\nsecond pass:", dict_identifiers)
        find_barcode_reads_for_QC.run_extract_reads_for_qc(directory, dict_identifiers)

        mp_func = partial(nanoplot, directory, threads_per_process)
        with Pool(processes=processes_to_start) as p:
            p.map(mp_func, samples)

    # # generate cat file
    try:
        cat_nanostats.get_nanostats_df(directory, samples)
    except:
        print(
            "NanoStats.txt not available, sequencing_summary_xxx file is required to run this."
        )


# time ~/SMART-CAMP/run_nanostat_analyses.py -d /mnt/usersData/Viral_CHO/  -t 10 -e "/home/james/SMART-CAMP/configs/viral_DNA_all.txt"


if __name__ == "__main__":
    directory = "/mnt/usersData/DNA/"
    threads = 5
    barcoded_samples = "/home/james/SMART-CAMP/configs/aDNA_all4.csv"

    # run_nanostat(directory, threads, barcoded_samples)
