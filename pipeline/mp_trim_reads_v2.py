#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 12:49:46 2021

@author: mangi
"""

import glob
import os
import re
import subprocess
import time
import json
import shutil
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from pipeline import dep_pre_chunker

encoding = "utf-8"


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


##############################################################################
# porechop -i /home/james/SequencingData/Direct_RNA/{sample}/**/fastq_pass/merged.fastq
# -o /home/james/SequencingData/Direct_RNA/analysis/{sample}/trimmed/trimmed_{sample}.fastq
# --threads 16 > ~/SequencingData/Direct_RNA/analysis/{sample}/trimmed/porechop_log.tsv
def multiprocess_sample(
    threads_per_process: str,
    directory: str,
    sample: str,
    chunked_files: str,
    fq_type: str,
    demux_temp_folder: str,
    barcoded_sample: bool,
    glob_fq: str,
) -> None:
    # chunked_files = "fastq_processed"
    # glob_fq = f"fastq_processed/pass/c*.fastq"
    filename = glob_fq.split("/")[-1].split(".fastq")[0]
    out_fq = f"{chunked_files}/{filename}_{fq_type}.fastq"
    if barcoded_sample:
        out_fq = f"{demux_temp_folder}/{filename}-T.fastq"

    if not Path(out_fq).is_file():
        print(f"Saving to: {out_fq}")
        porechop_run = subprocess.run(
            ["porechop", "-i", glob_fq, "-o", out_fq, "-t", threads_per_process],
            capture_output=True,
        )


def setup_mp(
    directory: str,
    sample: str,
    processes_to_start: int,
    threads_per_process: str,
    globbed_fq: list,
    chunked_files: str,
    fq_type: str,
    demux_temp_folder: str,
    barcoded_sample: bool,
) -> None:
    if len(globbed_fq) > 0:
        func = partial(
            multiprocess_sample,
            threads_per_process,
            directory,
            sample,
            chunked_files,
            fq_type,
            demux_temp_folder,
            barcoded_sample,
        )
        with Pool(processes=processes_to_start) as p:
            p.map(func, globbed_fq)  # process data_inputs iterable with pool

    else:
        print("No fastq files available for porechop to trim.")


def time_step(start_time: int, time_str: str, sample: str) -> (str, int):
    next_time_point = int(time.time())
    elapsed_time = next_time_point - start_time
    elapsed_time_min = elapsed_time / 60
    time_check = f"{time_str} runtime for {sample}: {elapsed_time} seconds or {elapsed_time_min} minutes.\n"
    print(time_check)

    return time_check, next_time_point


def save_to_trim(file: str, trimmed_files_to_glob: list) -> None:
    # check if trimmed fq file already exists
    if not Path(file).is_file():
        with open(file, "w") as outfile:
            for infilename in trimmed_files_to_glob:
                with open(infilename) as infile:
                    outfile.write(infile.read())
    print(f"Merged trimmed fastq: {file}")


def porechop_trim(
    directory: str,
    sample: str,
    processes_to_start: int,
    threads_per_process: str,
    meta_fastq_fail: bool,
    barcoded_sample: bool,
    epoch_time_start: int,
) -> None:

    # epoch_dir = f"{directory}/analysis/run/{str(epoch_time_start)}"
    epoch_dir_pc = f"{directory}/analysis/run/{str(epoch_time_start)}/porechop"

    os.makedirs(epoch_dir_pc, exist_ok=True)

    print("Trimming with porechop...")
    all_fq_dir = f"{directory}/{sample}/"
    all_fq = f"{all_fq_dir}/**/*q_p*s/F*.fastq"
    globbed_fq = glob.glob(all_fq)

    if len(globbed_fq) == 0:
        try:
            all_fq = f"{all_fq_dir}/**/*q_p*s/*.fastq"
            globbed_fq = glob.glob(all_fq)
            meta_fastq_fail = False
        except:
            print("No files found.")

    folder = ""
    chunked_files = ""
    fq_type_pass = ""
    fq_type_fail = ""

    # # requires some form of barcode check #
    if barcoded_sample:
        demux_temp_folder = f"{directory}/{sample}/**/fastq_pass/"
        demux_barcoded_fqs = f"{demux_temp_folder}/*fastq"
        globbed_fq = glob.glob(demux_barcoded_fqs)
        analysis = "Chunk and porechop"
        if len(globbed_fq) > 0:
            folder = "/".join(globbed_fq[0].split("/")[-3:-2])
            demux_temp_folder = f"{directory}/{sample}/{folder}/fastq_pass/"
            # dep_pre_chunker.chunk_fq_files(demux_temp_folder, demux_temp_folder, True, globbed_fq)
            chunk_run, epoch_time_chunk = time_step(
                epoch_time_start, "Pre-trim setup", sample
            )
            chunk_demux_glob = glob.glob(f"{demux_temp_folder}/2*.fastq")
            setup_mp(
                directory,
                sample,
                processes_to_start,
                threads_per_process,
                chunk_demux_glob,
                chunked_files,
                fq_type_pass,
                demux_temp_folder,
                barcoded_sample,
            )

            porechop_run, epoch_time_end = time_step(
                epoch_time_start, "Porechop post-demux", sample
            )

            # find all PORECHOP PROCESSED fastq files
            trimmed_files_to_glob = glob.glob(f"{demux_temp_folder}/2*-T.fastq")
            # stop replication of sample-T-T-....fastq
            trimmed_files_to_glob = [
                url for url in trimmed_files_to_glob if f"{sample}-T.fastq" in url
            ]
            file = f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
            save_to_trim(file, trimmed_files_to_glob)

            merge_trim_run, epoch_time_merge = time_step(
                epoch_time_end, "Merging of trimmed files", sample
            )

            times_taken = [porechop_run, merge_trim_run]

            with open(
                f"{epoch_dir_pc}/trim_time_{sample}.txt", "w", encoding="utf-8"
            ) as f:
                json.dump(times_taken, f)
        else:
            print("Porechop did not generate any trimmed files.")

    # # Non-barcoded samples
    if not barcoded_sample:
        # demux_temp_folder = f'{directory}/{sample}/**/fastq_pass/'
        analysis = "Chunking-only"
        if len(globbed_fq) > 0:
            folder = "/".join(globbed_fq[0].split("/")[:-2])
            chunked_files = f"{folder}/fastq_processed/"
            fq_type_pass = "pass"
            processed_pass = f"{folder}/fastq_processed/{fq_type_pass}"
            fq_pass = f"{folder}/fastq_pass/"
            os.makedirs(chunked_files, exist_ok=True)
            os.makedirs(processed_pass, exist_ok=True)
            dep_pre_chunker.chunk_fq_files(fq_pass, processed_pass, True, globbed_fq)

        if meta_fastq_fail:
            fail_fq = f"{all_fq_dir}/**/*q_f*/F*.fastq"
            fq_type_fail = "fail"
            globbed_fail_fq = glob.glob(fail_fq)

            if len(globbed_fail_fq) > 0:
                fq_fail = f"{folder}/fastq_fail/"
                processed_fail = f"{folder}/fastq_processed/{fq_type_fail}"
                os.makedirs(processed_fail, exist_ok=True)
                dep_pre_chunker.chunk_fq_files(
                    fq_fail, processed_fail, True, globbed_fail_fq
                )

        porechop_run, epoch_time_end = time_step(epoch_time_start, analysis, sample)


def main(
    threads: str,
    directory: str,
    sample: str,
    meta_fastq_fail: bool,
    barcoded_sample: bool,
    epoch_time_start: int,
) -> None:
    processes_to_start, threads_per_process = threads_n_processes(threads)
    porechop_trim(
        directory,
        sample,
        processes_to_start,
        threads_per_process,
        meta_fastq_fail,
        barcoded_sample,
        epoch_time_start,
    )


##############################################################################

##############################################################################
# these functions are called independently by mp_metagenomic_assessment_v2.py
def mp_run_porechop_sample_chunk(
    fastq_processed_dir: str,
    fq_type_pass: str,
    samples: list,
    threads_per_process: str,
    processes_to_start: int,
    top_dir: str,
) -> None:
    processed_pass = f"{fastq_processed_dir}/{fq_type_pass}"
    chunk_pass_glob = glob.glob(f"{processed_pass}/c*.fastq")

    check_samples = [
        cpg for cpg in chunk_pass_glob if any([True for s in samples if s in cpg])
    ]

    func = partial(
        multiprocess_sample_chunk, threads_per_process, top_dir, samples, fq_type_pass
    )
    with Pool(processes=processes_to_start) as p:
        p.map(func, check_samples)  # process data_inputs iterable with pool


def multiprocess_sample_chunk(
    threads_per_process: str, directory: str, samples: list, fq_type: str, glob_fq: str
) -> None:
    # chunked_files = "fastq_processed"
    # glob_fq = f"fastq_processed/pass/c*.fastq"
    sample = [sample for sample in samples if sample in glob_fq]

    if len(sample) > 0:
        sample = sample[0]

    filename = glob_fq.split("/")[-1].split(".fastq")[0]
    folder = "/".join(glob_fq.split("/")[:-2])
    # chunked_sample_out = f"{directory}/{sample}/"
    out_fq = f"{folder}/{filename}_{fq_type}.fastq"

    if not Path(out_fq).is_file():
        print(f"Opening file from: {glob_fq}, Saving to: {out_fq}")
        subprocess.run(
            ["porechop", "-i", glob_fq, "-o", out_fq, "-t", threads_per_process],
            capture_output=True,
        )


def main_two(threads: str, samples: list, directory: str, meta_fastq_fail: bool):
    processes_to_start, threads_per_process = threads_n_processes(threads)
    folder = f"{directory}/**/**/"
    chunked_files = f"{folder}/fastq_processed/"
    fq_type_pass = "pass"
    mp_run_porechop_sample_chunk(
        chunked_files,
        fq_type_pass,
        samples,
        threads_per_process,
        processes_to_start,
        directory,
    )

    if meta_fastq_fail:
        fq_type_fail = "fail"
        mp_run_porechop_sample_chunk(
            chunked_files,
            fq_type_fail,
            samples,
            threads_per_process,
            processes_to_start,
            directory,
        )

    trimmed_files_to_glob = []
    for sample in samples:
        all_fq_dir = f"{directory}/{sample}/"
        all_fq = f"{all_fq_dir}/**/fastq_processed/c*_*.fastq"
        trimmed_files_to_glob = glob.glob(all_fq)
        if len(trimmed_files_to_glob) > 0:
            file = f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
            save_to_trim(file, trimmed_files_to_glob)

        else:
            print("Porechop did not generate any trimmed files.")

    # delete intermediate chunked files in pass/fail folders to save space
    unique_id = re.compile(r"(\w+)\/fastq_pass")
    for sample in samples:
        if len(trimmed_files_to_glob) > 0:
            directory2delete = trimmed_files_to_glob[0]
            is_runid = unique_id.findall(directory2delete)
            if len(is_runid) > 0:
                process_pass = (
                    f"{directory}/{sample}/{is_runid}/fastq_processed/{fq_type_pass}"
                )
                process_fail = (
                    f"{directory}/{sample}/{is_runid}/fastq_processed/{fq_type_fail}"
                )
                if os.path.isdir(process_pass):
                    shutil.rmtree(process_pass)
                if os.path.isdir(process_fail):
                    shutil.rmtree(process_fail)


##############################################################################

if __name__ == "__main__":
    top_dir = "E:/SequencingData/Viral_CHO/"
    sample = "20210902_ssDNA_MVM100-1CHOK1-Filt-NoDNase-NoRapid_500000000CFU_16_36"
    meta_fastq_fail = True
    basecalled_already = False
    barcoded_sample = False
    threads = "20"
    epoch_time_start = 16841564

    main(threads, top_dir, sample, meta_fastq_fail, barcoded_sample, epoch_time_start)
    # main_two(threads, samples, directory, meta_fastq_fail)
