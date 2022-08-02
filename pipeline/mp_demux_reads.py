#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 12:49:46 2021

@author: mangi
"""

import glob
import os
import subprocess
import time
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from pipeline import dep_pre_chunker  # turn this on if script is a child

# import dep_pre_chunker # turn this on if script is parent

encoding = "utf-8"

######################################
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


######################################

# complete demultiplexing using qcat software
# cat ~/SequencingData/TCPA/S7B6/20210105_0309_MN34889_FAO93557_e5ac2be1/fastq_pass/*.fastq |
# qcat -b ~/SequencingData/TCPA/analysis/sample_data/S7B6/demux/ --threads 60 --min-read-length 500
# -k RBK004 --tsv > ~/SequencingData/TCPA/analysis/sample_data/S7B6/demux/qcat_log.tsv
def multiprocess_sample(
    threads_per_process: str,
    directory: str,
    experiment: str,
    chunked_files: str,
    kit: str,
    fq_type: str,
    glob_fq: str,
) -> None:
    # chunked_files = "fastq_processed"
    # glob_fq = f"fastq_processed/pass/c*.fastq"
    filename = glob_fq.split("/")[-1].split(".fastq")[0]  # will give e.g. chunk0

    # find identifier; not for each individual barcode but for experiment
    out_dir = (
        f"{directory}/analysis/sample_data/{experiment}/demux/{filename}-{fq_type}"
    )

    # # capture the log files
    out_devnull = f"{directory}/analysis/sample_data/{experiment}/demux/{filename}-{fq_type}/qcat_log.tsv"

    if Path(out_devnull).is_file():
        print(f"{out_devnull} already exists, skipping.")

    if not Path(out_devnull).is_file():
        print(
            f"merged_file: {glob_fq}; out_dir: {out_dir}; threads{threads_per_process}"
        )
        qcat_run = subprocess.run(
            [
                "qcat",
                "-f",
                glob_fq,
                "-b",
                out_dir,
                "--threads",
                threads_per_process,
                "--min-read-length",
                "200",
                "-k",
                kit,
                "--tsv",
            ],
            capture_output=True,
        )

        returncode = qcat_run.returncode
        if returncode != 0:
            print("No qcat demultiplexed barcoded reads / files found")
        with open(out_devnull, "w") as f:
            f.write(str(qcat_run.stdout, encoding))


def setup_mp(
    directory: str,
    sample: str,
    processes_to_start: int,
    threads_per_process: str,
    globbed_fq: list,
    chunked_files: str,
    kit: str,
    fq_type: str,
) -> None:
    if len(globbed_fq) > 0:
        func = partial(
            multiprocess_sample,
            threads_per_process,
            directory,
            sample,
            chunked_files,
            kit,
            fq_type,
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


def qcat_trim(
    directory: str,
    experiment: str,
    processes_to_start: int,
    threads_per_process: str,
    meta_fastq_fail: bool,
    kit: str,
) -> None:
    epoch_time_start = int(time.time())

    epoch_dir = f"{directory}/analysis/run/{str(epoch_time_start)}"
    epoch_dir_pc = f"{directory}/analysis/run/{str(epoch_time_start)}/porechop"
    os.makedirs(f"{directory}/analysis/run/", exist_ok=True)
    os.makedirs(epoch_dir, exist_ok=True)
    os.makedirs(epoch_dir_pc, exist_ok=True)

    print(f"De-multiplexing with qcat...{experiment}")
    all_fq_dir = f"{directory}/{experiment}/"
    all_fq = f"{all_fq_dir}/**/*q_p*s/F*.fastq"
    print(all_fq)
    globbed_fq = glob.glob(all_fq)
    print(globbed_fq)
    folder = ""
    chunked_files = ""
    if len(globbed_fq) > 0:
        folder = "/".join(globbed_fq[0].split("/")[:-2])
        chunked_files = f"{folder}/fastq_processed/"
        fq_type_pass = "pass"
        processed_pass = f"{folder}/fastq_processed/{fq_type_pass}"
        print(f"fastq pass: {processed_pass}")
        fq_pass = f"{folder}/fastq_pass/"
        os.makedirs(chunked_files, exist_ok=True)
        os.makedirs(processed_pass, exist_ok=True)
        dep_pre_chunker.chunk_fq_files(fq_pass, processed_pass, True, globbed_fq)
        chunk_run, epoch_time_chunk = time_step(
            epoch_time_start, "Chunking", experiment
        )

        chunk_pass_glob = glob.glob(f"{processed_pass}/c*.fastq")
        setup_mp(
            directory,
            experiment,
            processes_to_start,
            threads_per_process,
            chunk_pass_glob,
            chunked_files,
            kit,
            fq_type_pass,
        )

    if meta_fastq_fail:
        fail_fq = f"{all_fq_dir}/**/*q_f*/F*.fastq"
        fq_fail = f"{folder}/fastq_fail/"
        fq_type_fail = "fail"
        processed_fail = f"{folder}/fastq_processed/{fq_type_fail}"
        print(f"fastq fail: {processed_fail}; have you specified 'guppy-demux' in the Sample_original field??")
        os.makedirs(processed_fail, exist_ok=True)
        globbed_fail_fq = glob.glob(fail_fq)
        dep_pre_chunker.chunk_fq_files(fq_fail, processed_fail, True, globbed_fail_fq)
        chunk_run, epoch_time_chunk = time_step(
            epoch_time_start, "Chunking", experiment
        )

        chunk_fail_glob = glob.glob(f"{processed_fail}/c*.fastq")
        setup_mp(
            directory,
            experiment,
            processes_to_start,
            threads_per_process,
            chunk_fail_glob,
            chunked_files,
            kit,
            fq_type_fail,
        )

    qcat_run, epoch_time_end = time_step(epoch_time_start, "Qcat-demux", experiment)

    if len(globbed_fq) > 0:
        # find all PORECHOP PROCESSED fastq files
        chunked_demux_files = f"{directory}/analysis/sample_data/{experiment}/demux/"
        chunked_demux_files_to_glob = glob.glob(f"{chunked_demux_files}/**/b*.fastq")

        # barcode_names = [i.replace("\\","/").split("/")[-1].replace(".fastq", "") for i in chunked_demux_files_to_glob]
        barcode_names = list(
            set(
                [
                    i.split("/")[-1].replace(".fastq", "")
                    for i in chunked_demux_files_to_glob
                ]
            )
        )

        for barcode_name in barcode_names:
            sub_fq = [bc for bc in chunked_demux_files_to_glob if barcode_name in bc]
            file = f"{chunked_demux_files}/{barcode_name}.fastq"
            print(f"Saving merged barcode fq to: {file}")
            if not Path(file).is_file():
                with open(file, "w") as outfile:
                    for infilename in sub_fq:
                        with open(infilename) as infile:
                            outfile.write(infile.read())
            print(f"Barcoded fastq: {file}")

    else:
        print("Qcat did not generate any barcoded files.")

    merge_trim_run, epoch_time_merge = time_step(
        epoch_time_end, "Merging of demuxed files", experiment
    )

    times_taken = [qcat_run, merge_trim_run]

    # save sample name for later checking
    with open(f"{epoch_dir_pc}/demux_time_{experiment}.txt", "w") as text_file:
        for i in times_taken:
            text_file.write(f"{str(i)}\n")


def main(
    threads: str, directory: str, experiment: str, meta_fastq_fail: bool, kit: str
) -> None:
    processes_to_start, threads_per_process = threads_n_processes(threads)
    qcat_trim(
        directory,
        experiment,
        processes_to_start,
        threads_per_process,
        meta_fastq_fail,
        kit,
    )


if __name__ == "__main__":
    top_dir = "/mnt/usersData/Fungal_test/"
    experiment = "S8B5"
    meta_fastq_fail = True
    basecalled_already = False
    barcoded_sample = False
    threads = "20"
    kit = "RBK004"

    main(threads, top_dir, experiment, meta_fastq_fail, kit)
