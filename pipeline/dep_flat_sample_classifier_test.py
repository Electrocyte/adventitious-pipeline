#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 11:34:07 2022

@author: mangi
"""

from multiprocessing import Pool
from functools import partial
from pipeline import database_config, updated_blastn_interpreter_v2
from pathlib import Path
import pandas as pd
import tempfile
import subprocess
import argparse
import shutil
import glob
import time
import json
import os


encoding = "utf-8"
######################################
# handle threads and processes
def threads_n_processes(processes: int) -> (int, str):
    processes = int(processes)
    threads_per_process = 1

    # control thread count
    if processes > 20:
        threads_per_process = int(processes / processes)

    if processes <= 20:
        if processes > 5:
            threads_per_process = int(processes / 5)

    if processes <= 5:
        if processes > 1:
            threads_per_process = int(processes * 5)

    processes_to_start = processes
    if processes_to_start > 20:
        processes_to_start = 20

    print(processes_to_start, threads_per_process)
    return processes_to_start, threads_per_process


######################################


def gen_fasta(fasta: str, in_folder_glob: str) -> None:
    if not Path(fasta).is_file() or Path(fasta).stat().st_size == 0:
        fastq2fasta = subprocess.run(
            ["sed", "/^@/!d;s//>/;N", in_folder_glob], capture_output=True
        )
        returncode = fastq2fasta.returncode
        if returncode != 0:
            print("Unable to generate fasta from fastq")
        with open(fasta, "w") as f:
            print(f"Saving to: {fasta}")
            f.write(str(fastq2fasta.stdout, encoding))


def check_file_exists(file, existing_file) -> str:
    f_glob = glob.glob(file)
    if len(f_glob) > 0:
        return file
    else:
        return existing_file


# calculate number of reads in fasta
def file_len(fname: str) -> int:
    i = 0
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def read_data(file):
    try:
        df = pd.read_csv(file, delimiter="\t", skiprows=3)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()

    return df


def interpret_blastn(
    directory: str,
    sample: str,
    search_term: str,
    plots: bool,
    bN_name: str,
    no_host: str,
    github: str,
) -> int:

    print(directory, sample, search_term, False, bN_name, no_host, github)

    blastn_dir = f"{directory}/analysis/sample_data/{sample}/blastN/"

    # get all human, bacteria sequence ids
    seqid_dir = f"{github}/all_seqids.csv"
    if bN_name == "cviral":
        seqid_dir = f"{github}/C_RVDB_seqids.csv"

    # THIS WORKS FINE #############################
    # hs_classified_reads = HS_blastn_interpreter.main(blastn_dir, seqid_dir, plots, search_term, bN_name, no_host, sample)

    # THIS IS FOR TESTING #############################
    sseqids = updated_blastn_interpreter_v2.import_deanonymised_nomenclature(seqid_dir)
    print(sseqids)

    hs_classified_reads = 0
    try:
        hs_classified_reads = updated_blastn_interpreter_v2.main(
            blastn_dir, sseqids, plots, search_term, bN_name, no_host, sample, directory
        )
    except:
        print(f"Current sample contains no classified reads for BLAST: {sample}")
        hs_classified_reads = 0

    return hs_classified_reads


def interpret_centrifuge(
    sample: str, troubleshooting: str, directory: str, kingdom: str
) -> None:

    skip_describe = False

    fastq_file = f"{directory}/{sample}/**/sequencing_summ*"
    fastq_seq_summ_glob = glob.glob(fastq_file)
    print(f"Searching for qscores here: {fastq_seq_summ_glob}")
    qscores = pd.DataFrame()
    if len(fastq_seq_summ_glob) > 0:
        fastq_summary_df = pd.read_csv(fastq_seq_summ_glob[0], delimiter="\t")
        qscores = fastq_summary_df[["read_id", "mean_qscore_template"]]

    analysis_directory = f"{directory}/analysis/"
    summary_centrifuge_file = f"{analysis_directory}/sample_data/{sample}/centrifuge/describe_{kingdom}_out.csv"
    print(summary_centrifuge_file, sample)

    # for old files, will require turning this off

    if skip_describe:
        if not os.path.isfile(summary_centrifuge_file):
            print(f"Current file: {troubleshooting}")
            line_count = file_len(troubleshooting)
            if line_count > 1:
                df = pd.read_csv(troubleshooting, delimiter="\t")
                df = df.rename(columns={"readID": "read_id"})

                if len(qscores) > 0:
                    df = df.merge(qscores, on="read_id")

                df = df.loc[df["seqID"] != "unclassified"]
                if not df.empty:
                    try:
                        stats = df.groupby(["seqID", "taxID"])[
                            [
                                "score",
                                "2ndBestScore",
                                "hitLength",
                                "queryLength",
                                "numMatches",
                                "mean_qscore_template",
                            ]
                        ].describe()
                    except:
                        stats = df.groupby(["seqID", "taxID"])[
                            [
                                "score",
                                "2ndBestScore",
                                "hitLength",
                                "queryLength",
                                "numMatches",
                            ]
                        ].describe()
                    stats.columns = ["_".join(a) for a in stats.columns.to_flat_index()]
                    stats.reset_index(inplace=True)
                    stats["sample"] = sample
                    stats["sample"] = stats["sample"].replace(
                        to_replace=r"CFU", value="", regex=True
                    )
                    try:
                        stats.to_csv(summary_centrifuge_file, index=False)
                    except:
                        print("Failed to save centrifuge additional stats file.")
    else:
        print(f"Current file: {troubleshooting}")
        if os.path.isfile(troubleshooting):
            print(f"q-scores from sequencing summary: {len(qscores)}")
            print(f"File size: {os.stat(troubleshooting).st_size}")
            line_count = file_len(troubleshooting)
            if line_count > 1:
                df = pd.read_csv(troubleshooting, delimiter="\t")
                print(df)
                df = df.rename(columns={"readID": "read_id"})

                if len(qscores) > 0:
                    df = df.merge(qscores, on="read_id")

                df = df.loc[df["seqID"] != "unclassified"]
                if not df.empty:
                    try:
                        stats = df.groupby(["seqID", "taxID"])[
                            [
                                "score",
                                "2ndBestScore",
                                "hitLength",
                                "queryLength",
                                "numMatches",
                                "mean_qscore_template",
                            ]
                        ].describe()
                    except:
                        stats = df.groupby(["seqID", "taxID"])[
                            [
                                "score",
                                "2ndBestScore",
                                "hitLength",
                                "queryLength",
                                "numMatches",
                            ]
                        ].describe()
                    stats.columns = ["_".join(a) for a in stats.columns.to_flat_index()]
                    stats.reset_index(inplace=True)
                    stats["sample"] = sample
                    stats["sample"] = stats["sample"].replace(
                        to_replace=r"CFU", value="", regex=True
                    )
                    try:
                        stats.to_csv(summary_centrifuge_file, index=False)
                    except:
                        print("Failed to save centrifuge additional stats file.")


def run_centrifuge(
    directory: str,
    sample_classifier: str,
    threads: str,
    clf_idx: str,
    bN_name: str,
    new_fqs: str,
    no_host: str,
    epoch_dir_host: str,
) -> None:
    current_sample_classifier = sample_classifier.split("**")
    sample, classifier, kingdom = current_sample_classifier
    print(sample)
    print(classifier)

    new_fq = [s for s in new_fqs if sample in s]

    fastq4cent = (
        f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
    )
    fastq_folder = f"{directory}/analysis/sample_data/{sample}/trimmed/"
    read_count = fastq4cent
    if len(no_host) > 0:
        fastq4cent_putative = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/no_host_{sample}.fastq"
        )
        fastq4cent = check_file_exists(fastq4cent_putative, fastq4cent)

    if len(new_fq) > 0:
        fastq4cent = new_fq[0]
        fastq_folder = "/".join(fastq4cent.split("/")[:-1])

    print(
        f"\nRunning centrifuge for {sample}; directory: {directory}; threads: {threads}; clf_idx: {clf_idx}/{kingdom}"
    )
    print(fastq_folder, fastq4cent, os.path.isfile(fastq4cent))

    reads = 0
    classified_reads = 0
    if os.path.isdir(fastq_folder):
        in_folder_glob = glob.glob(fastq4cent)
        if len(in_folder_glob) > 0:
            in_folder_glob = in_folder_glob[0]
            read_count_glob = glob.glob(read_count)[0]
            reads = file_len(read_count_glob) / 4
            classifier_index = f"{clf_idx}/{kingdom}"
            cent_out = f"{directory}/analysis/sample_data/{sample}/centrifuge/"
            report = f"{cent_out}{sample}_{kingdom}{no_host}_centrifuge_report.tsv"
            troubleshooting = f"{cent_out}{sample}_{kingdom}{no_host}_centrifuge_troubleshooting_report.tsv"
            print(f"Saving to: {report}\n & {troubleshooting}\n")
            subprocess.run(
                [
                    "centrifuge",
                    "-q",
                    "-x",
                    classifier_index,
                    in_folder_glob,
                    "-p",
                    str(threads),
                    "--report-file",
                    report,
                    "-S",
                    troubleshooting,
                ],
                capture_output=True,
            )
            if os.path.isfile(report):
                report_df = pd.read_csv(report, delimiter="\t")
                classified_reads = report_df.numUniqueReads.sum()
                interpret_centrifuge(sample, troubleshooting, directory, kingdom)
                print(f"Centrifuge classified_reads count: {classified_reads}")
    print(f"reads, classified_reads: {reads}, {classified_reads}")
    Centrifuge_AA_list = [reads, 0]
    try:
        Centrifuge_AA_list = [int(reads), int(classified_reads)]
    except:
        Centrifuge_AA_list = [0, 0]

    print(f"{epoch_dir_host}/Centrifuge-AA-{sample}.json")
    with open(
        f"{epoch_dir_host}/Centrifuge-AA-{sample}.json", "w", encoding="utf-8"
    ) as f:
        json.dump(Centrifuge_AA_list, f)


def run_BLAST(
    directory: str,
    threads: int,
    sample_classifier: str,
    no_host: str,
    bN_name: str,
    hsbn_path_str: str,
    blastN_lib: str,
    epoch_dir_host: str,
    github: str,
    new_fqs: list,
) -> None:
    current_sample_classifier = sample_classifier.split("**")
    sample, classifier, kingdom = current_sample_classifier
    print(sample)
    print(classifier)

    new_fq = [s for s in new_fqs if sample in s]

    fastq = f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
    fasta = f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fasta"
    read_count = fastq

    if len(no_host) > 0:
        no_host_fq = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/no_host_{sample}.fastq"
        )
        fastq = check_file_exists(no_host_fq, fastq)
        fasta = fastq.replace("fastq", "fasta")
    hs_reads = 0
    hs_classified_reads = 0
    fastq_glob = glob.glob(fastq)

    if len(fastq_glob) > 0:
        fastq_glob = fastq_glob[0]
        read_count_glob = glob.glob(read_count)[0]
        hs_reads = file_len(read_count_glob) / 4
        print(hs_reads, sample)
        if len(new_fq) > 0:
            print(f"Host reads removed: {new_fq} ?")
            fastq = new_fq[0]
            fasta = fastq.replace("fastq", "fasta")

        gen_fasta(fasta, fastq_glob)
    fasta_glob = glob.glob(fasta)

    if len(fasta_glob):
        fasta_glob = fasta_glob[0]
        out_file = f"{directory}/analysis/sample_data/{sample}/blastN/trimmed_{sample}{no_host}_{bN_name}.hs_bn.tsv"
        cwd = os.getcwd()
        if os.path.isdir(f"{cwd}/hbndb"):
            shutil.rmtree(f"{cwd}/hbndb")
        print(
            f"Path to high speed blastn command: {hsbn_path_str}, \n'-outfmt', '6', '-task', 'megablast', '-evalue', '1e-50', '-out', \nSave location: {out_file}, \n'-max_target_seqs', '1', '-max_hsps', '1', '-num_threads', {threads}, \nInput fasta: {fasta_glob}, \nDatabase: {blastN_lib}\n"
        )
        subprocess.run(
            [
                hsbn_path_str,
                "-outfmt",
                "6",
                "-task",
                "megablast",
                "-evalue",
                "1e-50",
                "-out",
                out_file,
                "-max_target_seqs",
                "1",
                "-max_hsps",
                "1",
                "-num_threads",
                str(threads),
                fasta_glob,
                blastN_lib,
            ],
            capture_output=True,
        )
        hs_classified_reads = interpret_blastn(
            directory, sample, "hs_bn", False, bN_name, no_host, github
        )

    print(f"hs_reads: {hs_reads}; hs_classified_reads: {hs_classified_reads}, {sample}")

    BLAST_AA_list = [hs_reads, 0]
    try:
        BLAST_AA_list = [int(hs_reads), int(hs_classified_reads)]
    except:
        print(
            f"\n\nFAILED TO CONVERT NONETYPE TO INT FOR {sample}; hs_reads: {hs_reads}; hs_classified_reads: {hs_classified_reads}\n\n"
        )
        BLAST_AA_list = [0, 0]

    with open(f"{epoch_dir_host}/BLAST-AA-{sample}.json", "w", encoding="utf-8") as f:
        json.dump(BLAST_AA_list, f)


def clean_up_temp(path: str) -> None:
    try:
        shutil.rmtree(path)
    except Exception as e:
        print(f"Failed to delete: {path, e}")


def mp_function_ce_BN(
    data_directory: str,
    threads: int,
    no_host: str,
    bN_name: str,
    kingdom: str,
    new_fqs: str,
    hsbn_path_str: str,
    blastN_lib: str,
    clf_idx: str,
    epoch_dir_host: str,
    github: str,
    sample_classifier: str,
) -> None:
    with tempfile.TemporaryDirectory() as path:
        os.chdir(path)
        print(f"\n{path}")
        print(os.getcwd())
        print(sample_classifier)
        if "BLAST" in sample_classifier:
            run_BLAST(
                data_directory,
                threads,
                sample_classifier,
                no_host,
                bN_name,
                hsbn_path_str,
                blastN_lib,
                epoch_dir_host,
                github,
                new_fqs,
            )
            print(
                f"Finishing centrifuge and BLAST multiprocessing for {sample_classifier}"
            )
        if "Centrifuge" in sample_classifier:
            run_centrifuge(
                data_directory,
                sample_classifier,
                threads,
                clf_idx,
                bN_name,
                new_fqs,
                no_host,
                epoch_dir_host,
            )
            print(
                f"Finishing centrifuge and BLAST multiprocessing for {sample_classifier}"
            )
        clean_up_temp(path)


def main(
    data_directory: str,
    hsbn_path_str: str,
    clf_idx: str,
    bN_name: str,
    blastN_lib: str,
    processes_input: int,
    kingdom: str,
    no_host,
    new_fqs: list,
    epoch_dir_host: str,
    github: str,
    samples: list,
) -> None:
    metagenomic_classifiers = ["BLAST", "Centrifuge"]
    mp_variables = []
    for sample in samples:
        for metagenomic_classifier in metagenomic_classifiers:
            mp_variables.append(f"{sample}**{metagenomic_classifier}**{kingdom}")
    print(f"Number of samples to classify: {len(mp_variables)}")
    no_processes, threads = threads_n_processes(processes_input)
    print(no_processes, threads)

    # CENTRIFUGE AND BLAST
    partial_func = partial(
        mp_function_ce_BN,
        data_directory,
        threads,
        no_host,
        bN_name,
        kingdom,
        new_fqs,
        hsbn_path_str,
        blastN_lib,
        clf_idx,
        epoch_dir_host,
        github,
    )
    with Pool(processes=no_processes) as p:
        p.map(partial_func, mp_variables)


if __name__ == "__main__":

    blastN_libraries = database_config.blastN_databases()
    parser = argparse.ArgumentParser(
        description="Running metagenomic & sample multiprocessing"
    )
    # choose a blastN library
    parser.add_argument(
        "-b", "--BLASTn_lib", choices=blastN_libraries.keys(), default="cviral"
    )
    args = parser.parse_args()

    multiple_fastq = "/home/james/SMART-CAMP/configs/viral_DNA_all2.txt"
    print("Finding samples for processing.")
    samples = [line.rstrip("\n").split(",") for line in open(f"{multiple_fastq}")]
    samples = [item for sublist in samples for item in sublist]
    print(samples)

    new_fqs = [
        "/mnt/usersData/Viral_CHO/analysis/sample_data/20211210_ssDNA_MVM1-10CHOK1-Concentrate_500000CFU_23_36/trimmed/no_host_20211210_ssDNA_MVM1-10CHOK1-Concentrate_500000CFU_23_36.fastq"
        "/mnt/usersData/Viral_CHO/analysis/sample_data/20211210_ssDNA_MVM1-10CHOK1-Filtrate_500000CFU_23_36/trimmed/no_host_20211210_ssDNA_MVM1-10CHOK1-Filtrate_500000CFU_23_36.fastq"
        "/mnt/usersData/Viral_CHO/analysis/sample_data/20220110_ssDNA_MVM-10E6-1st2ndStrandSynth-zymo-no-synth_5000000CFU_24_36/trimmed/trimmed_20220110_ssDNA_MVM-10E6-1st2ndStrandSynth-zymo-no-synth_5000000CFU_24_36.fastq"
        "/mnt/usersData/Viral_CHO/analysis/sample_data/20220110_ssDNA_MVM-10E6-1st2ndStrandSynth-zymo-synth_5000000CFU_24_36/trimmed/trimmed_20220110_ssDNA_MVM-10E6-1st2ndStrandSynth-zymo-synth_5000000CFU_24_36.fastq"
    ]

    new_fqs = []

    data_directory = "/mnt/usersData/Viral_CHO/"
    hsbn_path_str = "/home/james/apps/hs-blastn"
    github = "/home/james/SMART-CAMP/"
    clf_idx = "/home/james/SequencingData/Centrifuge_libraries/"
    bN_name = args.BLASTn_lib
    blastN_lib = blastN_libraries[args.BLASTn_lib]
    kingdom = ""
    if (
        bN_name == "bacteria"
        or bN_name == "16S_16S"
        or bN_name == "16S_23S"
        or bN_name == "filter_bacteria"
    ):
        kingdom = "bacteria"
    if bN_name == "virus" or bN_name == "cviral":
        kingdom = "virus"
    if bN_name == "fungal_all":
        kingdom = "fungus"
    if bN_name == "v_f_b":
        kingdom = "v_f_b"

    no_host = "_no_host"
    processes_input = 10

    epoch_time = int(time.time())
    epoch_dir = f"{data_directory}/analysis/run/{str(epoch_time)}"
    epoch_dir_host = f"{epoch_dir}/host/"

    main(
        data_directory,
        hsbn_path_str,
        clf_idx,
        bN_name,
        blastN_lib,
        processes_input,
        kingdom,
        no_host,
        new_fqs,
        epoch_dir_host,
        github,
        samples,
    )


# time /home/james/SMART-CAMP/flat_sample_classifier_test.py -b "cviral"
# time /home/james/SMART-CAMP/flat_sample_classifier_test.py -b "filter_bacteria"
# time /home/james/SMART-CAMP/flat_sample_classifier_test.py -b "fungal_all"
