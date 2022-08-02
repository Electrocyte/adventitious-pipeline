#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 14:26:07 2021

@author: mangi
"""

import re
import os
from pathlib import Path
import pandas as pd
import glob
import json
import subprocess
from multiprocessing import Pool
from functools import partial

encoding = "utf-8"


class EmptyDirecory(Exception):
    pass


# BLAST libraries
def host_databases(host_dbs) -> dict:
    blastN_libraries = {
        "hsbn_human": host_dbs["human"],
        "hsbn_CHOK1": host_dbs["chinese_hamster"],
    }

    return blastN_libraries


def gen_fasta(fasta: str, in_folder_glob: str) -> None:
    if not Path(fasta).is_file():
        fastq2fasta = subprocess.run(
            ["sed", "/^@/!d;s//>/;N", in_folder_glob], capture_output=True
        )
        returncode = fastq2fasta.returncode
        if returncode != 0:
            print("Unable to generate fasta from fastq")
        with open(fasta, "w") as f:
            print(f"Saving to: {fasta}")
            f.write(str(fastq2fasta.stdout, encoding))


def human_ref_handler(row: str) -> str:
    compiled = re.compile(r"\|(\w+.\w+)")
    return compiled.findall(row)[0]


def split_name(df: pd.DataFrame) -> pd.DataFrame:
    df_split = pd.DataFrame(
        df.name.str.split(" ", 2).tolist(), columns=["genus", "species", "strain"]
    )
    df_split["name"] = df_split[["genus", "species", "strain"]].apply(
        lambda row: " ".join(row.values.astype(str)), axis=1
    )
    return df_split


# collect species genus, species and strain
def import_deanonymised_nomenclature(directory: str) -> pd.DataFrame:
    seqid_deanonymiser = pd.read_csv(
        directory, sep=",", header=None, names=["sseqid", "name"]
    )
    seqid_deanonymiser = seqid_deanonymiser.iloc[1:]
    seqid_deanonymised_names = split_name(seqid_deanonymiser)
    seqid_deanonymised_clean = seqid_deanonymiser.merge(
        seqid_deanonymised_names, on="name", how="inner"
    )
    seqid_deanonymised_clean = seqid_deanonymised_clean.drop_duplicates(
        subset=["sseqid"]
    )
    return seqid_deanonymised_clean


def blastn_import(directory: str, seqid_dir: str) -> pd.DataFrame:
    sseqids = import_deanonymised_nomenclature(seqid_dir)
    print(directory, seqid_dir)
    df = pd.read_csv(
        f"{directory}",
        sep="\t",
        header=None,
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
    print(df)
    if not df.empty:
        if df["sseqid"].str.contains("ref")[0]:
            df["sseqid"] = df["sseqid"].apply(human_ref_handler)
        seqid_deanonymised = df.merge(sseqids, on="sseqid", how="inner")
        return seqid_deanonymised


def fastq_formatter(new_fq: str, df: pd.DataFrame) -> None:
    lines = []
    for _index, row in df.iterrows():
        lines.append(row["readID"])
        lines.append(row["Sequence"])
        lines.append(row["Plus"])
        lines.append(row["Quality"])
    if not os.path.isfile(new_fq):
        print(f"\n\nSaving...{new_fq}")
        with open(new_fq, "w") as f:
            f.write("\n".join(lines))


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
    no_human_df = clean_df[~clean_df["readID_D"].isin(read_ids)]
    no_human_df = no_human_df.drop(["readID_D"], axis=1)
    return no_human_df


# import trimmed fastq files
def import_trimmed_reads(directory: str, sample: str, Direct_RNA: bool) -> tuple:
    folder = f"{directory}/analysis/sample_data/{sample}/trimmed/"
    reads = f"{folder}/trimmed*.fastq"
    if Direct_RNA:
        reads = f"{folder}/U2T_*.fastq"
    print(reads)
    fq_reads = glob.glob(f"{reads}", recursive=True)[0]
    return fq_reads


def import_centrifuge_data(directory):
    # check if file is empty or not
    if os.stat(directory).st_size > 0:
        df = pd.read_csv(directory, sep="\t")
        return df


# count number of lines in a file FAST (4X vs pandas)
# must subtract header line (-1)
# I would add i=-1 before for loop, since this code doesn't work for empty files.
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# clf_func = partial(run_centrifuge_n_blastn_host, directory, sample, threads_, cent_index, blastN_lib, bN_name, Direct_RNA, hsbn_path_str, epoch_dir, host)
def run_centrifuge_n_blastn_n_krakenuniq_host(
    directory: str,
    sample: str,
    threads: str,
    clf_idx: str,
    blastN_lib: str,
    bN_name: str,
    Direct_RNA: bool,
    hsbn_path_str: str,
    epoch_dir: str,
    cent_host: str,
    krakenuniq_command: str,
    classifier: str,
) -> None:

    if classifier == "BLAST":
        fastq = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
        )
        fasta = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fasta"
        )
        out_file = f"{directory}/analysis/sample_data/{sample}/blastN/{sample}_{bN_name}_host_removed.hsbn.tsv"
        if Direct_RNA:
            fastq = (
                f"{directory}/analysis/sample_data/{sample}/trimmed/U2T_{sample}.fastq"
            )
            fasta = (
                f"{directory}/analysis/sample_data/{sample}/trimmed/U2T_{sample}.fasta"
            )
            out_file = f"{directory}/analysis/sample_data/{sample}/blastN/U2T_{sample}_{bN_name}_host_removed.hsbn.tsv"

        host_read_count = 0
        if os.path.isfile(out_file):
            if Path(out_file).stat().st_size > 0:
                # subtract one for the header
                host_read_count = file_len(out_file) - 1

        print(fasta)
        print(f"Metagenomic analysis, blastN using {blastN_lib} library")
        fastq_glob = glob.glob(fastq)[0]

        gen_fasta(fasta, fastq_glob)
        fasta_glob = glob.glob(fasta)[0]

        print(f"Saving to: {out_file}")
        if not Path(out_file).is_file():
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
                    threads,
                    fasta_glob,
                    blastN_lib,
                ],
                capture_output=True,
            )

            if Path(out_file).stat().st_size > 0:
                # subtract one for the header
                host_read_count = file_len(out_file) - 1

        # save sample name for later checking /run/utc/host/
        BLAST_host_list = [out_file, str(host_read_count)]
        with open(f"{epoch_dir}/BLAST-{sample}.json", "w", encoding="utf-8") as f:
            json.dump(BLAST_host_list, f)

    if classifier == "Centrifuge":
        print(
            f"\nRunning centrifuge for {sample}; directory: {directory}; threads: {threads}; clf_idx: {clf_idx}"
        )

        fastq4cent = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
        )
        fastq_folder = f"{directory}/analysis/sample_data/{sample}/trimmed/"
        if Direct_RNA:
            fastq4cent = (
                f"{directory}/analysis/sample_data/{sample}/trimmed/U2T_{sample}.fastq"
            )

        if os.path.isdir(fastq_folder):
            in_folder_glob = glob.glob(fastq4cent)[0]

            # the host removal files for centrifuge need to be named in the following format; centrifuge used a prefix to identify the indices
            classifier_index = f"{clf_idx}/{cent_host}/{cent_host}"
            print(f"Classifier index: {classifier_index}")
            cent_out = f"{directory}/analysis/sample_data/{sample}/centrifuge/"
            report = f"{cent_out}{sample}_{cent_host}_centrifuge_report.tsv"
            troubleshooting = (
                f"{cent_out}{sample}_{cent_host}_centrifuge_troubleshooting_report.tsv"
            )

            host_read_count = 0

            print(f"Saving to: {report}\n & {troubleshooting}\n")
            if not os.path.isfile(troubleshooting):
                centrifuge_run = subprocess.run(
                    [
                        "centrifuge",
                        "-q",
                        "-x",
                        classifier_index,
                        in_folder_glob,
                        "-p",
                        threads,
                        "--report-file",
                        report,
                        "-S",
                        troubleshooting,
                    ],
                    capture_output=True,
                )
                returncode = centrifuge_run.returncode

                if returncode != 0:
                    print("Centrifuge analysis did not generate any classifications")

            if Path(report).stat().st_size > 0:
                ce_report = pd.read_csv(report, delimiter="\t")
                if not ce_report.empty:
                    host_read_count = ce_report.iloc[0]["numUniqueReads"]

            # save sample name for later checking /run/utc/host/
            centrifuge_host_list = [troubleshooting, str(host_read_count), report]
            with open(
                f"{epoch_dir}/Centrifuge-{sample}.json", "w", encoding="utf-8"
            ) as f:
                json.dump(centrifuge_host_list, f)

    if classifier == "Krakenuniq":

        host_read_count = 0
        full_fastq = (
            f"{directory}/analysis/sample_data/{sample}/trimmed/trimmed_{sample}.fastq"
        )
        if Direct_RNA:
            full_fastq = (
                f"{directory}/analysis/sample_data/{sample}/trimmed/U2T_{sample}.fastq"
            )

        database = "/mnt/usersData/krakenDB/vertebrate/"
        print(
            f"\nRunning krakenuniq for {sample}; directory: {directory}; threads: {threads}; database: {database}"
        )
        report_file = f"{directory}/analysis/sample_data/{sample}/kraken/{cent_host}_kraken_report.tsv"
        print(f"Saving to: {report_file}")
        output = f"{directory}/analysis/sample_data/{sample}/kraken/{cent_host}_kraken_classification_report.tsv"

        if os.path.isdir(database):
            if not os.path.isfile(report_file):
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
        else:
            print("Skipping kraken analysis.")

        if os.path.isfile(output):
            if Path(output).stat().st_size > 0:
                kr_report = pd.read_csv(
                    output,
                    delimiter="\t",
                    names=["classified", "readID", "taxID", "seq_len", "LCA-map"],
                )
                if not kr_report.empty:
                    host_read_count = len(
                        kr_report.loc[
                            (kr_report["taxID"] == 10029) | (kr_report["taxID"] == 9606)
                        ]
                    )

        # save sample name for later checking /run/utc/host/
        kr_host_list = [output, str(host_read_count)]
        with open(f"{epoch_dir}/Krakenuniq-{sample}.json", "w", encoding="utf-8") as f:
            json.dump(kr_host_list, f)


# Only works for U2T RNA samples
# new_fq, ce_host_read_count, hsbn_host_read_count = mp_host_remove.remove_host_reads(sample, top_dir, threads_, cent_index, Direct_RNA, hsbn_path_str, host, epoch_dir, github)
def remove_host_reads(
    sample: str,
    directory: str,
    threads_: str,
    cent_index: str,
    Direct_RNA: bool,
    hsbn_path_str: str,
    host: str,
    epoch_dir: str,
    github: str,
    krakenuniq_command: str,
    blastN_libraries: dict,
) -> (str, int, int, int):
    host_spp = host_databases(blastN_libraries)
    new_fq = ""
    hsbn_host_read_count = 0
    ce_host_read_count = 0
    kr_host_read_count = 0
    for k, v in host_spp.items():
        if host in k:
            blastN_lib = v
            bN_name = k
    print(host)
    if Path(f"{directory}/analysis/sample_data/{sample}/trimmed/").is_dir():
        if len(os.listdir(f"{directory}/analysis/sample_data/{sample}/trimmed/")) == 0:
            print(
                f"Directory ({directory}/{sample}/trimmed/) is empty; skipping sample."
            )
        else:
            print(
                f"Is folder a directory? {Path(f'{directory}/analysis/sample_data/{sample}').is_dir()}"
            )

            seqid_dir = f"{github}/host_seqids.csv"
            print(f"Human seqids: {seqid_dir}")

            new_fq = f"{directory}/analysis/sample_data/{sample}/trimmed/no_host_{sample}.fastq"
            if Direct_RNA:
                new_fq = f"{directory}/analysis/sample_data/{sample}/trimmed/U2T_no_host_{sample}.fastq"

            output = f"{directory}/analysis/sample_data/{sample}/kraken/{host}_kraken_classification_report.tsv"
            if os.path.isfile(output):
                if os.path.isfile(new_fq):

                    # hsblastn
                    out_file = f"{directory}/analysis/sample_data/{sample}/blastN/{sample}_{bN_name}_host_removed.hsbn.tsv"
                    if Direct_RNA:
                        out_file = f"{directory}/analysis/sample_data/{sample}/blastN/U2T_{sample}_{bN_name}_host_removed.hsbn.tsv"
                    if Path(out_file).stat().st_size > 0:
                        # subtract one for the header
                        # hsbn_host_read_count = file_len(out_file) - 1
                        hs_report = pd.read_csv(out_file, delimiter="\t")
                        if not hs_report.empty:
                            hsbn_host_read_count = len(hs_report)

                    # centrifuge
                    cent_out = f"{directory}/analysis/sample_data/{sample}/centrifuge/"
                    report = f"{cent_out}{sample}_{host}_centrifuge_report.tsv"
                    if Path(report).stat().st_size > 0:
                        # subtract one for the header
                        # ce_host_read_count = file_len(report) - 1
                        ce_report = pd.read_csv(report, delimiter="\t")
                        if not ce_report.empty:
                            ce_host_read_count = ce_report.iloc[0]["numUniqueReads"]

                    # krakenuniq
                    report_file = f"{directory}/analysis/sample_data/{sample}/kraken/{host}_kraken_report.tsv"
                    output = f"{directory}/analysis/sample_data/{sample}/kraken/{host}_kraken_classification_report.tsv"
                    if Path(output).stat().st_size > 0:
                        kr_report = pd.read_csv(
                            output,
                            delimiter="\t",
                            names=[
                                "classified",
                                "readID",
                                "taxID",
                                "seq_len",
                                "LCA-map",
                            ],
                        )
                        if not kr_report.empty:
                            kr_host_read_count = len(
                                kr_report.loc[
                                    (kr_report["taxID"] == 10029)
                                    | (kr_report["taxID"] == 9606)
                                ]
                            )

                    return (
                        new_fq,
                        ce_host_read_count,
                        hsbn_host_read_count,
                        kr_host_read_count,
                    )

            metagenomic_classifiers = ["BLAST", "Centrifuge", "Krakenuniq"]
            clf_func = partial(
                run_centrifuge_n_blastn_n_krakenuniq_host,
                directory,
                sample,
                threads_,
                cent_index,
                blastN_lib,
                bN_name,
                Direct_RNA,
                hsbn_path_str,
                epoch_dir,
                host,
                krakenuniq_command,
            )
            # concat all those chunked files.
            print("Beginning to concatenate chunked files")
            with Pool(processes=3) as p:
                p.map(clf_func, metagenomic_classifiers)

            # will require use of epoch_dir to generate data that can be later found.
            # import BLAST data
            with open(f"{epoch_dir}/BLAST-{sample}.json") as json_b_data:
                BLAST_json = json.loads(json_b_data.read())
            blastN_out, hsbn_host_read_count = BLAST_json  # (str, str)

            # import centrifuge data
            with open(f"{epoch_dir}/Centrifuge-{sample}.json") as json_c_data:
                Centrifuge_json = json.loads(json_c_data.read())
            (
                troubleshooting_file,
                ce_host_read_count,
                report,
            ) = Centrifuge_json  # (str, str, str)

            # import krakenuniq data
            with open(f"{epoch_dir}/Krakenuniq-{sample}.json") as json_k_data:
                kr_json = json.loads(json_k_data.read())
            kr_reads, kr_host_read_count = kr_json  # (str, str)

            cent_df = import_centrifuge_data(troubleshooting_file)
            kr_df = pd.DataFrame()
            if os.path.isfile(kr_reads):
                if Path(kr_reads).stat().st_size > 0:
                    kr_df = pd.read_csv(
                        kr_reads,
                        delimiter="\t",
                        names=["classified", "readID", "taxID", "seq_len", "LCA-map"],
                    )

            if host in "CHO":
                hostid = 10029
            else:
                hostid = 9606

            # krakenuniq
            kr_read_id_list = []
            unions = []
            if not kr_df.empty:
                classified_only = kr_df.loc[
                    (kr_df["taxID"] == 10029) | (kr_df["taxID"] == 9606)
                ]
                kr_read_id_list = classified_only["readID"].tolist()
                unions = kr_read_id_list

            # centrifuge
            human_only = cent_df[cent_df["taxID"] == hostid]
            cent_read_id_list = human_only["readID"].tolist()

            # load in readID db
            blastn_df = blastn_import(blastN_out, seqid_dir)

            if blastn_df is not None:
                blastn_df = blastn_df.rename(columns={"qseqid": "readID"})
                bN_read_id_list = blastn_df["readID"].tolist()
                unions = list(
                    set(bN_read_id_list)
                    .union(set(cent_read_id_list))
                    .union(set(kr_read_id_list))
                )

            print(f"directory: {directory}, sample: {sample}, Direct_RNA: {Direct_RNA}")
            all_reads = import_trimmed_reads(directory, sample, Direct_RNA)
            if not os.stat(all_reads).st_size == 0:
                out_df = generate_fastq_df(all_reads, unions)
                if not out_df.empty:
                    fastq_formatter(new_fq, out_df)

                    print(ce_host_read_count, hsbn_host_read_count, kr_host_read_count)
                    return (
                        new_fq,
                        ce_host_read_count,
                        hsbn_host_read_count,
                        kr_host_read_count,
                    )
    return new_fq, ce_host_read_count, hsbn_host_read_count, kr_host_read_count


if __name__ == "__main__":
    top_dir = "/home/james/SequencingData/Direct_RNA/"
    sample = "20210323_RNA_PAO1unfTCPA_38700000CFU_1st_2h"
    threads = "20"
    Direct_RNA = True
    cent_index = "~/SequencingData/Centrifuge_libraries/"
    hsbn_path_str = "/home/james/apps/hs-blastn"
    # host = "human"
    epoch_dir = "15635648416"
