#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 11:42:01 2021

@author: mangi
"""

from pathlib import Path
import os


# function that outputs dictionary of coverage databases, this can be added to if you custom databases are required.
def coverage_databases() -> dict:
    # # # 16S SEQUENCES
    # seq_16S_16S = f"{str(Path.home())}/SequencingData/Centrifuge_libraries/bacteria/library/16S_16S/"

    # # # 16S-23S SEQUENCES
    # seq_16S_23S = f"{str(Path.home())}/SequencingData/Centrifuge_libraries/bacteria/library/16S_23S/"

    # # SMALL SMART full bacterial sequence -> filtered for only samples not described as "plasmid" or "chromosome"; though in this context chromosome simply means full genome
    # small_smart_full_bacteria_seq = f"{str(Path.home())}/SequencingData/Centrifuge_libraries/bacteria/library/bacteria_small_smart/"

    # # full bacterial sequence
    # full_bacteria_seq = f"{str(Path.home())}/SequencingData/Centrifuge_libraries/bacteria/library/bacteria/"

    # # refseq viruses with bacteriophages
    # virus_seq = f"{str(Path.home())}/SequencingData/Centrifuge_libraries/viral/library/refseq_viral/"

    # # clustered viral genomes RVDB
    # cl_viral_seq = (
    #     f"{str(Path.home())}/SequencingData/Centrifuge_libraries/viral/library/C_RVDB/"
    # )

    # human sequence
    chinese_hamster_seq = "/mnt/usersData/krakenDB/vertebrate/vertebrate_mammalian/Chromosome/chinese_hamster/"

    # full human sequence
    full_human_seq = f"{str(Path.home())}/SequencingData/Centrifuge_libraries/human/"

    # # fungal genomes NCBI
    # seq_fungal = f"{str(Path.home())}/SequencingData/NCBI_fungal_full/genomic/"

    # print("Checking if coverage folders exist:")
    # print(f"16S: {os.path.isdir(seq_16S_16S)}")
    # print(f"16S-23S: {os.path.isdir(seq_16S_23S)}")
    # print(f"Full bacterial: {os.path.isdir(full_bacteria_seq)}")
    # print(
    #     f"Full bacterial (small non-random): {os.path.isdir(small_smart_full_bacteria_seq)}"
    # )
    # print(f"Full human: {os.path.isdir(full_human_seq)}")
    # print(f"clustered viral: {os.path.isdir(cl_viral_seq)}")
    # print(f"NCBI viral: {os.path.isdir(virus_seq)}")
    # print(f"full human_seq: {os.path.isdir(full_human_seq)}")
    # print(f"chinese_hamster_seq: {os.path.isdir(chinese_hamster_seq)}")
    # print(f"All fungal: {os.path.isdir(seq_fungal)}\n")

    coverage_libraries = {
        # "16S_16S": seq_16S_16S,
        # "16S_23S": seq_16S_23S,
        # "bacteria": full_bacteria_seq,
        # "filter_bacteria": small_smart_full_bacteria_seq,
        "human": full_human_seq,
        "chinese_hamster": chinese_hamster_seq,
        # "virus": virus_seq,
        # "cviral": cl_viral_seq,
        # "fungal_all": seq_fungal,
    }
    return coverage_libraries


def blastN_databases() -> dict:
    # viral-bacterial-fungal SEQUENCES (uses NCBI virus, fungal and smart small bacteria)
    v_f_b = (
        f"{str(Path.home())}/SequencingData/Centrifuge_libraries/v_f_b/v-f-b-small.fna"
    )

    # viral-bacterial-fungal SEQUENCES (uses NCBI virus, fungal and smart small bacteria)
    # includes clostridium sporogens and candida albicans 
    v_f_b2 = (
        f"{str(Path.home())}/SequencingData/Centrifuge_libraries/v_f_b/v-f-b-small2.fna"
    )
    
    # # 16S_23S SEQUENCES
    # seq_16S_16S = (
    #     f"{str(Path.home())}/SequencingData/NCBI_bacteria_16S_16S/fasta/16S_16S_db.fna"
    # )

    # # # 16S_23S SEQUENCES
    # seq_16S_23S = (
    #     f"{str(Path.home())}/SequencingData/NCBI_bacteria_16S_23S/fasta/16S_23S_db.fna"
    # )

    # # SMALL SMART full bacterial sequence -> filtered for only samples not described as "plasmid" or "chromosome"; though in this context chromosome simply means full genome
    # small_smart_full_bacteria_seq = f"{str(Path.home())}/SequencingData/NCBI_bacteria_full/fasta/small_smart_full_bacteria.fasta"

    # # full bacteria
    # full_bacteria_seq = (
    #     f"{str(Path.home())}/SequencingData/NCBI_bacteria_full/fasta/bact_full_db.fna"
    # )

    # # clustered viral
    # cl_viral_seq = (
    #     f"{str(Path.home())}/SequencingData/NCBI_RVDB/virus_RVDBs/C-RVDBv21.0.fasta"
    # )

    # # viral genomes NCBI
    # virus_seq = "/home/james/SequencingData/Centrifuge_libraries/viral/viral_input-sequences.fna"

    # CH sequence
    chinese_hamster_seq = "/mnt/usersData/krakenDB/vertebrate/vertebrate_mammalian/Chromosome/chinese_hamster/GCF_000223135.1_CHOK1_CriGri_1.0_genomic.fna"

    # human sequence
    full_human_seq = (
        f"{str(Path.home())}/SequencingData/NCBI_human/fasta/GRCh38_latest_genomic.fna"
    )

    # # fungal genomes NCBI
    # seq_fungal = (
    #     f"{str(Path.home())}/SequencingData/NCBI_fungal_full/fungal_input-sequences.fna"
    # )

    print("Checking if high speed blastN files exist:")
    print(f"v_f_b: {os.path.isfile(v_f_b)}")
    print(f"v_f_b2: {os.path.isfile(v_f_b2)}")    
    # print(f"16S: {os.path.isfile(seq_16S_16S)}")
    # print(f"16S-23S: {os.path.isfile(seq_16S_23S)}")
    # print(f"Bacterial all: {os.path.isfile(full_bacteria_seq)}")
    # print(
    #     f"Bacterial (small non-random) all: {os.path.isfile(small_smart_full_bacteria_seq)}"
    # )
    # print(f"RVDB clustered all: {os.path.isfile(cl_viral_seq)}")
    # print(f"NCBI viral: {os.path.isfile(virus_seq)}")
    print(f"Human all: {os.path.isfile(full_human_seq)}")
    print(f"CHO: {os.path.isfile(chinese_hamster_seq)}")
    # print(f"Fungal all: {os.path.isfile(seq_fungal)}\n")

    blastN_libraries = {
        "v_f_b": v_f_b,
        "v_f_b2": v_f_b2,
    #     "16S_16S": seq_16S_16S,
    #     "16S_23S": seq_16S_23S,
    #     "bacteria": full_bacteria_seq,
    #     "filter_bacteria": small_smart_full_bacteria_seq,
        "human": full_human_seq,
        "chinese_hamster": chinese_hamster_seq,
    #     "cviral": cl_viral_seq,
    #     "virus": virus_seq,
    #     "fungal_all": seq_fungal,
    }
    return blastN_libraries


# function that outputs dictionary of coverage databases, this can be added to if you custom databases are required.
def host_coverage_databases() -> dict:

    # CHO clean up
    # find ONLY the chromosome headers
    # rg -i ">.+chromosome" Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna
    # find all unplaced contigs
    # rg -i ">.+unplaced" Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna

    # delete all lines after given line - these are unassigned scaffolds
    # sed '28715648,$d' Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna > CHO_chromosomes_only.fna

    # delete lines after line number -- chromosome 1
    # sed '3446229,$d' CHO_chromosomes_only.fna > CHO_chromosomes_only_1A.fna
    # delete lines before line number -- all other chromosomes
    # sed '6876127,$!d' CHO_chromosomes_only.fna > CHO_chromosomes_only_1B.fna

    chinese_hamster_chromosomes = "/mnt/usersData/krakenDB/library/vertebrate_mammalian/Chromosome/Cricetulus_griseus_strain_17A_GY-tax10029-GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"

    # full human sequence
    full_human_seq = "/home/james/SequencingData/Centrifuge_libraries/human/human_input-sequences.fna"

    print("Checking if coverage folders exist:")
    print(f"Full human: {Path(full_human_seq).is_file()}")
    print(
        f"Chinese_hamster_seq (chromosomes): {Path(chinese_hamster_chromosomes).is_file()}"
    )

    host_coverage_libraries = {
        "human": full_human_seq,
        "chinese_hamster_chromosomes": chinese_hamster_chromosomes,
    }
    return host_coverage_libraries
