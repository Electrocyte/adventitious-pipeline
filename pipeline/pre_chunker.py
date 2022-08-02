#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 13:49:08 2021

@author: mangi
"""
# if merge file already exists, you can use bash:
# split --additional-suffix=.fastq -l 2000000 /mnt/usersData/Viral_CHO_chunk/merged.fastq out/chunk_
# command for chunking fastq files and adding both the prefix and suffix for fastq
# merged fastq file will be split by 2,000,000 lines

# generating multiple files of < 1GB size is more efficient when starting from individual
# small files, however it would be much simpler to complete this on python
# than bash, as it requires a step to check the size of the globbed files
# before merging them together

import glob
import os
from pathlib import Path
from multiprocessing import Pool
from functools import partial


def save_function(chunk_file: str, files_to_save: list) -> None:
    with open(chunk_file, "w") as outfile:
        for file2save in files_to_save:
            with open(file2save, "r") as infile:
                outfile.write(infile.read())


def chunk_fq_files(directory_in: str, directory_out, filename: str, save: bool) -> None:
    print(f"{directory_in}/{filename}")
    all_available_files = glob.glob(f"{directory_in}/{filename}")

    max_file = (
        Path(max(all_available_files, key=lambda x: os.stat(x).st_size)).stat().st_size
    )
    min_file = (
        Path(min(all_available_files, key=lambda x: os.stat(x).st_size)).stat().st_size
    )

    filesize_aggregate = 0
    file_number = 0
    chunk_number = 0
    files_to_save = []
    for available_file in all_available_files:

        if Path(available_file).is_file():
            filesize = Path(available_file).stat().st_size
            filesize_aggregate += filesize
            file_number += 1
            files_to_save.append(available_file)

            # handle situations where there are enough files to merge to give 1GB file
            if filesize_aggregate > 1e9:
                chunk_number += 1
                chunk_file = f"{directory_out}/chunk{chunk_number}.fastq"
                print(
                    f"\nWriting chunked file to memory and resetting filesize for next chunk. Number of files merged: {file_number}. Final file size: {'{:.2e}'.format(filesize_aggregate)}"
                )
                print(f"Saving fastq chunk to: {chunk_file}")
                if save:
                    save_function(chunk_file, files_to_save)

                # reset input parameters
                filesize_aggregate = 0
                file_number = 0
                files_to_save = []

            # handle situations where there are left over files from a sample that doesn't quite reach 1GB
            if chunk_number > 0 and available_file == all_available_files[-1]:
                chunk_number += 1
                chunk_file = f"{directory_out}/chunk{chunk_number}.fastq"
                print(
                    f"\nPreparing last small chunk. Number of files merged: {file_number}. Final file size: {'{:.2e}'.format(filesize_aggregate)}"
                )
                print(f"Saving fastq chunk to: {chunk_file}")
                if save:
                    save_function(chunk_file, files_to_save)

            # handle situations where there are too few files to reach 1GB at all
            if (
                chunk_number == 0
                and available_file == all_available_files[-1]
                and max_file != min_file
            ):
                chunk_file = f"{directory_out}/chunk{chunk_number}.fastq"
                print(
                    f"\nFewer than 1GB total files, saving as one file. Number of files merged: {file_number}. Final file size: {'{:.2e}'.format(filesize_aggregate)}"
                )
                print(f"Saving lonely fastq chunk to: {chunk_file}")
                if save:
                    save_function(chunk_file, files_to_save)

            # handle situations where there is only one file
            if len(all_available_files) == 1 and max_file == min_file:
                chunk_file = f"{directory_out}/chunk{chunk_number}.fastq"
                print(
                    f"\nOne file found, renaming. Number of files merged: {file_number}. Final file size: {'{:.2e}'.format(filesize_aggregate)}"
                )
                print(f"Saving fastq chunk to: {chunk_file}")
                if save:
                    save_function(chunk_file, files_to_save)


def test_all_inputs(
    test_many_fq: str,
    test_few_fq: str,
    test_one_fq: str,
    filename: str,
    directory_out: str,
    save: bool,
) -> None:
    print("\n\n>>Testing large number of fq files<<")
    chunk_fq_files(test_many_fq, directory_out, filename, save)
    print("\n\n>>Testing few fq files<<")
    chunk_fq_files(test_few_fq, directory_out, filename, save)
    print("\n\n>>Testing one fq file<<")
    chunk_fq_files(test_one_fq, directory_out, filename, save)


def multiprocess_sample(sample: str) -> None:
    filename = "F*.fastq"
    directory_in = f"/mnt/usersData/Viral_CHO/{sample}/**/fastq_pass/"
    directory_out = f"/mnt/usersData/Viral_CHO/analysis/{sample}/trimmed/"
    os.makedirs(directory_out, exist_ok=True)
    chunk_fq_files(directory_in, directory_out, filename, True)


if __name__ == "__main__":

    testing = False
    single_processing = True
    multi_processing = False
    filename = "F*.fastq"

    if testing:
        directory_test3 = "/mnt/usersData/Viral_CHO_chunk/20210825_ssDNA_MVM100-1CHOK1-QiagenNoDNaseNo95C_500000000CFU_15_36/20210825_1717_X1_FAQ74159_fa4e6eea/fastq_pass3/"
        directory_test2 = "/mnt/usersData/Viral_CHO_chunk/20210825_ssDNA_MVM100-1CHOK1-QiagenNoDNaseNo95C_500000000CFU_15_36/20210825_1717_X1_FAQ74159_fa4e6eea/fastq_pass2/"
        directory_test = "/mnt/usersData/Viral_CHO_chunk/20210825_ssDNA_MVM100-1CHOK1-QiagenNoDNaseNo95C_500000000CFU_15_36/20210825_1717_X1_FAQ74159_fa4e6eea/fastq_pass/"
        directory_out = "/mnt/usersData/Viral_CHO_chunk/analysis/20210825_ssDNA_MVM100-1CHOK1-QiagenNoDNaseNo95C_500000000CFU_15_36/trimmed/"

        test_all_inputs(
            directory_test,
            directory_test2,
            directory_test3,
            filename,
            directory_out,
            False,
        )

        print("\n\n>>Running sample...<<")
        chunk_fq_files(directory_test, directory_out, filename, False)

    if not testing:
        # samples = ["20210420_DNA_CHOK1_1000000CFU_3_42", "20210420_ssDNA_MVM_10000000000CFU_3_42", \
        #            "20210420_ssDNA_MVM100-1CHOK1_5000000000CFU_3_42", "20210420_ssDNA_MVM500-1CHOK1_25000000000CFU_3_42", \
        #                "20210420_ssDNA_MVM50-1CHOK1_2500000000CFU_3_42"]
        # samples = ['20210902_ssDNA_MVM100-1CHOK1-Filt-DNase-NoRapid_500000000CFU_16_36', '20210902_ssDNA_MVM100-1CHOK1-Filt-DNase-Rapid_500000000CFU_16_36', \
        #            '20210902_ssDNA_MVM100-1CHOK1-Filt-NoDNase-NoRapid_500000000CFU_16_36', '20210902_ssDNA_MVM100-1CHOK1-NoFilt-DNase-Rapid_500000000CFU_16_36', '20210902_ssDNA_MVM100-1CHOK1-NoFilt-NoDNase-NoRapid_500000000CFU_16_36']
        samples = ["S8B5"]

        # SINGLE PROCESSING
        if single_processing:
            for sample in samples:
                directory_in = f"/mnt/usersData/Fungal_test/{sample}/20210906_0659_MN34889_FAP20808_79383d5d/fastq_pass/"
                directory_out = f"/mnt/usersData/Fungal_test/{sample}/20210906_0659_MN34889_FAP20808_79383d5d/fastq_processed/"
                # directory_in = f"/mnt/usersData/Viral_CHO/{sample}/**/fastq_pass/"
                # directory_out = f"/mnt/usersData/Viral_CHO/analysis/{sample}/trimmed/"
                os.makedirs(directory_out, exist_ok=True)
                chunk_fq_files(directory_in, directory_out, filename, True)

        # MULTI PROCESSING
        if multi_processing:
            with Pool(processes=5) as p:
                p.map(
                    multiprocess_sample, samples
                )  # process data_inputs iterable with pool
