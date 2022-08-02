#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 11:36:49 2022

@author: mangi
"""

import glob
import os
import shutil
import gzip


def move_file(files: list, out: str) -> list:
    print(f"Moving fastq files to: {out}")
    for file in files:
        filename = file.replace("\\", "/").split("/")[-1:][0]
        print(f"cp", file, f"{out}/{filename}")
        shutil.copy(file, f"{out}/{filename}")


def de_gunzip_file(files: list) -> None:
    for file in files:
        with gzip.open(file, 'rb') as f_in:
            with open(file.replace(".gz", ""), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(file, file.replace(".gz", ""))
            

def main(folder_builder: list, directory: str) -> None:
    for (sample, experiment, barcode) in folder_builder:

        # save sequencing summary here
        os.makedirs(f"{directory}/{sample}/", exist_ok=True)
        os.makedirs(f"{directory}/{sample}/{barcode}/", exist_ok=True)

        sequencing_summary = glob.glob(
            f"{directory}/{experiment}/**/sequencing_summary*.txt"
        )

        if len(sequencing_summary) > 0:
            with open(sequencing_summary[0], "r") as f, open(
                f"{directory}/{sample}/{barcode}/sequencing_summary.txt", "w"
            ) as ff:
                for n, line in enumerate(f):
                    if n == 0:
                        ff.write(line)
                    if barcode in line:
                        ff.write(line)

        # save fastq files here from demultiplexed guppy folders
        # fq_files = glob.glob(f'{directory}/{experiment}/**/fastq_pass/{barcode}/*q')[0]
        # folder_name = "/".join(fq_files.replace("\\", "/").split("/")[-4:-3])

        pass_sample = f"{directory}/{sample}/{barcode}/fastq_pass/"
        fail_sample = f"{directory}/{sample}/{barcode}/fastq_fail/"

        os.makedirs(f"{directory}/{sample}/", exist_ok=True)
        os.makedirs(f"{directory}/{sample}/{barcode}/", exist_ok=True)
        os.makedirs(pass_sample, exist_ok=True)
        os.makedirs(fail_sample, exist_ok=True)
        
        # check for gunzip
        gz_fq_pass_files = glob.glob(
            f"{directory}/{experiment}/**/fastq_pass/{barcode}/*q*z"
        )
        gz_fq_fail_files = glob.glob(
            f"{directory}/{experiment}/**/fastq_fail/{barcode}/*q*z"
        )
        if len(gz_fq_fail_files) > 0:
            de_gunzip_file(gz_fq_fail_files)
        if len(gz_fq_pass_files) > 0:
            de_gunzip_file(gz_fq_pass_files)
        
        fq_pass_files = glob.glob(
            f"{directory}/{experiment}/**/fastq_pass/{barcode}/*q"
        )
        fq_fail_files = glob.glob(
            f"{directory}/{experiment}/**/fastq_fail/{barcode}/*q"
        )

        if len(fq_pass_files) > 0:
            move_file(fq_pass_files, pass_sample)
        if len(fq_fail_files) > 0:
            move_file(fq_fail_files, fail_sample)


if __name__ == "__main__":
    directory = "/mnt/usersData/CRAAM/"
    folder_builder = [
        [
            "20220321_DNA_working-bio-filter-chips-1-rapid-kit_10CFU_8_15",
            "S8B0",
            "barcode01",
        ],
        [
            "20220321_DNA_working-bio-filter-chips-2-rapid-kit_10CFU_8_15",
            "S8B0",
            "barcode02",
        ],
        [
            "20220321_DNA_unused-bio-filter-chips-2-rapid-kit_10CFU_8_15",
            "S8B0",
            "barcode03",
        ],
        ["20220324_DNA_sump-rapid-kit_10CFU_9_5", "S9B0", "barcode04"],
        ["20220324_DNA_drum-tank-1-rapid-kit_10CFU_9_5", "S9B0", "barcode05"],
        ["20220324_DNA_drum-tank-2-rapid-kit_10CFU_9_5", "S9B0", "barcode06"],
    ]
    main(folder_builder, directory)
