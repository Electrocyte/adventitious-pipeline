#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 13:26:34 2021

@author: mangi
"""


def check_sample_name(samples: list) -> None:
    for sample in samples:
        print(f"Sample: {sample}")
        sample_no = sample.split("_")
        if len(sample_no) != 6:
            print(f"Sample contains TOO MANY labels ({len(sample_no)}): {sample}")

        datapoints = {"date": 0, "batch": 4, "duration": 5}
        for k, v in datapoints.items():
            try:
                item = sample_no[v]
                if v == 5:
                    item = item.replace("h", "")
                check = int(item)

            except:
                print(f"ERROR! {k} is incorrect: {sample_no[v]}, \t\t{sample}")

        if not "CFU" in sample_no[3] and not "PFU" in sample_no[3]:
            print(
                f"ERROR! Concentration format incorrect! {sample_no[3]}, \t\t{sample}"
            )


if __name__ == "__main__":
    samples = [
        "20211f102_RNA_PA-BE-BA-BL_1000CFU_1_10h",
        "20211102_RNA_PA-CE-BA-BL_1000PFU_1_10h",
        "20211102_RNA_PA-CE-CA-CL_1000CFU_1_10h",
        "20211105_RNA_PA-CE-BA-BL_800CFU_1_10h",
        "20211105_RNA_PA-CE-BA-BL-PL_800CFU_1_10h",
        "20211105_RNA_PA-CE-BA-BL_800CFU_PL_1_10h",
        "20211105_RNA_PA-CE-CA-CBL_80_1_10h",
        "20211105_RNA_PA-CE-CA-CBL_8000000000000000000000000_1_10h",
        "20211105_RNA_PA-CE-CA-CL_800CFU_1_10",
        "20211105_RNA_PA-CE-CA-CL_800CFU_1_foo",
    ]
    check_sample_name(samples)
