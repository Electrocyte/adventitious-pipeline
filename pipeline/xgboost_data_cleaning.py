#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 11:28:04 2022

@author: mangi
"""

################## DATA CLEANING ##################

from pathlib import Path
import pandas as pd
import numpy as np
from functools import partial
import seaborn as sns
import matplotlib.pyplot as plt
from string import whitespace
import matplotlib.patches as mpatches


def load_in_data(
    BLAST: pd.DataFrame,
    Centrifuge: pd.DataFrame,
    Nanoplot: pd.DataFrame,
    kingdom: str,
    BLASTn_name: str,
) -> (list, pd.DataFrame, pd.DataFrame, pd.DataFrame):
    print(f"BLAST file: {BLAST}")
    print(f"Centrifuge file: {Centrifuge}")
    BLAST_df = pd.read_csv(BLAST)
    BLAST_df["db"] = BLASTn_name

    Centrifuge_df = pd.read_csv(Centrifuge)
    Centrifuge_df["db"] = kingdom

    Nanoplot_df = pd.read_csv(Nanoplot)
    Nanoplot_df = Nanoplot_df.rename(columns={"index": "sample"})
    Nanoplot_df["sample"] = Nanoplot_df["sample"].replace(
        to_replace=r"CFU", value="", regex=True
    )
    nano_cols = [
        "Meanreadlength",
        "Medianreadlength",
        "Numberofreads",
        "ReadlengthN50",
        "Totalbases",
    ]
    for n_col in nano_cols:
        Nanoplot_df[n_col] = Nanoplot_df[n_col].str.replace(",", "").astype(float)
    full_nano_cols = [
        "Activechannels",
        "Meanreadlength",
        "Meanreadquality",
        "Medianreadlength",
        "Medianreadquality",
        "Numberofreads",
        "ReadlengthN50",
        "Totalbases",
    ]

    return full_nano_cols, Centrifuge_df, Nanoplot_df, BLAST_df


def check_NaNs(df: pd.DataFrame, col: str) -> list:
    idx = set(np.where(df[col].notnull())[0])
    return list(idx)
    

def check_TPs(df: pd.DataFrame, species: dict):
    df.reset_index(inplace=True,drop=True)
    true_masks = []
    for k, vv in species.items():
        if isinstance(vv, list):
            for v in vv:
                true_mask = df.loc[
                    (df["strain"].str.contains(v)) & (df.name.str.contains(k))
                ].index
                df.loc[true_mask, "mask"] = "True_positive"
                true_masks.append(true_mask)
        else:
            true_mask = df.loc[
                (df["strain"].str.contains(vv)) & (df.name.str.contains(k))
            ].index
            df.loc[true_mask, "mask"] = "True_positive"
            true_masks.append(true_mask)
    true_masks = list(set([item for sublist in true_masks for item in sublist]))
    TPs = df.iloc[true_masks]

    return list(TPs["sample"].unique()), TPs
    

def std_val_func(col_type: str, row: str) -> int:
    counter = 1  
    cols = [i for i in row.index if col_type in i]
    for col in cols:
        if not pd.isnull(row[col]):
            counter += 1
    return counter  


def count_for_cols(ce_df: pd.DataFrame, hs_df: pd.DataFrame, col_type: str, col_out: str) -> \
                    (pd.DataFrame, pd.DataFrame):
    
    partial_func = partial(std_val_func, col_type)

    ce_df[col_out] = ce_df.apply(partial_func, axis = 1)
    hs_df[col_out] = hs_df.apply(partial_func, axis = 1)
    
    return ce_df, hs_df


def value_added(df: pd.DataFrame, cols: list, classifier_col: str) -> pd.DataFrame:
    vc = df[cols].value_counts()
    col = "-".join(cols)
    vc_i = vc.to_frame().reset_index().rename(columns={0: f"{col}-count"})
    df_merge = pd.merge(df, vc_i, on=cols, how='outer')
    df_merge[f"vc-{col}-fraction"] = df_merge[classifier_col] / (df_merge[classifier_col] + df_merge[f"{col}-count"])
    df_merge["read_qc"] = df_merge["mean_qscore_template_count"] / (df_merge["mean_qscore_template_count"] + df_merge["mean_qscore_template_max"])
    return df_merge
        
        
def clean_centrifuge(df: pd.DataFrame) -> pd.DataFrame:
    dfs = []
    for name, group in df.groupby(["sample","name"]):
        if len(group) < 5:
            dfs.append(group)
        else:
            group = group.nlargest(5, "score_count")
            group.reset_index(inplace=True,drop=True)
            dfs.append(group)
    cat_dfs = pd.concat(dfs)
    return cat_dfs


def subset_samples(df1: pd.DataFrame, df2: pd.DataFrame, df3: pd.DataFrame, \
                   negs: str, indep: str) -> \
    (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    
    # spike = "TC"
    # pure = "pure"
    # direct = "_DNA_"
    # amplicon = "_aDNA_"
    
    dfs = []
    
    if indep == "pure":
        for df in [df1, df2, df3]:
            NCs = [i for i in list(df1["sample"].unique()) if negs in i]
            TCs = [i for i in list(df1["sample"].unique()) if "TC" not in i]
            keep = list(set(NCs) | set(TCs))
            new_df = df.loc[df["sample"].isin(keep)]
            dfs.append(new_df)
        
    else:
        for df in [df1, df2, df3]:
            NCs = [i for i in list(df1["sample"].unique()) if negs in i]
            TCs = [i for i in list(df1["sample"].unique()) if indep in i]  
            keep = list(set(NCs) | set(TCs))
            new_df = df.loc[df["sample"].isin(keep)]
            dfs.append(new_df)
    
    df1_out, df2_out, df3_out = dfs
    
    return df1_out, df2_out, df3_out


def assess_quality(directory: str, BLASTn_name: str, 
                              kingdom: str, independent_var: str,
                              species: dict, plot: bool) -> None:    

    def clean_strings(string: str, cut_off: int) -> str:
        # from string import whitespace
        # string = "Cumulative frequency for detection of MVM reads per minute for batch19"
        len_str = len(string)
        mid_point = int(len_str / 2)
        whitespaces = [i for i, char in enumerate(string) if char in whitespace]
        cut_string = min(whitespaces, key=lambda x: abs(x - mid_point))
        new_string = f"{string[0:cut_string]}\n{string[cut_string:]}"
        
        subsplits = ""
        if len_str > cut_off:
            split_again = new_string.split("\n")
            for new_str in split_again:
                subsplit = clean_strings(new_str, cut_off)
                subsplits += subsplit
    
        if len(subsplits) > 0:
            len_str = len(subsplits)
            mid_point = int(len_str / 2)
            whitespaces = [i for i, char in enumerate(subsplits) if char in whitespace]
            cut_string = min(whitespaces, key=lambda x: abs(x - mid_point))
            new_string = f"{subsplits[0:cut_string]}\n{subsplits[cut_string:]}"
        
        return new_string
    
    Nanoplot = f"{directory}/nanoplot_summary_data.csv"
    BLAST =      f"{directory}/describe_{BLASTn_name}_no_host_all_agg.csv"
    Centrifuge = f"{directory}/centrifuge_{kingdom}_describe_meta.csv"
    
    full_nano_cols, Centrifuge_df, Nanoplot_df, BLAST_df = load_in_data(BLAST, Centrifuge, Nanoplot, kingdom, BLASTn_name)
    Centrifuge_df["strain"] = Centrifuge_df["sample"].str.split(r"_", expand=True)[2]
    BLAST_df['sample'] = BLAST_df['date'].astype(str)+"_"+BLAST_df['NA']+"_"+BLAST_df['strain']+"_"+BLAST_df['concentration_CFU'].astype(str)+"_"+BLAST_df['batch'].astype(str)+"_"+BLAST_df['duration_h'].astype(str)
   
    a_before_ce, ce_bef_df = check_TPs(Centrifuge_df, species)    
    a_before_hs, hs_bef_df = check_TPs(BLAST_df, species)

    BLAST_df_read_clean = BLAST_df.groupby(["sample", "genus_species"])["length_count"].mean()
    BLAST_df_read_clean_s = pd.DataFrame(BLAST_df_read_clean)
    BLAST_df_read_clean_s.reset_index(inplace=True)
    read_count_correction = BLAST_df_read_clean_s.groupby(["sample"]).sum()
    read_count_correction.to_csv(f"{_ML_out_}/blast-classified-read-count-mean.csv")
        
    BLAST_df = BLAST_df[list(BLAST_df.columns[-2:]) + list(BLAST_df.columns[:-2])]   
    
    #################### BLAST ####################

    hs_fps = BLAST_df.loc[~BLAST_df.index.isin(hs_bef_df.index)]   
    
    # only examine amplified DNA
    hs_bef_df = hs_bef_df.loc[hs_bef_df["sample"].str.contains(independent_var)]
    hs_fps = hs_fps.loc[hs_fps["sample"].str.contains(independent_var)]   
    
    h_fpshort = hs_fps.loc[hs_fps["pident_max"] < 83]
    h_fplong = hs_fps.loc[hs_fps["pident_max"] >= 83]
    h_tpshort = hs_bef_df.loc[hs_bef_df["pident_max"] < 83]
    h_tplong = hs_bef_df.loc[hs_bef_df["pident_max"] >= 83]

    h_fp_means = pd.DataFrame(hs_fps.groupby(["sample"])["pident_max"].mean())
    h_tp_means = pd.DataFrame(hs_bef_df.groupby(["sample"])["pident_max"].mean())
    
    #################### CENTRIFUGE ####################

    Centrifuge_df["abundance"] = Centrifuge_df["numUniqueReads"] / Centrifuge_df["totalUniqReads"]
    ce_fps = Centrifuge_df.loc[~Centrifuge_df.index.isin(ce_bef_df.index)]   
    # ce_fps = ce_fps.loc[ce_fps["abundance"] >= 1e-5]
    
    # only examine amplified DNA
    ce_bef_df["abundance"] = ce_bef_df["numUniqueReads"] / ce_bef_df["totalUniqReads"]
    ce_bef_df = ce_bef_df.loc[ce_bef_df["sample"].str.contains(independent_var)]
    # ce_bef_df = ce_bef_df.loc[ce_bef_df["abundance"] >= 1e-5]
    
    ce_fps = ce_fps.loc[ce_fps["sample"].str.contains(independent_var)]    
    
    c_fpshort = ce_fps.loc[(ce_fps["score_mean"] < 900)]
    # c_fpshort = c_fpshort.loc[(c_fpshort["abundance"] >= 1e-5)]
    
    c_fplong = ce_fps.loc[(ce_fps["score_mean"] >= 900)]
    # c_fplong = c_fplong.loc[(c_fplong["abundance"] >= 1e-5)]
    
    c_tpshort = ce_bef_df.loc[(ce_bef_df["score_mean"] < 900)]
    # c_tpshort = c_tpshort.loc[(c_tpshort["abundance"] >= 1e-5)]
    
    c_tplong = ce_bef_df.loc[(ce_bef_df["score_mean"] >= 900)]
    # c_tplong = c_tplong.loc[(c_tplong["abundance"] >= 1e-5)]    
    
    c_fp_means = pd.DataFrame(ce_fps.groupby(["sample"])["score_mean"].mean())
    c_tp_means = pd.DataFrame(ce_bef_df.groupby(["sample"])["score_mean"].mean())
    
    ####################
    
    years = ['2020','2021','2022']
    
    blue_patch = mpatches.Patch(color="blue", label='TP')
    black_patch = mpatches.Patch(color="black", label='FP')
    
    if plot:
        for part in years:
            h_year_fp = h_fp_means.loc[h_fp_means.index.str.contains(part)]
            h_year_tp = h_tp_means.loc[h_tp_means.index.str.contains(part)]
            
            f, ax = plt.subplots(figsize=(30, 15))
            
            sns.scatterplot(data=h_year_fp, x=h_year_fp.index, y="pident_max", color="black", s = 150)
            sns.scatterplot(data=h_year_tp, x=h_year_tp.index, y="pident_max", color="blue", s = 150)  
            
            plt.legend(
                loc=0, prop={"size": 25}, markerscale=4, handles=[blue_patch, black_patch])
                # ) bbox_to_anchor=(1.05, 0.75)
        
            ax.set_ylabel("Number of reads", size=40)
            ax.set_ylim([10, 100])
            
            title = f"Mean reads for sample plotted for true positive and false positive species ({part}) - BLAST"
            newline_title = clean_strings(title, 120)
        
            plt.title(newline_title, size=50)
            plt.tight_layout()
            plt.xticks(rotation=90)
            ax.tick_params(axis="x", labelsize=20)
            ax.tick_params(axis="y", labelsize=30)
            plt.show()
            ################################
    
            h_sh_year_tp = h_tpshort.loc[h_tpshort["sample"].str.contains(part)]
            h_lo_year_tp = h_tplong.loc[h_tplong["sample"].str.contains(part)]
            
            f, ax = plt.subplots(figsize=(30, 15))
        
            sns.scatterplot(data=h_sh_year_tp, x="sample", y="pident_max", color="black", s=150)
            sns.scatterplot(data=h_lo_year_tp, x="sample", y="pident_max", color="blue", s=150)
    
            ax.set_ylabel("Number of reads", size=40)
            ax.set_ylim([10, 100])
            
            title = f"All predicted species reads for sample plotted for true positive species ({part}) - BLAST"
            newline_title = clean_strings(title, 120)
        
            plt.title(newline_title, size=50)
            plt.tight_layout()
            plt.xticks(rotation=90)
            ax.tick_params(axis="x", labelsize=20)
            ax.tick_params(axis="y", labelsize=30)
            plt.show()
            ################################
    
            h_sh_year_fp = h_fpshort.loc[h_fpshort["sample"].str.contains(part)]
            h_lo_year_fp = h_fplong.loc[h_fplong["sample"].str.contains(part)]
            
            f, ax = plt.subplots(figsize=(30, 15))
        
            sns.scatterplot(data=h_sh_year_fp, x="sample", y="pident_max", color="black", s=150)
            sns.scatterplot(data=h_lo_year_fp, x="sample", y="pident_max", color="red", s=150)
    
            ax.set_ylabel("Number of reads", size=40)
            ax.set_ylim([10, 100])
            
            title = f"All predicted species reads for sample plotted for false positive species ({part}) - BLAST"
            newline_title = clean_strings(title, 120)
        
            plt.title(newline_title, size=50)
            plt.tight_layout()
            plt.xticks(rotation=90)
            ax.tick_params(axis="x", labelsize=20)
            ax.tick_params(axis="y", labelsize=30)
            plt.show()
    
        
        for part in years:
            c_year_fp = c_fp_means.loc[c_fp_means.index.str.contains(part)]
            c_year_tp = c_tp_means.loc[c_tp_means.index.str.contains(part)]
            
            f, ax = plt.subplots(figsize=(30, 15))
            
            sns.scatterplot(data=c_year_fp, x=c_year_fp.index, y="score_mean", color="black", s = 150)
            sns.scatterplot(data=c_year_tp, x=c_year_tp.index, y="score_mean", color="blue", s = 150)  
            
            plt.legend(
                loc=0, prop={"size": 25}, markerscale=4, handles=[blue_patch, black_patch])
                # ) bbox_to_anchor=(1.05, 0.75)
        
            ax.set_ylabel("Number of reads", size=40)
            ax.set_yscale("log")
            ax.set_ylim([10, 1e5])
            
            title = f"Mean reads for sample plotted for true positive and false positive species ({part}) - CENTRIFUGE"
            newline_title = clean_strings(title, 120)
        
            plt.title(newline_title, size=50)
            plt.tight_layout()
            plt.xticks(rotation=90)
            ax.tick_params(axis="x", labelsize=20)
            ax.tick_params(axis="y", labelsize=30)
            plt.show()
            ################################
    
            c_sh_year_tp = c_tpshort.loc[c_tpshort["sample"].str.contains(part)]
            c_lo_year_tp = c_tplong.loc[c_tplong["sample"].str.contains(part)]
            
            f, ax = plt.subplots(figsize=(30, 15))
        
            sns.scatterplot(data=c_sh_year_tp, x="sample", y="score_mean", color="black", s=150)
            sns.scatterplot(data=c_lo_year_tp, x="sample", y="score_mean", color="blue", s=150)
    
            ax.set_ylabel("Number of reads", size=40)
            ax.set_yscale("log")
            ax.set_ylim([10, 1e5])
            
            title = f"All predicted species reads for sample plotted for true positive species ({part}) - CENTRIFUGE"
            newline_title = clean_strings(title, 120)
        
            plt.title(newline_title, size=50)
            plt.tight_layout()
            plt.xticks(rotation=90)
            ax.tick_params(axis="x", labelsize=20)
            ax.tick_params(axis="y", labelsize=30)
            plt.show()
            ################################
    
            c_sh_year_fp = c_fpshort.loc[c_fpshort["sample"].str.contains(part)]
            c_lo_year_fp = c_fplong.loc[c_fplong["sample"].str.contains(part)]
            
            f, ax = plt.subplots(figsize=(30, 15))
        
            sns.scatterplot(data=c_sh_year_fp, x="sample", y="score_mean", color="black", s=150)
            sns.scatterplot(data=c_lo_year_fp, x="sample", y="score_mean", color="red", s=150)
    
            ax.set_ylabel("Number of reads", size=40)
            ax.set_yscale("log")
            ax.set_ylim([10, 1e5])
            
            title = f"All predicted species reads for sample plotted for false positive species ({part}) - CENTRIFUGE"
            newline_title = clean_strings(title, 120)
        
            plt.title(newline_title, size=50)
            plt.tight_layout()
            plt.xticks(rotation=90)
            ax.tick_params(axis="x", labelsize=20)
            ax.tick_params(axis="y", labelsize=30)
            plt.show()
        
    #################### CENTRIFUGE ####################
    c_tp_per = (len(c_tpshort) / (len(ce_bef_df)))*100
    c_fp_per = (len(c_fpshort) / (len(ce_fps)))*100
    c_lost_sample_tp_name = set(c_tp_means.index) - set(c_tplong["sample"].unique())
    
    z_tp_sample_lost = len(c_tpshort["sample"].unique())
    z_actual_tps_lost = len(set(c_tpshort["sample"].unique()) - set(c_lost_sample_tp_name))
    z_tp_sample_kept = len(c_tplong["sample"].unique())
    z_total_tp = len(c_tp_means)
    
    samples_lost = z_tp_sample_kept / z_total_tp * 100
    no_sa_lost = z_tp_sample_lost - z_actual_tps_lost

    print(f"CENTRIFUGE - Total true SAMPLES retained: {samples_lost:.2f}%; or {no_sa_lost}/{z_total_tp} samples lost")    
    print(f"TPs (total predicted spp) removed: {c_tp_per:.2f}%; \nFPs (total predicted spp) removed: {c_fp_per:.2f}")
    #################### CENTRIFUGE ####################
    
    h_tp_per = (len(h_tpshort) / (len(hs_bef_df)))*100
    h_fp_per = (len(h_fpshort) / (len(hs_fps)))*100
    h_lost_sample_tp_name = set(c_tp_means.index) - set(c_tplong["sample"].unique())
    
    z_tp_sample_lost = len(h_tpshort["sample"].unique())
    z_actual_tps_lost = len(set(h_tpshort["sample"].unique()) - set(h_lost_sample_tp_name))
    z_tp_sample_kept = len(h_tplong["sample"].unique())
    z_total_tp = len(h_tp_means)
    
    samples_lost = z_tp_sample_kept / z_total_tp * 100
    no_sa_lost = z_tp_sample_lost - z_actual_tps_lost

    print(f"BLAST - Total true SAMPLES retained: {samples_lost:.2f}%; or {no_sa_lost}/{z_total_tp} samples lost")    
    print(f"TPs (total predicted spp) removed: {h_tp_per:.2f}%; \nFPs (total predicted spp) removed: {h_fp_per:.2f}")
    

def main(directory: str,
        BLASTn_name: str,
        kingdom: str,
        species: dict, 
        _ML_out_: str,
        subset: bool,
        independent_var: str) -> (list, 
                                  list, 
                                  pd.DataFrame, 
                                  pd.DataFrame):    
    
    Nanoplot = f"{directory}/nanoplot_summary_data.csv"
    
    BLAST =      f"{directory}/describe_{BLASTn_name}_no_host_all_agg.csv"
    Centrifuge = f"{directory}/centrifuge_{kingdom}_describe_meta.csv"
    if Path(BLAST).is_file():
        full_nano_cols, Centrifuge_df, Nanoplot_df, BLAST_df = load_in_data(BLAST, Centrifuge, Nanoplot, kingdom, BLASTn_name)
        Centrifuge_df["strain"] = Centrifuge_df["sample"].str.split(r"_", expand=True)[2]
        BLAST_df['sample'] = BLAST_df['date'].astype(str)+"_"+BLAST_df['NA']+"_"+BLAST_df['strain']+"_"+BLAST_df['concentration_CFU'].astype(str)+"_"+BLAST_df['batch'].astype(str)+"_"+BLAST_df['duration_h'].astype(str)
        
        # a_before_ce, ce_bef_df = check_TPs(Centrifuge_df, species)
        # a_before_hs, hs_bef_df = check_TPs(BLAST_df, species)
        
        BLAST_df_read_clean = BLAST_df.groupby(["sample", "genus_species"])["length_count"].mean()
        BLAST_df_read_clean_s = pd.DataFrame(BLAST_df_read_clean)
        BLAST_df_read_clean_s.reset_index(inplace=True)
        read_count_correction = BLAST_df_read_clean_s.groupby(["sample"]).sum()
        read_count_correction.to_csv(f"{_ML_out_}/blast-classified-read-count-mean.csv")
        
        BLAST_df = BLAST_df[list(BLAST_df.columns[-2:]) + list(BLAST_df.columns[:-2])]   
    
        # this step may be required if training the dataset takes too long
        Centrifuge_df = Centrifuge_df.loc[Centrifuge_df["score_mean"] > 900]
        Centrifuge_df.reset_index(drop=True,inplace=True)
        
        # check if this is the lowest value in other data sets
        BLAST_df = BLAST_df.loc[BLAST_df["pident_max"] > 83]
        BLAST_df.reset_index(drop=True,inplace=True)
    
        Centrifuge_df, BLAST_df = count_for_cols(Centrifuge_df, BLAST_df, "std", "std_nans")
        Centrifuge_df = value_added(Centrifuge_df, ["name", "sample"], "score_count")
        BLAST_df = value_added(BLAST_df, ["name", "sample"], "length_count")
        Centrifuge_df = clean_centrifuge(Centrifuge_df)
        Centrifuge_df["abundance"] = Centrifuge_df["numUniqueReads"] / Centrifuge_df["totalUniqReads"]
        # Centrifuge_df = Centrifuge_df.loc[Centrifuge_df["abundance"] >= 1e-5]
        # Centrifuge_df.reset_index(drop=True,inplace=True)  
        
        if subset:
            Centrifuge_df, BLAST_df, Nanoplot_df = subset_samples(Centrifuge_df, BLAST_df, Nanoplot_df, "_0_", independent_var)
           
        cn_df = pd.merge(Centrifuge_df, Nanoplot_df, how='left', on=['sample'])  
        bn_df = pd.merge(BLAST_df, Nanoplot_df, how='left', on=['sample'])  
        
        # mostly removes Direct DNA samples.
        a_clean_ce, ce_clean_df = check_TPs(cn_df, species)
        a_clean_hs, hs_clean_df = check_TPs(bn_df, species)
        
        cd_data_cols = [ 'numReads', 'numUniqueReads', 'abundance',
           'totalUniqReads', 'score_count', 'score_mean',
           'score_std', 'score_min', 'score_25%', 'score_50%', 'score_75%',
           'score_max', '2ndBestScore_count', '2ndBestScore_mean',
           '2ndBestScore_std', '2ndBestScore_min', '2ndBestScore_25%',
           '2ndBestScore_50%', '2ndBestScore_75%', '2ndBestScore_max',
           'hitLength_count', 'hitLength_mean', 'hitLength_std', 'hitLength_min',
           'hitLength_25%', 'hitLength_50%', 'hitLength_75%', 'hitLength_max',
           'queryLength_count', 'queryLength_mean', 'queryLength_std',
           'queryLength_min', 'queryLength_25%', 'queryLength_50%',
           'queryLength_75%', 'queryLength_max', 'numMatches_count',
           'numMatches_mean', 'numMatches_std', 'numMatches_min', 'numMatches_25%',
           'numMatches_50%', 'numMatches_75%', 'numMatches_max',
           'mean_qscore_template_count', 'mean_qscore_template_mean',
           'mean_qscore_template_std', 'mean_qscore_template_min',
           'mean_qscore_template_25%', 'mean_qscore_template_50%',
           'mean_qscore_template_75%', 'mean_qscore_template_max', 
           'std_nans', 'name-sample-count', 'vc-name-sample-fraction', 'read_qc',
           'Activechannels', 'Meanreadlength', 'Meanreadquality',
           'Medianreadlength', 'Medianreadquality', 'Numberofreads',
           'ReadlengthN50', 'Totalbases']
        
        hd_data_cols = ['length_count', 'length_mean', 'length_std',
           'length_min', 'length_25%', 'length_50%', 'length_75%', 'length_max',
           'pident_count', 'pident_mean', 'pident_std', 'pident_min', 'pident_25%',
           'pident_50%', 'pident_75%', 'pident_max', 'bitscore_count',
           'bitscore_mean', 'bitscore_std', 'bitscore_min', 'bitscore_25%',
           'bitscore_50%', 'bitscore_75%', 'bitscore_max', 'mismatches_count',
           'mismatches_mean', 'mismatches_std', 'mismatches_min', 'mismatches_25%',
           'mismatches_50%', 'mismatches_75%', 'mismatches_max', 'gap_opens_count',
           'gap_opens_mean', 'gap_opens_std', 'gap_opens_min', 'gap_opens_25%',
           'gap_opens_50%', 'gap_opens_75%', 'gap_opens_max', 'evalue_count',
           'evalue_mean', 'evalue_std', 'evalue_min', 'evalue_25%', 'evalue_50%',
           'evalue_75%', 'evalue_max', 'mean_qscore_template_count',
           'mean_qscore_template_mean', 'mean_qscore_template_std',
           'mean_qscore_template_min', 'mean_qscore_template_25%',
           'mean_qscore_template_50%', 'mean_qscore_template_75%',
           'mean_qscore_template_max', 
           'b_score', 'std_nans', 'name-sample-count',
           'vc-name-sample-fraction', 'read_qc', 'Activechannels',
           'Meanreadlength', 'Meanreadquality', 'Medianreadlength',
           'Medianreadquality', 'Numberofreads', 'ReadlengthN50', 'Totalbases']
        
        for c_col in cd_data_cols:
            c_idx = check_NaNs(cn_df, c_col)
            cn_df = cn_df.iloc[c_idx]
    
        for b_col in hd_data_cols:
            b_idx = check_NaNs(bn_df, b_col)
            bn_df = bn_df.iloc[b_idx]
            
        # removes low quality amplicon-seq
        # checking reveals nothing useful lost by this
        a_non_nan_ce, ce_nonan_df = check_TPs(cn_df, species)
        a_non_nan_hs, hs_nonan_df = check_TPs(bn_df, species)
        
        # diff_clean_ce = list(set(a_before_ce) - set(a_clean_ce))
        # diff_nonan_ce = list(set(a_clean_ce) - set(a_non_nan_ce))
        # diff_clean_hs = list(set(a_before_hs) - set(a_clean_hs))
        # diff_nonan_hs = list(set(a_clean_hs) - set(a_non_nan_hs))
        
        # diff_clean_df_ce = cn_df.loc[cn_df["sample"].isin(diff_clean_ce)]
        # diff_nonan_df_ce = cn_df.loc[cn_df["sample"].isin(diff_nonan_ce)]
        # diff_clean_df_hs = bn_df.loc[bn_df["sample"].isin(diff_clean_hs)]
        # diff_nonan_df_hs = bn_df.loc[bn_df["sample"].isin(diff_nonan_hs)]    
        
        return hd_data_cols, cd_data_cols, cn_df, bn_df


################## DATA CLEANING ##################

if __name__ == "__main__":
    
    database_dict = {"v_f_b":"v_f_b"}
    BLASTn_name = list(database_dict.keys())[0]
    kingdom = list(database_dict.values())[0]
    species = {"Pseudomonas aeruginosa":["PA"], "Cutibacterium acnes":["Cacnes","Pacnes"], \
            "Escherichia coli":["EC"], "Klebsiella pneumoniae":["Klebpneu"], \
              "Candida albicans":["Calbicans"], "Staphylococcus aureus":["Saureus"], \
                "Bacillus subtilis": ["Bsubtilis"]}
    directory = "D:/SequencingData/Harmonisation/DNA/analysis"
    _ML_out_ = f"{directory}/ML_training-VFB/"
    subset = True
    independent_var = "_aDNA_" # amplicon
    
    hd_data_cols, \
    cd_data_cols, \
        cn_df, \
        bn_df = \
            main(
            directory,
            BLASTn_name,
            kingdom,
            species, 
            _ML_out_,
            subset,
            independent_var)

    assess_quality(directory, BLASTn_name, 
                              kingdom, independent_var,
                              species, False)




# from nonconformist.icp import IcpClassifier, IcpRegressor
# from nonconformist.nc import ClassifierNc, MarginErrFunc, ClassifierAdapter, RegressorNc, AbsErrorErrFunc
# # import copy  

# # https://towardsdatascience.com/how-to-add-uncertainty-estimation-to-your-models-with-conformal-prediction-a5acdb86ea05
# def classifier_calibration_curve(estimator: XGBClassifier, 
#                                  X: np.ndarray, 
#                                  y: np.ndarray, 
#                                  alphas =np.linspace(0,1,10, endpoint=True)):
#     errors = []
#     set_sizes = []
#     for a in alphas:
#         pred = estimator.predict(X, significance=a) # significance = p-value
#         set_sizes.append(np.mean([np.sum(set__) for set__ in pred]))
#         errors.append(1 - np.mean([set_[t] for set_, t in zip(pred, y)]))
#     return errors, set_sizes


# def classification_calibration_plot(dataset_used: str, 
#                                     metagenome_classifier_used: str, 
#                                     estimator: XGBClassifier, 
#                                     X:  np.ndarray, 
#                                     y:  np.ndarray, 
#                                     _type_: str,
#                                     alphas=np.linspace(0,1,10, endpoint=True)):
#     errors, sizes = classifier_calibration_curve(estimator,X,y,alphas)
#     fig, ax1 = plt.subplots(figsize=(15,15))
#     ax2 = ax1.twinx()
#     ax1.plot([0,1], [0,1])
    
#     ax1.plot(alphas, errors, 'o', color = 'black', linewidth=7.0)
#     ax2.plot(alphas, sizes,  '-', color = 'blue', linewidth=7.0)
    
#     ax1.tick_params(axis="x", labelsize=25)
#     ax1.tick_params(axis="y", labelsize=25)
#     ax2.tick_params(axis="y", labelsize=25)
    
#     ax1.set_xlabel('Significance', size=40)
#     ax1.set_ylabel('Error Rate', size=40)
#     ax2.set_ylabel('Avg. Set Size', size=40)
    
#     title = f'Classification Conformal Calibration Curve for {metagenome_classifier_used} on {dataset_used} data for {_type_}'
#     new_title = clean_strings(title, 120)
#     plt.title(new_title, size=40)
#     plt.legend(loc=0, prop={"size": 40}, markerscale=10)
#     plt.show()
    

# def conformal_analysis(dataset_used: str, 
#                        metagenome_classifier_used: str, 
#                        X_test: np.ndarray, 
#                        y_test: np.ndarray, 
#                        estimator: XGBClassifier, 
#                        X_train: np.ndarray, 
#                        y_train: np.ndarray,
#                        _type_: str,
#                        XGB_out: str):
#     import dill as pickle
#     ############ conformal ############
#     filename_save = f"icp_{_type_}-{metagenome_classifier_used}.sav"
#     if len(X_train) > 0:
#         X_calibration, X_test, y_calibration, y_test = train_test_split(X_test, y_test, test_size=0.4, random_state=736)
#         icp = IcpClassifier(ClassifierNc(ClassifierAdapter(estimator), MarginErrFunc()))
#         icp.fit(X_train, y_train)
#         icp.calibrate(X_calibration, y_calibration)
#         icp_model_save = f"{XGB_out}/{filename_save}"
#         print(f"Saving conformal icp model to: {icp_model_save}")
#         pickle.dump(icp, open(icp_model_save, "wb"))
#     else:
#         icp_model_save = f"{XGB_out}/{filename_save}"
#         if Path(icp_model_save).is_file():
#             icp = pickle.load(open(icp_model_save, "rb"))       
            
#     prediction03 = icp.predict(X_test, 0.3)    
#     prediction01 = icp.predict(X_test, 0.1)
#     prediction005 = icp.predict(X_test, 0.05)
#     prediction001 = icp.predict(X_test, 0.01)
#     for sig, prediction in {0.3: prediction03, 0.1: prediction01, 0.05: prediction005, 0.01: prediction001}.items():
#         df = pd.DataFrame(prediction)
#         # df['C'] = df[0].eq(df[1]) # label False / False; True / True = TRUE
#         df['p'] = df[0].ne(df[1]) # label False / True; True / False = TRUE
#         vals = df['p'].value_counts()
#         print(f"\nSignificance: {sig}; FALSE: {vals.loc[False]}; TRUE: {vals.loc[True]}")
#         print(f'Percent TRUE: {(vals.loc[True] / (vals.loc[True] + vals.loc[False])) * 100:.2f} %')
    
#     # from nonconformist.evaluation import ClassIcpCvHelper
#     # from nonconformist.evaluation import class_mean_errors
#     # from nonconformist import evaluation
#     # # evaluation requires updating of library import!!!!!!!!!!!!!!
#     # # from sklearn.model_selection import train_test_split
#     # # from sklearn.model_selection import StratifiedShuffleSplit
#     # # from sklearn.model_selection import KFold  
    
#     # icp_cv = ClassIcpCvHelper(icp)

#     # evaluation.cross_val_score(icp_cv, X_test, y_test, scoring_funcs=[class_mean_errors])
    
#     classification_calibration_plot(dataset_used, metagenome_classifier_used, icp, X_test, y_test, _type_)    
#     ############ conformal ############
    
#     return [prediction01, prediction005]


# metagenome_classifier_used = "centrifuge"
# dataset_used = "train-test"
# input_df = train_test_df_ce
# data_cols = cd_data_cols
# id_cols = ce_id_cols
# _type_ = "sample status"
# estimator = XGB_classifier_model
        