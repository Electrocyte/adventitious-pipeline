#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:26:48 2022

@author: mangi
"""

#################

from pathlib import Path
import pandas as pd
from pipeline import preprocess_unseen_data, ML_assessment  
    
unseen_dir = "D:/Dropbox/AA SMART/fourier_data/"
new_samples = 20220701

debug = True
# model location
directory = "D:/SequencingData/Harmonisation/DNA/analysis/"
_ML_out_ = f"{directory}/ML_training-VFB/"
XGB_out = f"{_ML_out_}/OCS-XGBoost-VFB/"
vees = 'PA-Cacnes-Pacnes-EC-Klebpneu-Calbicans-Saureus-Bsubtilis'

database_dict = {"v_f_b":"v_f_b"}
BLASTn_name = list(database_dict.keys())[0]
kingdom = list(database_dict.values())[0]

# # SUBSET SAMPLES
subset = True
independent_var = "_aDNA_" # amplicon

###############################################################################

Nanoplot = f"{unseen_dir}/nanoplot_summary_data.csv"
BLAST =      f"{unseen_dir}/describe_{BLASTn_name}_no_host_all_agg.csv"
Centrifuge = f"{unseen_dir}/centrifuge_{kingdom}_describe_meta.csv"

if Path(BLAST).is_file():
    full_nano_cols, Centrifuge_df, Nanoplot_df, BLAST_df = preprocess_unseen_data.load_in_data(BLAST, Centrifuge, Nanoplot, kingdom, BLASTn_name)
    
    Centrifuge_df["strain"] = Centrifuge_df["sample"].str.split(r"_", expand=True)[2]
    BLAST_df['sample'] = BLAST_df['date'].astype(str)+"_"+BLAST_df['NA']+"_"+BLAST_df['strain']+"_"+BLAST_df['concentration_CFU'].astype(str)+"_"+BLAST_df['batch'].astype(str)+"_"+BLAST_df['duration_h'].astype(str)
    
    BLAST_df_samples = [x for x in BLAST_df["sample"].unique() if int(x.split("_")[0]) > new_samples]
    BLAST_df = BLAST_df.loc[BLAST_df['sample'].isin(BLAST_df_samples)]
    Centrifuge_df_samples = [x for x in Centrifuge_df["sample"].unique() if int(x.split("_")[0]) > new_samples]
    Centrifuge_df = Centrifuge_df.loc[Centrifuge_df['sample'].isin(Centrifuge_df_samples)]

    BLAST_df_read_clean = BLAST_df.groupby(["sample", "genus_species"])["length_count"].mean()
    BLAST_df_read_clean_s = pd.DataFrame(BLAST_df_read_clean)
    BLAST_df_read_clean_s.reset_index(inplace=True)
    read_count_correction = BLAST_df_read_clean_s.groupby(["sample"]).sum()
    # read_count_correction.to_csv(f"{_ML_out_}/blast-classified-read-count-mean.csv")
    
    BLAST_df = BLAST_df[list(BLAST_df.columns[-2:]) + list(BLAST_df.columns[:-2])]   

    # this step may be required if training the dataset takes too long
    Centrifuge_df = Centrifuge_df.loc[Centrifuge_df["score_mean"] > 900]
    Centrifuge_df.reset_index(drop=True,inplace=True)
    
    # check if this is the lowest value in other data sets
    BLAST_df = BLAST_df.loc[BLAST_df["pident_max"] > 83]
    BLAST_df.reset_index(drop=True,inplace=True)

    Centrifuge_df, BLAST_df = preprocess_unseen_data.count_for_cols(Centrifuge_df, BLAST_df, "std", "std_nans")
    Centrifuge_df = preprocess_unseen_data.value_added(Centrifuge_df, ["name", "sample"], "score_count")
    BLAST_df = preprocess_unseen_data.value_added(BLAST_df, ["name", "sample"], "length_count")
    Centrifuge_df = preprocess_unseen_data.clean_centrifuge(Centrifuge_df)
    Centrifuge_df["abundance"] = Centrifuge_df["numUniqueReads"] / Centrifuge_df["totalUniqReads"]
    
    if subset:
        Centrifuge_df, BLAST_df, Nanoplot_df = preprocess_unseen_data.subset_samples(Centrifuge_df, BLAST_df, Nanoplot_df, "_0_", independent_var)
       
    cn_df = pd.merge(Centrifuge_df, Nanoplot_df, how='left', on=['sample'])  
    bn_df = pd.merge(BLAST_df, Nanoplot_df, how='left', on=['sample'])  
    
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
        c_idx = preprocess_unseen_data.check_NaNs(cn_df, c_col)
        cn_df = cn_df.iloc[c_idx]

    for b_col in hd_data_cols:
        b_idx = preprocess_unseen_data.check_NaNs(bn_df, b_col)
        bn_df = bn_df.iloc[b_idx]
     
ce_id_cols = ["name", "sample", "seqID"]
bl_id_cols = ["name", "sample", "sseqid"]

################## MACHINE LEARNING ##################
################## TRAINING + TESTING ##################

pd.options.mode.chained_assignment = None  # default='warn'

unseen_df_ce = cn_df.drop_duplicates(subset=ce_id_cols)
unseen_df_ce.reset_index(inplace=True, drop=True)
unseen_df_bn = bn_df.drop_duplicates(subset=bl_id_cols)
unseen_df_bn.reset_index(inplace=True, drop=True)

samp_pivot_ce, sample_status_ce = preprocess_unseen_data.run_unseen_sample_analysis("centrifuge", vees, XGB_out, ce_id_cols, debug, unseen_df_ce)
samp_pivot_hs, sample_status_hs = preprocess_unseen_data.run_unseen_sample_analysis("BLAST", vees, XGB_out, bl_id_cols, debug, unseen_df_bn)

genuine_pivot_ce, genuine_status_ce = preprocess_unseen_data.run_unseen_genuine_analysis("centrifuge", vees, XGB_out, ce_id_cols, debug, unseen_df_ce, sample_status_ce)
genuine_pivot_hs, genuine_status_hs = preprocess_unseen_data.run_unseen_genuine_analysis("BLAST", vees, XGB_out, bl_id_cols, debug, unseen_df_bn, sample_status_hs)

# ML_assessment.calculate_twin_model_accuracy() requires implementation for the model to correctly predict the sterility status.
#
#

if debug:
    species = {"Pseudomonas aeruginosa":["PA", "9027"], "Cutibacterium acnes":["Cacnes","Pacnes"], \
                "Escherichia coli":["EC"], "Klebsiella pneumoniae":["Klebpneu"], \
                  "Candida":["Calbicans"], "Staphylococcus aureus":["Saureus"], \
                    "Bacillus subtilis": ["Bsubtilis"]}    
        
    sample_status_ce = preprocess_unseen_data.debug_true_labels(XGB_out, sample_status_ce)
    ML_assessment.run_class_accuracy(sample_status_ce, "centrifuge", "Correct", "unseen-sample", "sample_encoded_mask", "XGB_sample_prediction")
    
    sample_status_hs = preprocess_unseen_data.debug_true_labels(XGB_out, sample_status_hs)
    ML_assessment.run_class_accuracy(sample_status_hs, "BLAST", "Correct", "unseen-sample", "sample_encoded_mask", "XGB_sample_prediction")
    
    genuine_status_ce = preprocess_unseen_data.apply_mask(genuine_status_ce, species)
    ML_assessment.run_class_accuracy(genuine_status_ce, "centrifuge", "Correct", "unseen-genuine", "mask", "XGB_genuine_prediction")
    
    genuine_status_hs = preprocess_unseen_data.apply_mask(genuine_status_hs, species)
    ML_assessment.run_class_accuracy(genuine_status_hs, "BLAST", "Correct", "unseen-genuine", "mask", "XGB_genuine_prediction")














