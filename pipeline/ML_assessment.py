#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 13:13:01 2022

@author: mangi
"""

import math    
import pandas as pd
import numpy as np
from functools import partial


def calculate_accuracy(df: pd.DataFrame, col: str, 
                       _type_: str, mg: str) -> float:
    s = df[col].value_counts()
    correct_predictions = s.iloc[0]
    accuracy = (correct_predictions / len(df[col])) 

    return accuracy

# classification error
def calculate_classification_accuracy_interval(acc: float, 
                                               df: pd.DataFrame, 
                                               _type_: str, 
                                               mg: str):
    error = 1 - acc
    n = len(df)
    zs = {"95%": 1.96}
    # zs = {"95%": 1.96, "98%": 2.33, "99%": 2.58}
    
    for percent_ci, z in zs.items():
        e_interval = z * math.sqrt( (error * (1 - error)) / n)
        plus_ci = acc*100 + e_interval*100
        nega_ci = acc*100 - e_interval*100
        print(f"\nConfidence interval is {acc*100:.2f}% +/- {e_interval*100:.2f}% [{nega_ci:.2f} : {plus_ci:.2f}], at {percent_ci} confidence; Data type: {_type_}; Metagenomic classifier: {mg}.")
    print("\n")


def make_decision(row):
    if (
        row["XGB_gen_con_prediction"] == "True_positive"
        and row["XGB_sample_prediction"] == "True_positive"
    ):
        return "True_positive"
    if (
        row["XGB_gen_con_prediction"] == "True_negative"
        and row["XGB_sample_prediction"] == "True_positive"
    ):
        return "False_negative"
    if (
        row["XGB_gen_con_prediction"] == "False_positive"
        and row["XGB_sample_prediction"] == "True_negative"
    ):
        return "True_negative"
    if (
        row["XGB_gen_con_prediction"] == "False_positive"
        and row["XGB_sample_prediction"] == "True_positive"
    ):
        return "False_positive"
    if (
        row["XGB_gen_con_prediction"] == "True_positive"
        and row["XGB_sample_prediction"] == "True_negative"
    ):
        return "False_negative"


def decision_mask(row):
    if row["mask"] == "True_positive" and row["sample_true_mask"] == "True_positive":
        return "True_positive"
    if row["mask"] == "True_negative" and row["sample_true_mask"] == "True_positive":
        return "False_negative"
    if row["mask"] == "Background" and row["sample_true_mask"] == "True_positive":
        return "False_positive"
    if row["mask"] == "Background" and row["sample_true_mask"] == "True_negative":
        return "False_positive"
    if row["mask"] == "False_positive" and row["sample_true_mask"] == "True_negative":
        return "True_negative"
    if row["mask"] == "False_positive" and row["sample_true_mask"] == "True_positive":
        return "False_positive"
    if row["mask"] == "True_positive" and row["sample_true_mask"] == "True_negative":
        return "False_negative"
    

def sample_decision(row):
    # TPs, FPs, TNs, FNs
    if (
        not pd.isna(row["False_positive"])
        and not pd.isna(row["False_negative"])
        and not pd.isna(row["True_negative"])
        and not pd.isna(row["True_positive"])
    ):
        return "Sample contains a lot of noise, could be sterile or contaminated"
    # FPs
    elif (
        pd.isna(row["False_negative"])
        and pd.isna(row["True_negative"])
        and pd.isna(row["True_positive"])
    ):
        if row["False_positive"] > 2:
            return "Potential contamination"
        else:
            return "Sterile, background detected"
    # FNs
    elif (
        pd.isna(row["False_positive"])
        and pd.isna(row["True_negative"])
        and pd.isna(row["True_positive"])
    ):
        return "Contaminated"
    # TNs
    elif (
        pd.isna(row["False_positive"])
        and pd.isna(row["False_negative"])
        and pd.isna(row["True_positive"])
    ):
        return "Sterile"
    # TPs
    elif (
        pd.isna(row["False_positive"])
        and pd.isna(row["False_negative"])
        and pd.isna(row["True_negative"])
    ):
        return "Contaminated"
    # TNs, FPs
    elif pd.isna(row["False_negative"]) and pd.isna(row["True_positive"]):
        return "Sterile"
    # TPs, FPs
    elif pd.isna(row["False_negative"]) and pd.isna(row["True_negative"]):
        return "Potential contamination"
    # FNs, TNs
    elif pd.isna(row["False_positive"]) and pd.isna(row["True_positive"]):
        if row["False_negative"] - row["True_negative"] > -5:
            return "Potential contamination"
        else:
            return "Sterile, background detected"
    # FNs, TPs
    elif pd.isna(row["False_positive"]) and pd.isna(row["True_negative"]):
        return "Contaminated"
    # FNs, FPs, TNs
    elif (
        pd.isna(row["True_positive"])
        and not pd.isna(row["False_positive"])
        and not pd.isna(row["True_negative"])
        and not pd.isna(row["False_negative"])
    ):
        return "Sterile, background detected"
    # TPs, FPs, TNs
    elif (
        pd.isna(row["False_negative"])
        and not pd.isna(row["False_positive"])
        and not pd.isna(row["True_negative"])
        and not pd.isna(row["True_positive"])
    ):
        if row["True_negative"] > row["True_positive"] and row["False_positive"] > row["True_positive"]:
            return "Sterile"
        else:
            return "Contaminated"
    # TPs, FPs, FNs
    elif (
        pd.isna(row["True_negative"])
        and not pd.isna(row["False_positive"])
        and not pd.isna(row["False_negative"])
        and not pd.isna(row["True_positive"])
    ):
        return "Potential contamination"
    # TPs, FPs, FNs
    elif pd.isna(row["True_negative"]) and pd.isna(row["True_positive"]):
        return "Inconclusive"
    else:
        return "Unknown"
    
    
def check_final_decision(x):
    return True if x["Decision_mask"] == x["Decision"] else False


def check_unseen_decision(mask_col, not_mask_col, row):
    return True if row[mask_col] == row[not_mask_col] else False


def run_class_accuracy(df: pd.DataFrame,
                       metagenomic_cla: str,
                       col_correct: str,
                       model_type: str,
                       mask_col,
                       not_mask_col) -> None:
    partial_func = partial(check_unseen_decision, mask_col, not_mask_col)
    
    df[col_correct] = df.apply(partial_func, axis=1)
    acc = calculate_accuracy(df, col_correct, model_type, metagenomic_cla)
    calculate_classification_accuracy_interval(acc, df, model_type, metagenomic_cla)
    

def generate_decision_df(df: pd.DataFrame) -> pd.DataFrame:
    decision_df = (
        df.groupby(["sample", "Decision"]).size().reset_index()
    )
    decision_pivot = decision_df.pivot_table(0, ["sample"], "Decision")
    if "False_positive" not in decision_pivot.columns:
        decision_pivot["False_positive"] = np.nan
    if "False_negative" not in decision_pivot.columns:
        decision_pivot["False_negative"] = np.nan
    if "True_negative" not in decision_pivot.columns:
        decision_pivot["True_negative"] = np.nan
    if "True_positive" not in decision_pivot.columns:
        decision_pivot["True_positive"] = np.nan
        
    decision_pivot["Sterility"] = decision_pivot.apply(sample_decision, axis=1)
    decision_pivot.reset_index(inplace=True)
    
    return decision_pivot


def check_tts_status(row):
    if row["Dec_Correct"] == True and row["sample-status"] == "eval":
        return "TRUE-eval"
    if row["Dec_Correct"] == False and row["sample-status"] == "eval":
        return "FALSE-eval"
    if row["Dec_Correct"] == True and row["sample-status"] == "unknowns":
        return "TRUE-unknowns"
    if row["Dec_Correct"] == False and row["sample-status"] == "unknowns":
        return "FALSE-unknowns"


def get_cleaned_data(samp_df: pd.DataFrame, genu_df: pd.DataFrame) -> pd.DataFrame:
    req_columns=['sample', 'name', 'sample_true_mask', 'XGB_sample_prediction', 'sample-status']
    req_columns2=['gen-samp', 'gen-name', 'mask', 'XGB_gen_con_prediction', 'gen-con-pred-correct']
    genuine_status_ce = genu_df.rename(columns = {'name':'gen-name','sample':'gen-samp','XGB_genuine_prediction':'XGB_gen_con_prediction', 'Correct':'gen-con-pred-correct'})
    samp_clean = samp_df[req_columns]
    genu_clean = genuine_status_ce[req_columns2]
    model_compare_df = pd.concat([samp_clean, genu_clean], axis = 1)
    return model_compare_df


def calculate_twin_model_accuracy(df: pd.DataFrame, mg: str, _type_: str) -> (pd.DataFrame, pd.DataFrame):
    encoding_dict = {"True_negative": 0, "True_positive": 1}
    encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
    eval_gen_con_status = df.replace({"XGB_sample_prediction": encoding_dict_rev})
            
    eval_gen_con_status["Decision"] = eval_gen_con_status.apply(make_decision, axis=1)
    eval_gen_con_status["Decision_mask"] = eval_gen_con_status.apply(decision_mask, axis=1)
    eval_gen_con_status["Dec_Correct"] = eval_gen_con_status.apply(check_final_decision, axis=1)
    
    eval_gen_con_status["check-tts-status-decision"] = eval_gen_con_status.apply(check_tts_status, axis=1)
    
    cg_acc_dec = calculate_accuracy(eval_gen_con_status, "check-tts-status-decision", f"{_type_}-sample-genuine", mg)
    
    calculate_classification_accuracy_interval(cg_acc_dec, eval_gen_con_status, f"{_type_}-sample-genuine", mg)
    
    decision_pivot = generate_decision_df(eval_gen_con_status)
    
    decision_pivot["cfu"] = decision_pivot["sample"].str.split("_", expand = True)[3]
    
    return decision_pivot, eval_gen_con_status