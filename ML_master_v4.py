#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 11:35:31 2022

@author: mangi
"""

import pandas as pd
import numpy as np
import time
import json
import os
from collections.abc import Iterable
from typing import Tuple
from pipeline import xgboost_data_cleaning


################## PARAMETERS ##################

# DATA INPUT
bacteria_fungus_human = True

# USE FEATURE WIZ AND SAVE OUTPUTS
save = True
featurewiz_CBC = True
featurewiz_XGB = True
run_gridsearch = False

data_cleaning = True
feature_engineering = True
run_models = True

# DATABASE PARAMETERS
epoch_time = str(int(time.time()))
database_dict = {"v_f_b":"v_f_b"}
samples_to_subset = []
# github = "D:/GitHub/SMART-CAMP/"
github = "/mnt/d/GitHub/SMART-CAMP/"

# # SUBSET SAMPLES
subset = True
independent_var = "_aDNA_" # amplicon

# INDIVIDUAL EXPERIMENT INPUTS - SPIKE SPECIES, LABELS, TRAINING DATABASE, SAMPLES TO EXCLUDE FOR EVALUATION.

if bacteria_fungus_human:
    # directory = "D:/SequencingData/Harmonisation/DNA/analysis"
    directory = "/mnt/d/SequencingData/Harmonisation/DNA/analysis"
    json_mask = "catboost_decision_true_mask_bact_fungal-human-dictionary-824.json"
    species = {"Pseudomonas aeruginosa":["PA"], "Cutibacterium acnes":["Cacnes","Pacnes"], \
                "Escherichia coli":["EC"], "Klebsiella pneumoniae":["Klebpneu"], \
                  # "Candida albicans":["Calbicans"], "Staphylococcus aureus":["Saureus"], \
                 "Candida":["Calbicans"], "Staphylococcus aureus":["Saureus"], \
                   "Bacillus subtilis": ["Bsubtilis"]}

_ML_out_ = f"{directory}/ML_training-VFB/"
out_OCS = f"{_ML_out_}/OneClassSVM-VFB/"
cat_out = f"{_ML_out_}/OCS-catboost-VFB/"
XGB_out = f"{_ML_out_}/OCS-XGBoost-VFB/"
os.makedirs(_ML_out_, exist_ok=True)
os.makedirs(cat_out, exist_ok=True) 
os.makedirs(XGB_out, exist_ok=True)   
os.makedirs(out_OCS, exist_ok=True)

BLASTn_name = list(database_dict.keys())[0]
kingdom = list(database_dict.values())[0]

################## PARAMETERS ##################
#
#
#
################## DATA CLEANING ##################

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el

vs = []
for k,v in species.items():
    vs.append(v)
vs = list(flatten(vs))
vees = "-".join(vs)


hd_data_cols, \
cd_data_cols, \
    cn_df, \
    bn_df = \
        xgboost_data_cleaning.main(
        directory,
        BLASTn_name,
        kingdom,
        species, 
        _ML_out_,
        subset,
        independent_var)
        
################## DATA CLEANING ##################
#
#
#
################## FEATURE ENGINEERING / SELECTION ##################

# add mask for TP & TN vs other
def apply_mask(df: pd.DataFrame, species: dict) -> pd.DataFrame:
    df["strain"] = df["sample"].str.split("_", expand=True)[2]

    true_masks = []
    for k, vv in species.items():
        if isinstance(vv, list):
            for v in vv:
                true_mask = df.loc[
                    (df.strain.str.contains(v)) & (df.name.str.contains(k))
                ].index
                df.loc[true_mask, "mask"] = "True_positive"
                true_masks.append(true_mask)
        else:
            true_mask = df.loc[
                (df.strain.str.contains(vv)) & (df.name.str.contains(k))
            ].index
            df.loc[true_mask, "mask"] = "True_positive"
            true_masks.append(true_mask)
    true_masks = list(set([item for sublist in true_masks for item in sublist]))

    false_mask = list(set(list(df.index)) - set(true_masks))
    df.loc[false_mask, "mask"] = "False_positive"
    df.loc[true_masks, "mask"] = "True_positive"
    return df


masked_ce = apply_mask(cn_df, species)
masked_bl = apply_mask(bn_df, species)

################## FEATURE ENGINEERING / SELECTION ##################
#
#
#
################## DATA SPLIT - TRAIN + TEST / EVALUATION ##################

from sklearn.model_selection import train_test_split

# SPLIT DATA INTO TRAIN, TEST, EVALUATION
def split_data(df: pd.DataFrame, columns: list) -> Tuple[pd.DataFrame, pd.DataFrame]:
    X = df[columns].values
    indices = df.index.values
    
    train_test, evaluation, indices_train_test, indices_evaluation = train_test_split(
        X, indices, test_size=0.25, random_state=736
    )
    
    train_test_df = df.iloc[indices_train_test]
    evaluation_df = df.iloc[indices_evaluation]
    return train_test_df, evaluation_df

train_test_df_ce, evaluation_df_ce = split_data(masked_ce, cd_data_cols)
train_test_df_bn, evaluation_df_bn = split_data(masked_bl, hd_data_cols)

print(f"Split counts for training, centrifuge: \n{train_test_df_ce['mask'].value_counts()}")
print(f"Split counts for evaluation, centrifuge: \n{evaluation_df_ce['mask'].value_counts()}")
print(f"Split counts for training, BLAST: \n{train_test_df_bn['mask'].value_counts()}")
print(f"Split counts for evaluation, BLAST: \n{evaluation_df_bn['mask'].value_counts()}")

################## DATA SPLIT - TRAIN + TEST / EVALUATION ##################
#
#
#
################## MACHINE LEARNING ##################
################## TRAINING + TESTING ##################

from functools import partial
from featurewiz import featurewiz
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from string import whitespace
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import statistics
from collections import Counter
from sklearn.metrics import classification_report, confusion_matrix
from pathlib import Path
from sklearn.metrics import accuracy_score, recall_score, precision_score
from matplotlib import pyplot
    

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


def prepare_data_sample(true_mask: dict, 
                        input_df: pd.DataFrame,
                        data_cols: list
                        ) -> Tuple[pd.DataFrame, list, 
                              pd.DataFrame, list, 
                              pd.DataFrame, dict,
                              pd.DataFrame]:

    def custom_label_encoder(encoding_dict: dict, row: int):
        for name, encoded in encoding_dict.items():
            if name in row["sample_true_mask"]:
                return encoded
    
    
    def one_hot_encoder(encoding_dict: dict, row: int):
        for name, encoded in encoding_dict.items():
            if encoded == row["sample_encoded_mask"]:
                if encoded == 0:
                    return (1, 0)
                if encoded == 1:
                    return (0, 1)
                else:
                    return (0, 0)                              
                              
    true_masked_dfs = []
    for k, v in true_mask.items():
        true_mask_df = pd.DataFrame.from_dict(
            v, columns=["sample_true_mask"], orient="index"
        )
        true_mask_df["db"] = k
        true_masked_dfs.append(true_mask_df)
    cat_true_mask_df = pd.concat(true_masked_dfs)
    cat_true_mask_df.reset_index(inplace=True)
    cat_true_mask_df = cat_true_mask_df.rename(columns={"index": "sample"})
    cat_true_mask_df = cat_true_mask_df.replace("CFU", "", regex=True)
    true_mask_df_merge = pd.merge(
        
        ###############################################
        input_df,
        ###############################################
        
        cat_true_mask_df,
        how="left",
        left_on=["sample", "db"],
        right_on=["sample", "db"],
    )
    
    true_mask_df_merge = true_mask_df_merge.drop(
        true_mask_df_merge.index[
            true_mask_df_merge.loc[
                true_mask_df_merge["sample_true_mask"] == "False_positive"
            ].index
        ]
    )
    true_mask_df_merge.reset_index(inplace=True, drop=True)
    
    encoding_dict = {"True_negative": 0, "True_positive": 1}
    target_columns = ["sample-mask-TP", "sample-mask-TN"]
    
    encoded_ML = true_mask_df_merge[data_cols + ["sample_true_mask"]]
    
    label_partial_func = partial(custom_label_encoder, encoding_dict)
    encoded_ML["sample_encoded_mask"] = encoded_ML.apply(label_partial_func, axis=1)
    
    enc_partial_func = partial(one_hot_encoder, encoding_dict)
    encoded_ML["temp-mask"] = encoded_ML.apply(enc_partial_func, axis=1)
    
    encoded_ML[target_columns] = pd.DataFrame(
        encoded_ML["temp-mask"].tolist(), index=encoded_ML.index
    )
    
    encoded_ML = encoded_ML.drop(["temp-mask"], axis=1)
    
    all_masks = [
        "sample_true_mask",
        "sample_encoded_mask",
        "sample-mask-TN",
        "sample-mask-TP",
    ]
    
    mask_to_drop = list(set(all_masks) - set(target_columns))
    other_mask_df = encoded_ML[mask_to_drop]
    encoded_ML = encoded_ML.drop(mask_to_drop, axis=1)
    
    return encoded_ML, target_columns, true_mask_df_merge, all_masks, other_mask_df, encoding_dict, true_mask_df
    

def check_prediction(x):
    return True if x["sample_true_mask"] == x["XGB_sample_prediction"] else False


def draw_heatmap(fw_heatmap_df: pd.DataFrame, 
                 XGB_out: str,
                 dataset_used: str,
                 metagenome_classifier_used: str,
                 title: str) -> None:
    featurecorr = fw_heatmap_df.corr()
    featurecorr.index = featurecorr.index.str.replace(r"_", " ")
    featurecorr.columns = featurecorr.columns.str.replace(r"_", " ")
    
    mask = np.zeros_like(featurecorr)
    mask[np.triu_indices_from(mask)] = True
    
    f, ax = plt.subplots(figsize=(15, 15))
    ax = sns.heatmap(featurecorr, mask=mask)

    new_title = clean_strings(title, 120)
    plt.title(new_title, size=40)
    plt.tight_layout()
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    # use matplotlib.colorbar.Colorbar object
    cbar = ax.collections[0].colorbar
    # here set the labelsize by 20
    cbar.ax.tick_params(labelsize=30)
    plt.savefig(f"{XGB_out}/{title}-decision.png", dpi=300, bbox_inches="tight")
    plt.show()
              

def check_mask_correct_pred(x):
    return True if x["mask"] == x["XGB_gen_con_prediction"] else False


def custom_samp_label_encoder(encoding_dict: dict, row: int):
    for name, encoded in encoding_dict.items():
        if name in row["XGB_sample_prediction"]:
            return encoded


def one_hot_samp_encoder(encoding_dict: dict, row: int):
    for name, encoded in encoding_dict.items():
        if encoded == row["XGB_sample_prediction"]:
            if encoded == 0:
                return (1, 0)
            if encoded == 1:
                return (0, 1)
            else:
                return (0, 0)
                
 
def custom_contaminant_label_encoder(encoding_dict: dict, row: int):
    for name, encoded in encoding_dict.items():
        if name in row["mask"]:
            return encoded


def one_hot_contaminant_encoder(encoding_dict: dict, row: int):
    for name, encoded in encoding_dict.items():
        if encoded == row["encoded_mask"]:
            if encoded == 0:
                return (1, 0)
            if encoded == 1:
                return (0, 1)
            else:
                return (0, 0)
                 

pd.options.mode.chained_assignment = None  # default='warn'

json_file = f"{cat_out}{json_mask}"

with open(json_file) as json_data:
    true_mask = json.loads(json_data.read())
    
ce_id_cols = ["name", "sample", "seqID"]
bl_id_cols = ["name", "sample", "sseqid"]

################## MACHINE LEARNING ##################
################## TRAINING + TESTING ##################

train_test_df_ce = train_test_df_ce.drop_duplicates(subset=ce_id_cols)
train_test_df_ce.reset_index(inplace=True, drop=True)
train_test_df_bn = train_test_df_bn.drop_duplicates(subset=bl_id_cols)
train_test_df_bn.reset_index(inplace=True, drop=True)

train_test_df_bn.to_csv(f"{_ML_out_}/test-hs-dataset.csv", index=False)


def evaluate_model_performance(y_pred, 
                               XGB_classifier_model, 
                               X_test, 
                               X_train, 
                               y_train, 
                               y_test, 
                               metagenome_classifier_used: str, 
                               ml_model: str,
                               _type_: str, 
                               encoding_dict: dict):    
    # 5 folds, scored on accuracy
    # https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter
    if _type_ == "test":
        cvs1 = cross_val_score(
            XGB_classifier_model, X_train, y_train, cv=10, scoring="accuracy"
        )
    cvs2 = cross_val_score(
        XGB_classifier_model, X_test, y_test, cv=10, scoring="accuracy"
    )
    
    e_d_list = list({f"{k}-{v}" for k,v in encoding_dict.items()})
    
    print(
        "-------------------------------------------------------------------------------"
    )
    print(f"XGBoost - {ml_model} status - {_type_} - {metagenome_classifier_used}\nEncoding dictionary: {e_d_list[0]}, {e_d_list[1]}")
    print("  ")
    confusion = confusion_matrix(y_test, y_pred)
    TP = confusion[1, 1]
    TN = confusion[0, 0]
    FP = confusion[0, 1]
    FN = confusion[1, 0]
    
    # Classification Accuracy
    print("Classification Accuracy:", f"{(TP + TN) / float(TP + TN + FP + FN):.2f}")
    # print("Classification Accuracy:", f"{accuracy_score(y_test, y_pred):.2f}")
    # Classification Error
    classification_error = (FP + FN) / float(TP + TN + FP + FN)
    print("Classification Error:", f"{classification_error:.2f}")
    # print("Classification Error:", f"{1 - accuracy_score(y_test, y_pred):.2f}")    
    # Sensitivity
    sensitivity = TP / float(FN + TP)
    print(f"Sensitivity: {sensitivity:.2f}")
    # print(f"Sensitivity: {recall_score(y_test, y_pred):.2f}")
    # Specificity
    specificity = TN / (TN + FP)
    print(f"Specificity: {specificity:.2f}")
    # False Positive Rate
    false_positive_rate = FP / float(TN + FP)
    print(f"FPR: {false_positive_rate:.2f}")
    # print(f"FPR: {1 - specificity:.2f}")
    # Precision
    precision = TP / float(TP + FP)
    print(f"Precision: {precision:.2f}")
    # print(f"Precision: {precision_score(y_test, y_pred):.2f}")   
    
    jaccard_index = TP / (TP + FP + TN)
    print(f"Jaccard index / critical success index: {jaccard_index:.2f}")
    
    print(confusion)
    print(classification_report(y_test, y_pred))  # Output
    print(
        "-------------------------------------------------------------------------------"
    )

    if _type_ == "test":    
        print(f"Cross-val #1 mean: {statistics.mean(cvs1):.2f}")
        print(f"Cross-val #1 stdev: {statistics.stdev(cvs1):.2f}")
    
    print(f"Cross-val #2 mean: {statistics.mean(cvs2):.2f}")
    print(f"Cross-val #2 stdev: {statistics.stdev(cvs2):.2f}")
    

def plot_error(dataset_used: str, 
               metagenome_classifier_used: str, 
               y_pred: np.ndarray, 
               y_test: np.ndarray, 
               XGB_classifier_model: XGBClassifier): 

    predictions = [round(value) for value in y_pred]
    accuracy = accuracy_score(y_test, predictions)
    print(f"\nAccuracy score for test sample: {accuracy:.2f}\n")
    results = XGB_classifier_model.evals_result()
    epochs = len(results["validation_0"]["error"])
    x_axis = range(0, epochs)

    # plot log loss
    fig, ax = pyplot.subplots(figsize=(15,15))
    ax.plot(x_axis, results["validation_0"]["logloss"], label="Train", linestyle="--",linewidth=7.0)
    ax.plot(x_axis, results["validation_1"]["logloss"], label="Test", linestyle="--",linewidth=7.0)
    ax.legend()
    pyplot.xlabel("Epochs", size=40)
    pyplot.ylabel("Log Loss", size=40)
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    title = f"XGBoost Log Loss for sample contaminant status for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used}"
    new_title = clean_strings(title, 120)
    pyplot.title(new_title, size=40)
    plt.legend(loc=0, prop={"size": 40}, markerscale=10)
    pyplot.show()

    # plot classification error = (FP + FN) / float(TP + TN + FP + FN)
    fig, ax = pyplot.subplots(figsize=(15,15))
    ax.plot(x_axis, results["validation_0"]["error"], label="Train", linestyle="--",linewidth=7.0)
    ax.plot(x_axis, results["validation_1"]["error"], label="Test", linestyle="--",linewidth=7.0)
    ax.legend()
    pyplot.xlabel("Epochs", size=40)
    pyplot.ylabel("Classification Error", size=40)
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    title = f"XGBoost Classification Error for sample contaminant status for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used}"
    new_title = clean_strings(title, 120)
    pyplot.title(new_title, size=40)
    plt.legend(loc=0, prop={"size": 40}, markerscale=10)
    pyplot.show()
 
    # plot auc
    fig, ax = pyplot.subplots(figsize=(15,15))
    ax.plot(x_axis, results["validation_0"]["auc"], label="Train", linestyle="--",linewidth=7.0)
    ax.plot(x_axis, results["validation_1"]["auc"], label="Test", linestyle="--",linewidth=7.0)
    ax.legend()
    pyplot.xlabel("Epochs", size=40)
    pyplot.ylabel("True Positive Rate", size=40)
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    title = f"AUC for sample contaminant status for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used}"
    new_title = clean_strings(title, 120)
    pyplot.title(new_title, size=40)
    plt.legend(loc=0, prop={"size": 40}, markerscale=10)
    pyplot.show()
    
    
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

def run_xgboost_sample_status(metagenome_classifier_used: str,
                              dataset_used: str,
                              true_mask: dict,
                              input_df: pd.DataFrame,
                              data_cols: list, 
                              id_cols: list, 
                              vees: str, 
                              XGB_out: str,
                              save: bool,
                              run_gridsearch: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    filename_save = f"xgboost_sample_status-{metagenome_classifier_used}.sav"
    
    encoded_ML, target_columns, true_mask_df_merge, all_masks, \
        other_mask_df, encoding_dict, true_mask_df = \
        prepare_data_sample(true_mask, input_df, data_cols)
    
    # find optimal features
    features = featurewiz(encoded_ML, target=target_columns, corr_limit=0.70, verbose=2)
    feature_names, feature_df = features
    fw_heatmap_df = feature_df
    human_readable = feature_df
    
    title = f"XGBoost-sample-contaminant-status -- Heat map depicting important features for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used} using featurewiz"
    draw_heatmap(fw_heatmap_df, 
                     XGB_out,
                     dataset_used,
                     metagenome_classifier_used,
                     title)

    feat_labels = feature_names
    
    with open(
        f"{XGB_out}/feature-wiz-xgboost-features-{vees}-{metagenome_classifier_used}-sample-status.json",
        "w", encoding="utf-8") as f:
        json.dump(feat_labels, f)
    
    instance_names_df = true_mask_df_merge[id_cols]
    
    feature_wiz_df = pd.merge(
        instance_names_df, human_readable, how="left", left_index=True, right_index=True
    )
    feature_wiz_df = pd.merge(
        feature_wiz_df, other_mask_df, how="left", left_index=True, right_index=True
    )
    feature_wiz_df = feature_wiz_df.drop_duplicates(subset=id_cols)
    feature_wiz_df.reset_index(inplace=True, drop=True)
    
    ######################################
    y = feature_wiz_df["sample_encoded_mask"].values  
    indices = feature_wiz_df.index.values
    feature_wiz_clean = feature_wiz_df.drop(
        all_masks + id_cols, axis=1
    )
    X = feature_wiz_clean.iloc[:, :].values
    
    true_mask_df_merge = true_mask_df_merge.drop_duplicates(subset=id_cols)
    true_mask_df_merge.reset_index(inplace=True, drop=True)
    
    X_train, X_test, indices_train, indices_test = train_test_split(
        X, indices, test_size=0.25, random_state=736
    )
    y_train, y_test = y[indices_train], y[indices_test]
    ######################################
    
    sc = StandardScaler()
    
    training_X = sc.fit_transform(X_train)
    
    test_X = sc.transform(X_test)

    eval_metric = ["logloss", "auc", "error"]
    eval_set = [(training_X, y_train), (test_X, y_test)]
    
    ######################################
    XGB_classifier_model = XGBClassifier(
        verbosity=0, random_state=736, alpha=0, gamma=0, max_depth=5, 
        subsample=0.5, colsample_by_tree = 0.5, objective = "binary:logistic",
    )
    
    XGB_classifier_model.fit(training_X, y_train, eval_metric=eval_metric, eval_set=eval_set, verbose=False)
    ######################################
    
    # conformal_XGB = copy.deepcopy(XGB_classifier_model)
    print("\n\nConformal analysis")
    # conformal_analysis(dataset_used, metagenome_classifier_used, test_X, y_test, XGB_classifier_model, training_X, y_train, "sample status", XGB_out)
    
    y_pred = XGB_classifier_model.predict(test_X)

    print(f"\nAccuracy score for test sample: {accuracy_score(y_test, y_pred):.2f}\n")
    
    plot_error(dataset_used, metagenome_classifier_used, y_pred, y_test, XGB_classifier_model)
    
    # save the scaler
    if save:
        pickle.dump(sc, open(f"{XGB_out}/sc-sample-status-XGB-scaler-{metagenome_classifier_used}.pkl", "wb"))
    
    ############ GRID SEARCH ##############
    
    # Method 1
    if run_gridsearch:
        grid = dict()
        grid["silent"] = [1]
        grid["gamma"] = [0.001, 0.01, 0.1, 0]
        grid["max_depth"] = [4, 5, 6]
        grid["lambda"] = [0, 1, 5]
        grid["alpha"] = [0, 1, 5]
    
        # Instantiate GridSearchCV
        xgb = XGBClassifier(verbosity=0, random_state=736)
        gscv = GridSearchCV(estimator=xgb, param_grid=grid, scoring="accuracy", cv=5)
    
        # fit the model
        gscv.fit(training_X, y_train)
    
        # returns the estimator with the best performance
        print(gscv.best_estimator_)
    
        # returns the best score
        print(gscv.best_score_)
    
        # returns the best parameters
        print(gscv.best_params_)
    ############ GRID SEARCH ##############
    
    # Get predicted probabilities for each class
    y_preds_proba = XGB_classifier_model.predict_proba(test_X)
    y_preds_train = XGB_classifier_model.predict_proba(training_X)
    
    # roc curve for models
    fpr1, tpr1, thresh1 = roc_curve(y_test, y_preds_proba[:, 1], pos_label=1)
    fpr2, tpr2, thresh2 = roc_curve(y_train, y_preds_train[:, 1], pos_label=1)
    
    random_probs = [0 for i in range(len(y_test))]
    p_fpr, p_tpr, _ = roc_curve(y_test, random_probs, pos_label=1)
    
    # auc scores
    auc_score1 = roc_auc_score(y_test, y_preds_proba[:, 1])
    auc_score2 = roc_auc_score(y_train, y_preds_train[:, 1])
    
    plotting_data = [
        auc_score1,
        auc_score2,
        fpr1.tolist(),
        tpr1.tolist(),
        fpr2.tolist(),
        tpr2.tolist(),
        ]
    
    if save:
        with open(f"{XGB_out}/ROC-sample-status-XGB-plotting-data-{metagenome_classifier_used}.json", "w", encoding="utf-8") as f:
            json.dump(plotting_data, f)
    
    plt.style.use("seaborn")
    f, ax = plt.subplots(figsize=(15, 15))
    # plot roc curves
    plt.plot(fpr1, tpr1, linestyle="--", color="red", label="XGB-test")
    plt.plot(fpr2, tpr2, linestyle="--", color="orange", label="XGB-train")
    # plt.plot(p_fpr, p_tpr, linestyle="--", color="blue")
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    # title
    title = f"AUC-ROC curve for sample contaminant status for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used}"
    new_title = clean_strings(title, 120)
    plt.title(new_title, size=40)
    # x label
    plt.xlabel("False Positive Rate", size=40)
    # y label
    plt.ylabel("True Positive rate", size=40)
    plt.text(
        0.8,
        0.2,
        f"AUC-test: {auc_score1:.2f}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        size=30,
    )
    plt.text(
        0.8,
        0.25,
        f"AUC-train: {auc_score2:.2f}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        size=30,
    )
    plt.legend(loc=4, prop={"size": 40}, markerscale=10)
    if save:
        plt.savefig(f"{XGB_out}/{title}.png", dpi=300)
    plt.show()
    
    
    XGB_model_save = f"{XGB_out}/{filename_save}"
    print(f"Saving one XGBoost classifier model to: {XGB_model_save}")
    if save:
        pickle.dump(XGB_classifier_model, open(XGB_model_save, "wb"))


    evaluate_model_performance(y_pred, XGB_classifier_model, 
                                   X_test, X_train, y_train, 
                                   y_test, metagenome_classifier_used,
                                   "sample contaminant",
                                   "test", encoding_dict)
        
    # to get convert dictionary to list
    # count number of true labels and samples
    t_positives = 0
    no_samples = len(true_mask_df)
    for db, dict_ in true_mask.items():
        words = list(true_mask[db].values())
        keyss = Counter(words).keys()  # equals to list(set(words)))
        valuess = Counter(words).values()  # counts the elements' frequency)
        index_key = 0
        if "True_positive" in keyss:
            index_key = list(keyss).index("True_positive")
            t_positive = list(valuess)[index_key]
            if t_positive > 0:
                t_positives = t_positive
    
    print(
        f"Positive labels: {t_positives}/{no_samples}; {t_positives/no_samples*100:.2f} %"
        )
    
    sample_status = true_mask_df_merge[
        ["sample", "name", "mask", "sample_true_mask"]
    ]
    sample_status.loc[indices_test, "XGB_sample_prediction"] = XGB_classifier_model.predict(
        test_X
    )
    sample_status.loc[
        indices_train, "XGB_sample_prediction"
    ] = XGB_classifier_model.predict(training_X)
    
    sample_status.loc[indices_test, "sample-status"] = "test"
    sample_status.loc[indices_train, "sample-status"] = "train"
    
    encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
    sample_status = sample_status.replace({"XGB_sample_prediction": encoding_dict_rev})
    
    sample_status["sample-pred-correct"] = sample_status.apply(check_prediction, axis=1)
    
    return sample_status, true_mask_df
   

ce_tt_sample_status, true_mask_df = run_xgboost_sample_status("centrifuge", "train-test", true_mask, train_test_df_ce, cd_data_cols, ce_id_cols, vees, XGB_out, save, run_gridsearch)
hs_tt_sample_status, _ = run_xgboost_sample_status("BLAST", "train-test", true_mask, train_test_df_bn, hd_data_cols, bl_id_cols, vees, XGB_out, save, run_gridsearch)

                
def run_xgboost_genuine_contaminant_status(metagenome_classifier_used: str,
                              dataset_used: str,
                              true_mask: dict,
                              input_df: pd.DataFrame,
                              input_df_xgboost: pd.DataFrame,
                              data_cols: list, 
                              id_cols: list, 
                              vees: str,
                              true_mask_df: pd.DataFrame, 
                              XGB_out: str,
                              save: bool,
                              run_gridsearch: bool) -> pd.DataFrame:

    filename_save = f"xgboost_genuine_contaminant_status-{metagenome_classifier_used}.sav"
    
    genuine_contaminant_df = pd.merge(input_df, input_df_xgboost[['sample_true_mask', 'XGB_sample_prediction', 'sample-status',
           'sample-pred-correct']], left_index=True, right_index=True)
    
    # convert prediction T/F into encoded 0,1
    encoding_dict_samp = {"True_positive": 1, "True_negative": 0}
    target_columns_samp = ["sample-pred-TP", "sample-pred-TN"]
    label_partial_func = partial(custom_samp_label_encoder, encoding_dict_samp)
    genuine_contaminant_df["XGB_sample_prediction"] = genuine_contaminant_df.apply(label_partial_func, axis=1)
    
    enc_samp_partial_func = partial(one_hot_samp_encoder, encoding_dict_samp)
    genuine_contaminant_df["temp-mask"] = genuine_contaminant_df.apply(enc_samp_partial_func, axis=1)    
    
    genuine_contaminant_df[target_columns_samp] = pd.DataFrame(
        genuine_contaminant_df["temp-mask"].tolist(), index=genuine_contaminant_df.index
    )
    
    genuine_contaminant_df = genuine_contaminant_df.drop(["temp-mask"], axis=1)
    
    encoded_ML = genuine_contaminant_df[data_cols + ["mask"] + target_columns_samp]

    encoding_dict = {"False_positive": 0, "True_positive": 1}
    target_columns = ["mask-FP", "mask-TP"]
    
    label_partial_func = partial(custom_contaminant_label_encoder, encoding_dict)
    encoded_ML["encoded_mask"] = encoded_ML.apply(label_partial_func, axis=1)
    
    enc_partial_func = partial(one_hot_contaminant_encoder, encoding_dict)
    encoded_ML["temp-mask"] = encoded_ML.apply(enc_partial_func, axis=1)
    encoded_ML[target_columns] = pd.DataFrame(
        encoded_ML["temp-mask"].tolist(), index=encoded_ML.index
    )
    encoded_ML = encoded_ML.drop(["temp-mask"], axis=1)
    
    all_masks = ["mask", "encoded_mask", "mask-FP", "mask-TP"]
    
    mask_to_drop = list(set(all_masks) - set(target_columns))
    other_mask_df = encoded_ML[mask_to_drop]
    encoded_ML = encoded_ML.drop(mask_to_drop, axis=1)
    
    # find optimal features
    features = featurewiz(encoded_ML, target=target_columns, corr_limit=0.70, verbose=2)
    feature_names, feature_df = features
    fw_heatmap_df = feature_df
    human_readable = feature_df
    
    title = f"XGBoost-genuine-contaminant-status -- Heat map depicting important features for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used} using featurewiz"
    draw_heatmap(fw_heatmap_df, 
                     XGB_out,
                     dataset_used,
                     metagenome_classifier_used,
                     title)
    
    feat_labels = feature_names
    if save:
        with open(
            f"{XGB_out}/feature-wiz-xgboost-features-{vees}-genuine-contaminant-{metagenome_classifier_used}.json",
            "w", encoding="utf-8") as f:
            json.dump(feat_labels, f)

    instance_names_df = genuine_contaminant_df[id_cols]
    
    feature_wiz_df = pd.merge(
        instance_names_df, human_readable, how="left", left_index=True, right_index=True
    )
    feature_wiz_df = pd.merge(
        feature_wiz_df, other_mask_df, how="left", left_index=True, right_index=True
    )
    feature_wiz_df = feature_wiz_df.drop_duplicates(subset=id_cols)
    feature_wiz_df.reset_index(inplace=True, drop=True)
    
    y = feature_wiz_df["encoded_mask"].values  
    indices = feature_wiz_df.index.values
    feature_wiz_clean = feature_wiz_df.drop(
        all_masks + id_cols, axis=1
    )
    X = feature_wiz_clean.iloc[:, :].values
    
    instance_names_df = instance_names_df.drop_duplicates(subset=id_cols)
    instance_names_df.reset_index(inplace=True, drop=True)
    
    X_train, X_test, indices_train, indices_test = train_test_split(
        X, indices, test_size=0.25, random_state=736
    )
    y_train, y_test = y[indices_train], y[indices_test]
    
    sc = StandardScaler()
    
    training_X = sc.fit_transform(X_train)
    
    test_X = sc.transform(X_test)
    
    XGB_classifier_model = XGBClassifier(
        verbosity=0, random_state=736, alpha=0, gamma=0, max_depth=5, subsample=0.5, colsample_by_tree = 0.5
    )
    
    XGB_classifier_model.fit(training_X, y_train)
    
    print("\n\nConformal analysis")
    # conformal_analysis(dataset_used, metagenome_classifier_used, test_X, y_test, XGB_classifier_model, training_X, y_train, "genuine status", XGB_out)
    
    y_pred = XGB_classifier_model.predict(test_X)

    print(f"\nAccuracy score for test sample: {accuracy_score(y_test, y_pred):.2f}\n")

    # save the scaler
    if save:
        pickle.dump(sc, open(f"{XGB_out}/sc-genuine-contaminant-XGB-scaler-{metagenome_classifier_used}.pkl", "wb"))
    
    ############ GRID SEARCH ##############
    
    # Method 1
    if run_gridsearch:
        grid = dict()
        grid["silent"] = [1]
        grid["gamma"] = [0.001, 0.01, 0.1, 0]
        grid["max_depth"] = [4, 5, 6]
        grid["lambda"] = [0, 1, 5]
        grid["alpha"] = [0, 1, 5]
    
        # Instantiate GridSearchCV
        xgb = XGBClassifier(verbosity=0, random_state=736)
        gscv = GridSearchCV(estimator=xgb, param_grid=grid, scoring="accuracy", cv=5)
    
        # fit the model
        gscv.fit(training_X, y_train)
    
        # returns the estimator with the best performance
        print(gscv.best_estimator_)
    
        # returns the best score
        print(gscv.best_score_)
    
        # returns the best parameters
        print(gscv.best_params_)
    ############ GRID SEARCH ##############
    
    y_preds_proba = XGB_classifier_model.predict_proba(test_X)
    y_preds_train = XGB_classifier_model.predict_proba(training_X)
    
    # roc curve for models
    fpr1, tpr1, thresh1 = roc_curve(y_test, y_preds_proba[:, 1], pos_label=1)
    fpr2, tpr2, thresh2 = roc_curve(y_train, y_preds_train[:, 1], pos_label=1)
    
    random_probs = [0 for i in range(len(y_test))]
    p_fpr, p_tpr, _ = roc_curve(y_test, random_probs, pos_label=1)
    
    # auc scores
    auc_score1 = roc_auc_score(y_test, y_preds_proba[:, 1])
    auc_score2 = roc_auc_score(y_train, y_preds_train[:, 1])
    
    plotting_data = [
        auc_score1,
        auc_score2,
        fpr1.tolist(),
        tpr1.tolist(),
        fpr2.tolist(),
        tpr2.tolist(),
        ]
    
    if save:
        with open(f"{XGB_out}/ROC-adventitious-agent-XGB-plotting-data-{metagenome_classifier_used}.json", "w", encoding="utf-8") as f:
            json.dump(plotting_data, f)
    
    plt.style.use("seaborn")
    f, ax = plt.subplots(figsize=(15, 15))
    # plot roc curves
    plt.plot(fpr1, tpr1, linestyle="--", color="red", label="XGB-test")
    plt.plot(fpr2, tpr2, linestyle="--", color="orange", label="XGB-train")
    # plt.plot(p_fpr, p_tpr, linestyle="--", color="blue")
    ax.tick_params(axis="x", labelsize=25)
    ax.tick_params(axis="y", labelsize=25)
    # title
    title = f"AUC-ROC curve for genuine contaminant status for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used}"
    new_title = clean_strings(title, 120)
    plt.title(new_title, size=40)
    # x label
    plt.xlabel("False Positive Rate", size=40)
    # y label
    plt.ylabel("True Positive rate", size=40)
    plt.text(
        0.8,
        0.2,
        f"AUC-test: {auc_score1:.2f}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        size=30,
    )
    plt.text(
        0.8,
        0.25,
        f"AUC-train: {auc_score2:.2f}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        size=30,
    )
    plt.legend(loc=4, prop={"size": 40}, markerscale=10)
    if save:
        plt.savefig(f"{XGB_out}/{title}.png", dpi=300)
    plt.show()
    
    
    XGB_model_save = f"{XGB_out}/{filename_save}"
    print(f"Saving one XGBoost classifier model to: {XGB_model_save}")
    if save:
        pickle.dump(XGB_classifier_model, open(XGB_model_save, "wb"))

    evaluate_model_performance(y_pred, XGB_classifier_model, 
                                       X_test, X_train, y_train, 
                                       y_test, metagenome_classifier_used,
                                       "genuine contaminant",
                                       "test",
                                       encoding_dict)

    # to get convert dictionary to list
    # count number of true labels and samples
    t_positives = 0
    no_samples = len(true_mask_df)
    for db, dict_ in true_mask.items():
        words = list(true_mask[db].values())
        keyss = Counter(words).keys()  # equals to list(set(words)))
        valuess = Counter(words).values()  # counts the elements' frequency)
        index_key = 0
        if "True_positive" in keyss:
            index_key = list(keyss).index("True_positive")
            t_positive = list(valuess)[index_key]
            if t_positive > 0:
                t_positives = t_positive
    
    print(
        f"Positive labels: {t_positives}/{no_samples}; {t_positives/no_samples*100:.2f} %"
        )
    
    gen_con_status = genuine_contaminant_df[
        ["sample", "name", "mask", "sample_true_mask", "XGB_sample_prediction"]
    ]
    gen_con_status.loc[indices_test, "XGB_gen_con_prediction"] = XGB_classifier_model.predict(
        test_X
    )
    gen_con_status.loc[
        indices_train, "XGB_gen_con_prediction"
    ] = XGB_classifier_model.predict(training_X)
    
    gen_con_status.loc[indices_test, "sample-status"] = "test"
    gen_con_status.loc[indices_train, "sample-status"] = "train"
    
    encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
    gen_con_status = gen_con_status.replace({"XGB_gen_con_prediction": encoding_dict_rev})
    
    gen_con_status["gen-con-pred-correct"] = gen_con_status.apply(check_mask_correct_pred, axis=1)
    
    return gen_con_status


ce_tt_gen_con_status = run_xgboost_genuine_contaminant_status("centrifuge", "train-test", true_mask, train_test_df_ce, ce_tt_sample_status, cd_data_cols, ce_id_cols, vees, true_mask_df, XGB_out, save, run_gridsearch)
hs_tt_gen_con_status = run_xgboost_genuine_contaminant_status("BLAST", "train-test", true_mask, train_test_df_bn, hs_tt_sample_status, hd_data_cols, bl_id_cols, vees, true_mask_df, XGB_out, save, run_gridsearch)

################## TRAINING + TESTING ##################
################## MACHINE LEARNING ##################
#
#
#
################## MACHINE LEARNING ##################
################## EVALUATION ##################

evaluation_df_ce = evaluation_df_ce.drop_duplicates(subset=ce_id_cols)
evaluation_df_ce.reset_index(inplace=True, drop=True)
evaluation_df_bn = evaluation_df_bn.drop_duplicates(subset=bl_id_cols)
evaluation_df_bn.reset_index(inplace=True, drop=True)


def run_xgboost_sample_status_evaluate(metagenome_classifier_used: str,
                              dataset_used: str,
                              true_mask: dict,
                              input_df: pd.DataFrame,
                              data_cols: list, 
                              id_cols: list, 
                              vees: str, 
                              XGB_out: str,
                              save: bool,
                              run_gridsearch: bool) -> pd.DataFrame:

    encoded_ML, target_columns, true_mask_df_merge, all_masks, \
        other_mask_df, encoding_dict, true_mask_df = \
        prepare_data_sample(true_mask, input_df, data_cols)
        
    with open(
        f"{XGB_out}/feature-wiz-xgboost-features-{vees}-{metagenome_classifier_used}-sample-status.json"
    ) as json_sample_XGB:
        feat_labels = json.loads(json_sample_XGB.read())
        
    sc = pickle.load(open(f"{XGB_out}/sc-sample-status-XGB-scaler-{metagenome_classifier_used}.pkl", "rb"))
    
    filename_save = f"xgboost_sample_status-{metagenome_classifier_used}.sav"
    XGB_model_save = f"{XGB_out}/{filename_save}"
    if Path(XGB_model_save).is_file():
        XGB_classifier_model = pickle.load(open(XGB_model_save, "rb"))
    
        fw_heatmap_eval_df = encoded_ML[feat_labels + target_columns]
        title = f"XGBoost-sample-contaminant-status -- Heat map depicting important features for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used} using featurewiz"
        draw_heatmap(fw_heatmap_eval_df, 
                         XGB_out,
                         dataset_used,
                         metagenome_classifier_used,
                         title)
        
        instance_names_df = true_mask_df_merge[id_cols]
        
        feature_wiz_df = pd.merge(
            instance_names_df, fw_heatmap_eval_df, how="left", left_index=True, right_index=True
        )
        feature_wiz_df = pd.merge(
            feature_wiz_df, other_mask_df, how="left", left_index=True, right_index=True
        )
        feature_wiz_df = feature_wiz_df.drop_duplicates(subset=id_cols)
        feature_wiz_df.reset_index(inplace=True, drop=True)
        
        y = feature_wiz_df["sample_encoded_mask"].values  
        indices = feature_wiz_df.index.values
        feature_wiz_clean = feature_wiz_df.drop(
            all_masks + id_cols, axis=1
        )
        X = feature_wiz_clean.iloc[:, :].values
        
        X_eval = sc.fit_transform(X)
        
        y_pred = XGB_classifier_model.predict(X_eval)
    
        # print("\n\nConformal analysis")
        # conformal_analysis(dataset_used, metagenome_classifier_used, X_eval, y, XGB_classifier_model, [], [], "sample status", XGB_out)
    
        print(f"\nAccuracy score for evaluation sample: {accuracy_score(y, y_pred):.2f}\n")
    
        # Get predicted probabilities for each class
        y_preds_proba = XGB_classifier_model.predict_proba(X_eval)
    
        try:
            with open(f"{XGB_out}/ROC-sample-status-XGB-plotting-data-{metagenome_classifier_used}.json") as json_data:
                plotting_data = json.loads(json_data.read())
            auc_score2, auc_score3, fpr2, tpr2, fpr3, tpr3 = plotting_data
            fpr2, tpr2, fpr3, tpr3 = (
                np.array(fpr2),
                np.array(tpr2),
                np.array(fpr3),
                np.array(tpr3),
            )
    
            # roc curve for models
            fpr1, tpr1, thresh1 = roc_curve(y, y_preds_proba[:, 1], pos_label=1)
    
            random_probs = [0 for i in range(len(y))]
            p_fpr, p_tpr, _ = roc_curve(y, random_probs, pos_label=1)
    
            # auc scores
            auc_score1 = roc_auc_score(y, y_preds_proba[:, 1])
    
            centrifuge_plot = "viral, fungal and bacterial database"
            centrifuge_label = "viral-fungal-bacterial"
            
            plt.style.use("seaborn")
            f, ax = plt.subplots(figsize=(15, 15))
            # plot roc curves
            plt.plot(
                fpr1,
                tpr1,
                linestyle="--",
                color="red",
                label=f"XGB-eval",
            )
            plt.plot(fpr2, tpr2, linestyle="--", color="green", label="XGB-test")
            plt.plot(fpr3, tpr3, linestyle="--", color="orange", label="XGB-train")
            plt.plot(p_fpr, p_tpr, linestyle="--", color="blue")
            ax.tick_params(axis="x", labelsize=25)
            ax.tick_params(axis="y", labelsize=25)
            # title
            title = f"AUC-ROC curve for sample contaminant status for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used} using the {centrifuge_plot}"
            new_title = clean_strings(title, 120)
            plt.title(new_title, size=40)
            # x label
            plt.xlabel("False Positive Rate", size=40)
            # y label
            plt.ylabel("True Positive rate", size=40)
    
            plt.text(
                0.8,
                0.45,
                f"AUC-eval: {auc_score1:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                size=30,
            )
            plt.text(
                0.8,
                0.4,
                f"AUC-test: {auc_score2:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                size=30,
            )
            plt.text(
                0.8,
                0.35,
                f"AUC-train: {auc_score3:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                size=30,
            )
    
            plt.legend(loc=4, prop={"size": 30}, markerscale=10)
            plt.savefig(f"{XGB_out}/{title}.png", dpi=300)
            plt.show()
    
        except:
            print("Missing examples of class X, cannot plot")
    

        evaluate_model_performance(y_pred, XGB_classifier_model, 
                                           X_eval, _, _, 
                                           y, metagenome_classifier_used,
                                           "sample contaminant",
                                           "evaluate",
                                           encoding_dict)
    
        sample_status = true_mask_df_merge[
            ["sample", "name", "mask", "sample_true_mask"]
        ]
        sample_status.loc[indices, "XGB_sample_prediction"] = XGB_classifier_model.predict(
            X_eval
        )
        
        sample_status.loc[indices, "sample-status"] = "eval"
    
        encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
        sample_status = sample_status.replace({"XGB_sample_prediction": encoding_dict_rev})
        
        sample_status["sample-pred-correct"] = sample_status.apply(check_prediction, axis=1)
        
        return sample_status
    

ce_eval_sample_status = run_xgboost_sample_status_evaluate("centrifuge", "evaluation", true_mask, evaluation_df_ce, cd_data_cols, ce_id_cols, vees, XGB_out, save, run_gridsearch)
hs_eval_sample_status = run_xgboost_sample_status_evaluate("BLAST", "evaluation", true_mask, evaluation_df_bn, hd_data_cols, bl_id_cols, vees, XGB_out, save, run_gridsearch)


def run_xgboost_genuine_contaminant_status_evaluate(metagenome_classifier_used: str,
                              dataset_used: str,
                              true_mask: dict,
                              input_df: pd.DataFrame,
                              input_df_xgboost: pd.DataFrame,
                              data_cols: list, 
                              id_cols: list, 
                              vees: str,
                              true_mask_df: pd.DataFrame, 
                              XGB_out: str,
                              save: bool,
                              run_gridsearch: bool) -> pd.DataFrame:

    filename_save = f"xgboost_genuine_contaminant_status-{metagenome_classifier_used}.sav"
    
    genuine_contaminant_df = pd.merge(input_df, input_df_xgboost[['sample_true_mask', 'XGB_sample_prediction', 'sample-status',
           'sample-pred-correct']], left_index=True, right_index=True)    
    
    # convert prediction T/F into encoded 0,1
    encoding_dict_samp = {"True_positive": 1, "True_negative": 0}
    target_columns_samp = ["sample-pred-TP", "sample-pred-TN"]
    label_partial_func = partial(custom_samp_label_encoder, encoding_dict_samp)
    genuine_contaminant_df["XGB_sample_prediction"] = genuine_contaminant_df.apply(label_partial_func, axis=1)
    
    enc_samp_partial_func = partial(one_hot_samp_encoder, encoding_dict_samp)
    genuine_contaminant_df["temp-mask"] = genuine_contaminant_df.apply(enc_samp_partial_func, axis=1)    
    
    genuine_contaminant_df[target_columns_samp] = pd.DataFrame(
        genuine_contaminant_df["temp-mask"].tolist(), index=genuine_contaminant_df.index
    )
    
    genuine_contaminant_df = genuine_contaminant_df.drop(["temp-mask"], axis=1)
    
    encoded_ML = genuine_contaminant_df[data_cols + ["mask"] + target_columns_samp]
    
    encoding_dict = {"False_positive": 0, "True_positive": 1}
    target_columns = ["mask-FP", "mask-TP"]
    
    label_partial_func = partial(custom_contaminant_label_encoder, encoding_dict)
    encoded_ML["encoded_mask"] = encoded_ML.apply(label_partial_func, axis=1)
    
    enc_partial_func = partial(one_hot_contaminant_encoder, encoding_dict)
    encoded_ML["temp-mask"] = encoded_ML.apply(enc_partial_func, axis=1)
    encoded_ML[target_columns] = pd.DataFrame(
        encoded_ML["temp-mask"].tolist(), index=encoded_ML.index
    )
    encoded_ML = encoded_ML.drop(["temp-mask"], axis=1)
    
    all_masks = ["mask", "encoded_mask", "mask-FP", "mask-TP"]
    
    mask_to_drop = list(set(all_masks) - set(target_columns))
    other_mask_df = encoded_ML[mask_to_drop]
    encoded_ML = encoded_ML.drop(mask_to_drop, axis=1)
    
    with open(
        f"{XGB_out}/feature-wiz-xgboost-features-{vees}-genuine-contaminant-{metagenome_classifier_used}.json"
    ) as json_sample_XGB:
        feat_labels = json.loads(json_sample_XGB.read())
        
    sc = pickle.load(open(f"{XGB_out}/sc-genuine-contaminant-XGB-scaler-{metagenome_classifier_used}.pkl", "rb"))
    
    filename_save = f"xgboost_genuine_contaminant_status-{metagenome_classifier_used}.sav"
    XGB_model_save = f"{XGB_out}/{filename_save}"
    if Path(XGB_model_save).is_file():
        XGB_classifier_model = pickle.load(open(XGB_model_save, "rb"))
    
        fw_heatmap_eval_df = encoded_ML[feat_labels + target_columns]
    
        title = f"XGBoost-genuine-contaminant-status -- Heat map depicting important features for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used} using featurewiz"
        draw_heatmap(fw_heatmap_eval_df, 
                         XGB_out,
                         dataset_used,
                         metagenome_classifier_used,
                         title)
    
        instance_names_df = genuine_contaminant_df[id_cols]
        
        feature_wiz_df = pd.merge(
            instance_names_df, fw_heatmap_eval_df, how="left", left_index=True, right_index=True
        )
        feature_wiz_df = pd.merge(
            feature_wiz_df, other_mask_df, how="left", left_index=True, right_index=True
        )
        feature_wiz_df = feature_wiz_df.drop_duplicates(subset=id_cols)
        feature_wiz_df.reset_index(inplace=True, drop=True)
        
        y = feature_wiz_df["encoded_mask"].values  
        indices = feature_wiz_df.index.values
        feature_wiz_clean = feature_wiz_df.drop(
            all_masks + id_cols, axis=1
        )
        X = feature_wiz_clean.iloc[:, :].values
        
        X_eval = sc.fit_transform(X)
        
        y_pred = XGB_classifier_model.predict(X_eval)

        # print("\n\nConformal analysis")
        # conformal_analysis(dataset_used, metagenome_classifier_used, X_eval, y, XGB_classifier_model, [], [], "genuine status", XGB_out)
    
        print(f"\nAccuracy score for evaluation sample: {accuracy_score(y, y_pred):.2f}\n")
    
        # Get predicted probabilities for each class
        y_preds_proba = XGB_classifier_model.predict_proba(X_eval)
        
        try:
            with open(f"{XGB_out}/ROC-adventitious-agent-XGB-plotting-data-{metagenome_classifier_used}.json") as json_data:
                plotting_data = json.loads(json_data.read())
            auc_score2, auc_score3, fpr2, tpr2, fpr3, tpr3 = plotting_data
            fpr2, tpr2, fpr3, tpr3 = (
                np.array(fpr2),
                np.array(tpr2),
                np.array(fpr3),
                np.array(tpr3),
            )
    
            # roc curve for models
            fpr1, tpr1, thresh1 = roc_curve(y, y_preds_proba[:, 1], pos_label=1)
    
            random_probs = [0 for i in range(len(y))]
            p_fpr, p_tpr, _ = roc_curve(y, random_probs, pos_label=1)
    
            # auc scores
            auc_score1 = roc_auc_score(y, y_preds_proba[:, 1])
    
            centrifuge_plot = "viral, fungal and bacterial database"
            # centrifuge_label = "viral-fungal-bacterial"
            
            plt.style.use("seaborn")
            f, ax = plt.subplots(figsize=(15, 15))
            # plot roc curves
            plt.plot(
                fpr1,
                tpr1,
                linestyle="--",
                color="red",
                label=f"XGB-eval",
            )
            plt.plot(fpr2, tpr2, linestyle="--", color="green", label="XGB-test")
            plt.plot(fpr3, tpr3, linestyle="--", color="orange", label="XGB-train")
            plt.plot(p_fpr, p_tpr, linestyle="--", color="blue")
            ax.tick_params(axis="x", labelsize=25)
            ax.tick_params(axis="y", labelsize=25)
            # title
            title = f"AUC-ROC curve for genuine contaminant status for dataset - {dataset_used} and metagenome classifier - {metagenome_classifier_used} using the {centrifuge_plot}"
            new_title = clean_strings(title, 120)
            plt.title(new_title, size=40)
            # x label
            plt.xlabel("False Positive Rate", size=40)
            # y label
            plt.ylabel("True Positive rate", size=40)
    
            plt.text(
                0.8,
                0.45,
                f"AUC-eval: {auc_score1:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                size=30,
            )
            plt.text(
                0.8,
                0.4,
                f"AUC-test: {auc_score2:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                size=30,
            )
            plt.text(
                0.8,
                0.35,
                f"AUC-train: {auc_score3:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                size=30,
            )
    
            plt.legend(loc=4, prop={"size": 30}, markerscale=10)
            plt.savefig(f"{XGB_out}/{title}.png", dpi=300)
            plt.show()
    
        except:
            print("Missing examples of class X, cannot plot")
        
        evaluate_model_performance(y_pred, XGB_classifier_model, 
                                           X_eval, _, _, 
                                           y, metagenome_classifier_used,
                                           "genuine contaminant",
                                           "evaluate",
                                           encoding_dict)
    
        gen_con_status = genuine_contaminant_df[
            ["sample", "name", "mask", "sample_true_mask", "XGB_sample_prediction"]
        ]
        gen_con_status.loc[indices, "XGB_gen_con_prediction"] = XGB_classifier_model.predict(
            X_eval
        )
    
        gen_con_status.loc[indices, "sample-status"] = "eval"
    
        encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
        gen_con_status = gen_con_status.replace({"XGB_gen_con_prediction": encoding_dict_rev})
        
        gen_con_status["gen-con-pred-correct"] = gen_con_status.apply(check_mask_correct_pred, axis=1)
        
        return gen_con_status

        
ce_eval_gen_con_status = run_xgboost_genuine_contaminant_status_evaluate("centrifuge", "evaluation", true_mask, evaluation_df_ce, ce_eval_sample_status, cd_data_cols, ce_id_cols, vees, true_mask_df, XGB_out, save, run_gridsearch)
hs_eval_gen_con_status = run_xgboost_genuine_contaminant_status_evaluate("BLAST", "evaluation", true_mask, evaluation_df_bn, hs_eval_sample_status, hd_data_cols, bl_id_cols, vees, true_mask_df, XGB_out, save, run_gridsearch)

################## EVALUATION ##################
################## MACHINE LEARNING ##################
#
#
#
################## ASSESSMENT OF EFFICACY ##################

def calculate_accuracy(df: pd.DataFrame, col: str, 
                       _type_: str, mg: str) -> float:
    s = df[col].value_counts()
    correct_predictions = s.iloc[0]
    accuracy = (correct_predictions / len(df[col])) 
    # accuracy_percent = accuracy * 100
    # print(f"Data type: {_type_}; Metagenomic classifier: {mg}; Accuracy: {accuracy_percent:.1f}%")
    # print(f"{_type_}-{mg}; Accuracy: {accuracy_percent:.1f}%")
    return accuracy
    
cs_acc = calculate_accuracy(ce_eval_sample_status, "sample-pred-correct", "evaluation-sample", "centrifuge")
hs_acc = calculate_accuracy(hs_eval_sample_status, "sample-pred-correct", "evaluation-sample", "BLAST")
cg_acc = calculate_accuracy(ce_eval_gen_con_status, "gen-con-pred-correct", "evaluation-genuine", "centrifuge")
hg_acc = calculate_accuracy(hs_eval_gen_con_status, "gen-con-pred-correct", "evaluation-genuine", "BLAST")

import math    
# from statsmodels.stats.proportion import proportion_confint

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


calculate_classification_accuracy_interval(cs_acc, ce_eval_sample_status, "evaluation-sample", "centrifuge")
calculate_classification_accuracy_interval(hs_acc, hs_eval_sample_status, "evaluation-sample", "BLAST")
calculate_classification_accuracy_interval(cg_acc, ce_eval_gen_con_status, "evaluation-genuine", "centrifuge")
calculate_classification_accuracy_interval(hg_acc, hs_eval_gen_con_status, "evaluation-genuine", "BLAST")


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
        decision_pivot["True_negative"] = np.nan()
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


def calculate_twin_model_accuracy(df: pd.DataFrame, mg: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    encoding_dict = {"True_negative": 0, "True_positive": 1}
    encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
    eval_gen_con_status = df.replace({"XGB_sample_prediction": encoding_dict_rev})
            
    eval_gen_con_status["Decision"] = eval_gen_con_status.apply(make_decision, axis=1)
    eval_gen_con_status["Decision_mask"] = eval_gen_con_status.apply(decision_mask, axis=1)
    eval_gen_con_status["Dec_Correct"] = eval_gen_con_status.apply(check_final_decision, axis=1)
    
    eval_gen_con_status["check-tts-status-decision"] = eval_gen_con_status.apply(check_tts_status, axis=1)
    
    cg_acc_dec = calculate_accuracy(eval_gen_con_status, "check-tts-status-decision", "evaluation-sample-genuine", mg)
    calculate_classification_accuracy_interval(cg_acc_dec, eval_gen_con_status, "evaluation-sample-genuine", mg)
    
    decision_pivot = generate_decision_df(eval_gen_con_status)
    
    decision_pivot["cfu"] = decision_pivot["sample"].str.split("_", expand = True)[3]
    
    return decision_pivot, eval_gen_con_status


decision_pivot_ce, eval_gen_con_status_ce = calculate_twin_model_accuracy(ce_eval_gen_con_status, "centrifuge")
decision_pivot_hs, eval_gen_con_status_hs = calculate_twin_model_accuracy(hs_eval_gen_con_status, "BLAST")  
