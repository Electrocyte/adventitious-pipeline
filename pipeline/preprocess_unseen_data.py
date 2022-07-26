#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:02:49 2022

@author: mangi
"""

from pathlib import Path
import pandas as pd
import numpy as np
from functools import partial
import pickle
import json
from typing import Tuple

import seaborn as sns
import matplotlib.pyplot as plt
from string import whitespace


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


def draw_heatmap(fw_heatmap_df: pd.DataFrame, 
                  XGB_out: str,
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


from sklearn.model_selection import cross_val_score
import statistics
from sklearn.metrics import classification_report, confusion_matrix
# from sklearn.metrics import accuracy_score, recall_score, precision_score


def evaluate_model_performance(y_pred, 
                               XGB_classifier_model, 
                               X_eval, 
                               y_test, 
                               metagenome_classifier_used: str, 
                               ml_model: str,
                               _type_: str, 
                               encoding_dict: dict):    
    # 5 folds, scored on accuracy
    # https://scikit-learn.org/stable/modules/model_evaluation.html#scoring-parameter

    cvs2 = cross_val_score(
        XGB_classifier_model, X_eval, y_test, cv=10, scoring="accuracy"
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
    
    print(f"Cross-val #2 mean: {statistics.mean(cvs2):.2f}")
    print(f"Cross-val #2 stdev: {statistics.stdev(cvs2):.2f}")
  
    
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
        

# def check_TPs(df: pd.DataFrame, species: dict):
#     df["strain"] = df["sample"].str.split("_", expand=True)[2]
#     df.reset_index(inplace=True,drop=True)
#     true_masks = []
#     for k, vv in species.items():
#         if isinstance(vv, list):
#             for v in vv:
#                 true_mask = df.loc[
#                     (df["strain"].str.contains(v)) & (df.name.str.contains(k))
#                 ].index
#                 df.loc[true_mask, "mask"] = "True_positive"
#                 true_masks.append(true_mask)
#         else:
#             true_mask = df.loc[
#                 (df["strain"].str.contains(vv)) & (df.name.str.contains(k))
#             ].index
#             df.loc[true_mask, "mask"] = "True_positive"
#             true_masks.append(true_mask)
#     true_masks = list(set([item for sublist in true_masks for item in sublist]))
#     TPs = df.iloc[true_masks]
#     # df.loc[TPs, "sample-mask"] = True
#     TP_samples = list(TPs["sample"].unique())
#     print(TP_samples)
#     def add_TP(row):

#         if row["sample"] in TP_samples:
#             row["sample-mask"] = True
#             return row
#         else:
#             return row
    
#     df = df.apply(add_TP, axis = 1)
#     return df


def prepare_data_sample(true_mask: dict, 
                        input_df: pd.DataFrame,
                        ) -> (pd.DataFrame, list, 
                              pd.DataFrame, list, 
                              pd.DataFrame, dict,
                              pd.DataFrame):

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
        left_on=["sample"],
        right_on=["sample"],
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
    # target_columns = ["sample-mask-TP", "sample-mask-TN"]
    
    label_partial_func = partial(custom_label_encoder, encoding_dict)
    true_mask_df_merge["sample_encoded_mask"] = true_mask_df_merge.apply(label_partial_func, axis=1)
    
    return true_mask_df_merge


def debug_true_labels(XGB_out: str,
                      df: pd.DataFrame):
    json_mask = "2022-07.json"
    json_file = f"{XGB_out}{json_mask}"
    
    with open(json_file) as json_data:
        true_mask = json.loads(json_data.read())
        
    return prepare_data_sample(true_mask, df)


def count_predictions(sample_status: pd.DataFrame, col: str) -> pd.DataFrame:    
    predictions_per_sample = sample_status.groupby(["sample"])[col].value_counts()
    predictions_per_sample = pd.DataFrame(predictions_per_sample)    
    predictions_per_sample = predictions_per_sample.rename(columns={col: "counts"})
    
    predictions_per_sample.reset_index(inplace=True)
    predictions_per_sample.set_index('sample', inplace=True)
    pps_pivot = pd.pivot_table(predictions_per_sample, values='counts', index=['sample'], columns = [col])
    pps_pivot = pps_pivot.fillna(0)
    return pps_pivot


def run_unseen_sample_analysis(metagenome_classifier_used: str,
                               vees: str,
                               XGB_out: str, 
                               id_cols: list,
                               debug: bool,
                               unseen_df: pd.DataFrame) -> Tuple[pd.DataFrame, 
                                                                 pd.DataFrame]:

    with open(
        f"{XGB_out}/feature-wiz-xgboost-features-{vees}-{metagenome_classifier_used}-sample-status.json"
    ) as json_sample_XGB:
        feat_labels = json.loads(json_sample_XGB.read())
    
    sc = pickle.load(open(f"{XGB_out}/sc-sample-status-XGB-scaler-{metagenome_classifier_used}.pkl", "rb"))
    
    filename_save = f"xgboost_sample_status-{metagenome_classifier_used}.sav"
    XGB_model_save = f"{XGB_out}/{filename_save}"
    
    if Path(XGB_model_save).is_file():
        XGB_classifier_model = pickle.load(open(XGB_model_save, "rb"))
    
        unseen_fw_df = unseen_df[feat_labels]
        
        if debug:
            title = f"XGBoost-sample-contaminant-status -- unseen data - {metagenome_classifier_used} using featurewiz"
            draw_heatmap(unseen_fw_df, 
                              XGB_out,
                              metagenome_classifier_used,
                              title)
    
        instance_names_df = unseen_df[id_cols]
    
        feature_wiz_df = pd.merge(
            instance_names_df, unseen_fw_df, how="left", left_index=True, right_index=True
        )
        
        feature_wiz_df = feature_wiz_df.drop_duplicates(subset=id_cols)
        feature_wiz_df.reset_index(inplace=True, drop=True)
    
        indices = feature_wiz_df.index.values
        uns_clean = feature_wiz_df.drop(
            id_cols, axis=1
        )
        X = uns_clean.iloc[:, :].values
    
        X_eval = sc.fit_transform(X)
        
        y_pred = XGB_classifier_model.predict(X_eval)
        
        # Get predicted probabilities for each class
        # y_preds_proba = XGB_classifier_model.predict_proba(X_eval)   
    
        if debug:
    
            encoding_dict = {"True_negative": 0, "True_positive": 1}
            
            y_mask = debug_true_labels(XGB_out, feature_wiz_df)
            
            # json_mask = "2022-07.json"
            # json_file = f"{XGB_out}{json_mask}"
        
            # with open(json_file) as json_data:
            #     true_mask = json.loads(json_data.read())
        
            # y_mask = prepare_data_sample(true_mask, feature_wiz_df)
                
            y = y_mask["sample_encoded_mask"].values
            
            # handle only positives; need check before running this
            print(y_mask["sample_encoded_mask"].unique())
            if 0 not in y_mask["sample_encoded_mask"].unique():
                print("Adding dummy vars")
                newrow = np.zeros(X_eval.shape[1])
                X_evalA = X_eval
                for i in range(3):
                    X_evalA = np.vstack([X_evalA, newrow])
                y_evalA = np.pad(y, (0, 3), 'constant')
                y_predA = np.pad(y_pred, (0, 3), 'constant')
                evaluate_model_performance(y_predA, XGB_classifier_model, 
                                                   X_evalA, 
                                                   y_evalA, metagenome_classifier_used,
                                                   "sample contaminant",
                                                   "unknown",
                                                   encoding_dict)  
            else:
                "Skipping addition of dummy vars"
                evaluate_model_performance(y_pred, XGB_classifier_model, 
                                                   X_eval, 
                                                   y, metagenome_classifier_used,
                                                   "sample contaminant",
                                                   "unknown",
                                                   encoding_dict)   
        sample_status = feature_wiz_df[
            ["sample", "name"]
        ]
        sample_status.loc[indices, "XGB_sample_prediction"] = XGB_classifier_model.predict(
            X_eval
        )
        
        sample_status.loc[indices, "sample-status"] = "unknowns"
    
        encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
        sample_status = sample_status.replace({"XGB_sample_prediction": encoding_dict_rev})
    
        samp_pivot = count_predictions(sample_status, "XGB_sample_prediction")
        return samp_pivot, sample_status
    

def custom_samp_label_encoder(encoding_dict: dict, row: int):
    for name, encoded in encoding_dict.items():
        if name in row["XGB_sample_prediction"]:
            return encoded


def one_hot_samp_encoder(encoding_dict: dict, row: int):
    for name, encoded in encoding_dict.items():
        if encoded == row["XGB_sample_prediction"]:
            if encoded == 1:
                return (1, 0)
            if encoded == 0:
                return (0, 1)
            else:
                return (0, 0)
         

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

         
def run_unseen_genuine_analysis(metagenome_classifier_used: str,
                               vees: str,
                               XGB_out: str, 
                               id_cols: list,
                               debug: bool,
                               unseen_df: pd.DataFrame,
                               sample_status: pd.DataFrame) -> Tuple[pd.DataFrame, 
                                                                 pd.DataFrame]:
    with open(
        f"{XGB_out}/feature-wiz-xgboost-features-{vees}-genuine-contaminant-{metagenome_classifier_used}.json"
    ) as json_sample_XGB:
        feat_labels = json.loads(json_sample_XGB.read()) 
        
    sc = pickle.load(open(f"{XGB_out}/sc-genuine-contaminant-XGB-scaler-{metagenome_classifier_used}.pkl", "rb"))
    
    filename_save = f"xgboost_genuine_contaminant_status-{metagenome_classifier_used}.sav"
    XGB_model_save = f"{XGB_out}/{filename_save}"

    encoding_dict_samp = {"True_positive": 1, "True_negative": 0}
    target_columns_samp = ["sample-pred-TP", "sample-pred-TN"]

    label_partial_func = partial(custom_samp_label_encoder, encoding_dict_samp)
    sample_status["XGB_sample_prediction"] = sample_status.apply(label_partial_func, axis=1)
    
    enc_samp_partial_func = partial(one_hot_samp_encoder, encoding_dict_samp)
    sample_status["temp-mask"] = sample_status.apply(enc_samp_partial_func, axis=1)   

    sample_status[target_columns_samp] = pd.DataFrame(
        sample_status["temp-mask"].tolist(), index=sample_status.index
    )    

    sample_status = sample_status.drop(["temp-mask"], axis=1)
    
    unseen_df_clean = unseen_df.drop(["sample", "name"], axis=1)
    feature_wiz_df = pd.merge(
        sample_status, unseen_df_clean, how="left", left_index=True, right_index=True
    )
    
    if Path(XGB_model_save).is_file():
        XGB_classifier_model = pickle.load(open(XGB_model_save, "rb"))
    
        unseen_fw_df = feature_wiz_df[feat_labels]
        
        if debug:
            # plot heat map
            title = f"XGBoost-genuine-contaminant-status -- unseen data - {metagenome_classifier_used} using featurewiz"
            draw_heatmap(unseen_fw_df, 
                              XGB_out,
                              metagenome_classifier_used,
                              title)
    
        instance_names_df = unseen_df[id_cols]
    
        feature_wiz_df = pd.merge(
            instance_names_df, unseen_fw_df, how="left", left_index=True, right_index=True
        )
        
        feature_wiz_df = feature_wiz_df.drop_duplicates(subset=id_cols)
        feature_wiz_df.reset_index(inplace=True, drop=True)
    
        indices = feature_wiz_df.index.values
        uns_clean = feature_wiz_df.drop(
            id_cols, axis=1
        )
        X = uns_clean.iloc[:, :].values
    
        X_eval = sc.fit_transform(X)
        
        y_pred = XGB_classifier_model.predict(X_eval)
        
        # Get predicted probabilities for each class
        y_preds_proba = XGB_classifier_model.predict_proba(X_eval)   
        encoding_dict = {"False_positive": 0, "True_positive": 1}
        
        if debug:
            # investigate samples ONLY IF THE MASK IS KNOWN
            species = {"Pseudomonas aeruginosa":["PA", "9027"], "Cutibacterium acnes":["Cacnes","Pacnes"], \
                        "Escherichia coli":["EC"], "Klebsiella pneumoniae":["Klebpneu"], \
                          "Candida":["Calbicans"], "Staphylococcus aureus":["Saureus"], \
                            "Bacillus subtilis": ["Bsubtilis"]}
            masked_df = apply_mask(feature_wiz_df, species)
            
            def custom_contaminant_label_encoder(encoding_dict: dict, row: int):
                for name, encoded in encoding_dict.items():
                    if name in row["mask"]:
                        return encoded


            def one_hot_contaminant_encoder(encoding_dict: dict, row: int):
                for name, encoded in encoding_dict.items():
                    if encoded == row["encoded_mask"]:
                        if encoded == 1:
                            return (1, 0)
                        if encoded == 0:
                            return (0, 1)
                        else:
                            return (0, 0)

            target_columns = ["mask-FP", "mask-TP"]
            
            label_partial_func = partial(custom_contaminant_label_encoder, encoding_dict)
            masked_df["encoded_mask"] = masked_df.apply(label_partial_func, axis=1)
            
            enc_partial_func = partial(one_hot_contaminant_encoder, encoding_dict)
            masked_df["temp-mask"] = masked_df.apply(enc_partial_func, axis=1)
            masked_df[target_columns] = pd.DataFrame(
                masked_df["temp-mask"].tolist(), index=masked_df.index
            )
            masked_df = masked_df.drop(["temp-mask"], axis=1) 
               
            y = masked_df["encoded_mask"].values
            
            # handle only positives; need check before running this
            print(masked_df["encoded_mask"].unique())
            if 0 not in masked_df["encoded_mask"].unique():
                print("Adding dummy vars")
                newrow = np.zeros(X_eval.shape[1])
                X_evalA = X_eval
                for i in range(3):
                    X_evalA = np.vstack([X_evalA, newrow])
                y_evalA = np.pad(y, (0, 3), 'constant')
                y_predA = np.pad(y_pred, (0, 3), 'constant')
                evaluate_model_performance(y_predA, XGB_classifier_model, 
                                                    X_evalA, 
                                                    y_evalA, metagenome_classifier_used,
                                                    "genuine contaminant",
                                                    "unknown",
                                                    encoding_dict)  
            else:
                "Skipping addition of dummy vars"
                evaluate_model_performance(y_pred, XGB_classifier_model, 
                                                    X_eval, 
                                                    y, metagenome_classifier_used,
                                                    "genuine contaminant",
                                                    "unknown",
                                                    encoding_dict)   
        cont_genu_status = feature_wiz_df[
            ["sample", "name"]
        ]
        cont_genu_status.loc[indices, "XGB_genuine_prediction"] = XGB_classifier_model.predict(
            X_eval
        )
        
        cont_genu_status.loc[indices, "sample-status"] = "unknowns"
    
        encoding_dict_rev = {value: key for (key, value) in encoding_dict.items()}
        cont_genu_status = cont_genu_status.replace({"XGB_genuine_prediction": encoding_dict_rev})
    
        samp_pivot = count_predictions(cont_genu_status, "XGB_genuine_prediction")
        return samp_pivot, cont_genu_status