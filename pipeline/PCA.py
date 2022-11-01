#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 14:46:27 2022

@author: mangi
"""

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.decomposition import PCA
from typing import Tuple


def relabel(df):
    df.loc[df["sample"].str.contains("Calbicans"), "short"] = "C. albicans"
    df.loc[df["sample"].str.contains("Saureus"), "short"] = "S. aureus"
    df.loc[df["sample"].str.contains("Cacnes"), "short"] = "C. acnes"
    df.loc[df["sample"].str.contains("clost"), "short"] = "C. sporogenes"
    df.loc[df["sample"].str.contains("PAO1"), "short"] = "P. a PAO1"
    df.loc[df["sample"].str.contains("9027"), "short"] = "P. a 9027"
    df.loc[df["sample"].str.contains("Klebpneu"), "short"] = "K. pneumoniae"
    df.loc[df["sample"].str.contains("Tcells"), "short"] = "Negative"
    df.loc[df["sample"].str.contains("Cell-free-medium"), "short"] = "Negative"
    df.loc[df["sample"].str.contains("Plain"), "short"] = "Negative"
    df.loc[df["sample"].str.contains("Medium"), "short"] = "Negative"
    return df


def troubleshoot(X_all: np.ndarray, df: pd.DataFrame):
    ts_df = pd.DataFrame(X_all)
    ts_cat = pd.concat([ts_df, df["sample"]], axis=1)
    ts_cat = relabel(ts_cat)
    def get_cmap(n, name='hsv'):
        '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
        RGB color; the keyword argument name must be a standard mpl colormap name.'''
        return plt.cm.get_cmap(name, n)
        
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    
    cmap = get_cmap(len(ts_cat["short"].unique()))
    for n, i in enumerate(ts_cat["short"].unique()):
        sample = ts_cat.loc[ts_cat["short"] == i]
        xs = sample[0]
        ys = sample[1]
        zs = sample[2]
        ax.scatter(xs, ys, zs, c=cmap(n), marker='^', alpha = 0.5, label = i)

        if n+1 == len(ts_cat["short"].unique()):
            ax.scatter(xs, ys, zs, c="black", marker='o', alpha = 1, s = 50, label = i)
    ax.legend()
    plt.show()


def process_data(df: pd.DataFrame, 
                 feat_labels: list,
                 sc,
                 metagenome_classifier_used: str, 
                 ctrl_idxs: list) -> np.ndarray:
    clean_df = df[feat_labels]
    
    X = clean_df.iloc[:, :].values
    
    # Applying PCA
    pca = PCA(n_components = 3)

    X = sc.fit_transform(X)
    X_all = pca.fit_transform(X)    
    
    return X_all


def get_xyz(X_all: np.ndarray, ctrl_idxs: list) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:    
    xs = X_all[:,0]
    ys = X_all[:,1]
    zs = X_all[:,2]
    
    try:
        TNs = [X_all[i-1] for i in ctrl_idxs]
        TNs = np.vstack(TNs)
        TNxs = TNs[:,0]
        TNys = TNs[:,1]
        TNzs = TNs[:,2]
        return TNxs, TNys, TNzs, xs, ys, zs
    
    except:
        return 0, 0, 0, xs, ys, zs


def run(df: pd.DataFrame, 
        metagenome_classifier_used: str, 
        unseen: bool, data_aug: bool,
        XGB_out: str):

    # XGB_out = 'D:/SequencingData/Harmonisation/DNA/analysis//ML_training-VFB//OCS-XGBoost-VFB/'
    
    filename = f"{XGB_out}/data-aug-{metagenome_classifier_used}-sample-status.csv"
    
    if data_aug:
        df = pd.read_csv(filename)
    
    vees = 'PA-Cacnes-Pacnes-EC-Klebpneu-Calbicans-Saureus-Bsubtilis'
    
    with open(
        f"{XGB_out}/feature-wiz-xgboost-features-{vees}-{metagenome_classifier_used}-sample-status.json"
    ) as json_sample_XGB:
        feat_labels = json.loads(json_sample_XGB.read())
    
    controls = ["Cell-free-medium", "Tcells", "Plain-medium", \
                "CellfreeMedia", "PlainMedia", "TCell", \
                    "TC-G-REX-CD19", "TC-Well-Plate-CD19"]

    df.reset_index(inplace=True, drop=True)
    ctrl_idxs = []
    for ctrl in controls:
        ctrl_idx = df.loc[df["sample"].str.contains(ctrl)].index
        ctrl_idxs.extend(ctrl_idx)
    
    df.loc[df.index.isin(ctrl_idxs), "label"] = "0" # TN
    df.loc[~df.index.isin(ctrl_idxs), "label"] = "1" # TP  

    if not unseen:
        min_class = df["label"].value_counts().min()    
        tn_sample = df.loc[df["label"] == "0"].sample(min_class)
        tp_sample = df.loc[df["label"] == "1"].sample(min_class)
        tn_sample.to_csv(f"{XGB_out}/TN-subsample-{metagenome_classifier_used}-PCA.csv", index=False)
        tp_sample.to_csv(f"{XGB_out}/TP-subsample-{metagenome_classifier_used}-PCA.csv", index=False)
    else:
        TN_subset = pd.read_csv(f"{XGB_out}/TN-subsample-{metagenome_classifier_used}-PCA.csv")
        TP_subset = pd.read_csv(f"{XGB_out}/TP-subsample-{metagenome_classifier_used}-PCA.csv")
        standard_df = pd.concat([TN_subset, TP_subset])
        standard_df.reset_index(inplace=True, drop=True)
        std_idx = standard_df.loc[standard_df["label"] == 0].index

    sc = pickle.load(open(f"{XGB_out}/sc-sample-status-XGB-scaler-{metagenome_classifier_used}.pkl", "rb"))  
    
    X_all = process_data(df, feat_labels, sc, metagenome_classifier_used, ctrl_idxs)
    X_all_standard = process_data(standard_df, feat_labels, sc, metagenome_classifier_used, std_idx)

    TNxs, TNys, TNzs, xs, ys, zs = get_xyz(X_all, ctrl_idxs)
    std_TNxs, std_TNys, std_TNzs, std_xs, std_ys, std_zs = get_xyz(X_all_standard, std_idx)

    troubleshoot(X_all, df)  

    if unseen:
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection='3d')
        
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        
        ax.scatter(xs, ys, zs, c="blue", marker='^', alpha = 0.5)
        ax.scatter(std_xs, std_ys, std_zs, c="grey", marker='^', alpha = 0.1)
        
        try:
            ax.scatter(TNxs, TNys, TNzs, c="red", marker='o', s = 50, alpha = 0.5)
        except: 
            pass
        
        ax.scatter(std_TNxs, std_TNys, std_TNzs, c="black", marker='o', s = 50, alpha = 0.1)
        
        plt.show()
 
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection='3d')
        
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        
        ax.scatter(xs, ys, zs, c="blue", marker='^', alpha = 0.5)

        try:
            ax.scatter(TNxs, TNys, TNzs, c="red", marker='o', s = 50, alpha = 0.5)
        except: 
            pass
        
        plt.show()   
 
    if not unseen:
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(projection='3d')
        
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        
        ax.scatter(xs, ys, zs, c="blue", marker='^', alpha = 0.5)
        ax.scatter(TNxs, TNys, TNzs, c="red", marker='o', s = 50, alpha = 0.5)
        
        plt.show()