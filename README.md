# Adventitious Agent Analysis Pipeline(s)
The pipeline is composed of three tools: metagenomic classification pipeline, read quality analysis and machine learning classification pipeline.
<br /><br />

# Table of Contents  
[Metagenomic classification pipeline](https://github.com/Electrocyte/adventitious-pipeline/blob/main/README.md#metagenomic-classification-pipeline)  
[Read quality analysis](https://github.com/Electrocyte/adventitious-pipeline/blob/main/README.md#read-quality-analysis)  
[Machine learning pipeline](https://github.com/Electrocyte/adventitious-pipeline/blob/main/README.md#machine-learning-pipeline)  
[Pipeline inputs](https://github.com/Electrocyte/adventitious-pipeline/blob/main/README.md#pipeline-inputs)
<br /><br />

## Metagenomic classification pipeline
* Reads generated by the Oxford Nanopore are pre-processed using the [Guppy](https://github.com/nanoporetech) basecaller, 
* Demultiplexed ([qcat](https://github.com/nanoporetech/qcat)) and adapters removed ([porechop](https://github.com/rrwick/Porechop)). 
* Optionally if a host genome is provided these host reads can be removed from the total reads to improve throughput runtime and reduce noise from the host reads.
* Metagenomic classifiication is completed using [centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml), [high-speed BLASTN](https://github.com/chenying2016/queries) and [krakenuniq](https://github.com/fbreitwieser/krakenuniq). Each classifier is run in parallel with the sample, runtime will be dependent on processing power.
* Databases used for metagenomic classification are derived from [NCBI refseq](https://www.ncbi.nlm.nih.gov/refseq/) and the [Reference Viral Database](https://rvdb.dbi.udel.edu/).
* Statistics for each predicted species are generated using data from centrifuge (troubleshooting file) and HS-BLASTn, which have the [pandas `describe()`](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.describe.html) function applied to them.
* Summary statistics for read quality are added to each metagenomic output using readIDs from sequencing summary text file.
<br /><br />

Run Pipeline -  Example run commands for metagenomic pipeline for bacterial, fungal and viral databases
```
/adventitious-pipeline/mp_metagenomic_assessment_v4.py -d /path/to/data/  
    -t 10 -ci /location/of/Centrifuge_libraries/ -c -qf 
    -bl "v_f_b" -m -rs 1 
    -fd "/adventitious-pipeline/configs/viral_examples.txt" 
    -hn "TC,Jurkat"
```
<br /><br />

## Read quality analysis
* Run [Nanoplot](https://github.com/wdecoster/NanoPlot) for list of samples from config file.
* Extract quality scores for barcoded samples from aggregate fastq read file.
* Concatenate NanoStats.txt from all samples into single csv file.
<br />
 	
~~~~
/adventitious-pipeline/run_nanostat_analyses.py -d /path/to/data/ -t 10 -e "/adventitious-pipeline/configs/viral_examples.txt" 
~~~~

<br /><br />

## Machine learning pipeline
* Data & data cleaning
    * Only metagenomic classification data from centrifuge and HS-BLASTn are used.
    * Read quality statistics are incorporated during initial metagenomic pipeline run.
    * Overall experiment run data appended (from NanoStats.txt).
    * Host read analysis, classified and unclassified read counts appended.
    * HS-BLASTn reads are filtered for predictions with max percent identity less than 83%; this value was found to be below the minimum value for a true positive spike species
    * Centrifuge reads are filtered using a minimum classification mean score greater than 900.
    * This reduces a large number of off-target predictions without filtering out potential adventitious agents.
    * Columns containing standard deviation with NaNs are dropped.
    * Rows with NaN value for read quality are dropped.
    * Databases used are defined by a BLAST:centrifuge dictionary key:value pair: 
    ```
    database_dict = {"v_f_b":"v_f_b"}
    ```
<br />

* Data split - train, test and evaluation datasets
    * Initial data input are split based on user defined samples to keep back.
    * Currently, this is defined using the batch number for the experiment run e.g. `[83,815,26]`.
    * Later during model building, the training data (i.e. not the data kept back) will be split using [sklearn `train_test_split()`](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.train_test_split.html).
<br />

* Label choice
    * Binary classifiers were designed to handle the classification of the subset adventitious agent candidates into (1) contaminated sample and (2) genuine contaminant predictions.
    * Label encoding was as follows:
    * (1) `True_positive: 1, True_negative: 0`.
    * (2) `False_positive: 0, True_positive: 1`.
    * XGBoost Classifier model was used for training, testing, and evaluation.
    * Quality control is assessed using receiver operating characteristic curves alongside cross validation (`cv=5, scoring='accuracy'`), confusion matrices and classification reports. 
<br />

* XGBoost Classifer - Binary classification - two questions posed:
    * Is the sample contaminated?
     ```
     XGBClassifier(random_state = 736, 
                                    alpha = 0, 
                                    gamma = 0, 
                                    max_depth = 5)   


    ```
    <br />

    * Is the predicted contaminant a genuine adventitious agent?
    ```
    XGBClassifier(random_state = 736, 
                                    alpha = 0, 
                                    gamma = 0.001, 
                                    max_depth = 4)    


    ```
    * Generated models are saved using pickle and reloaded for evaluation data.
    * Important features from featurewiz are similarly saved during training and imported for evaluation.
    * Standard scalers are saved from the training step and imported into the evaluation model to regularise the dataset.
    <br />
* Decision matrix
    *  A decision matrix is used to generate a prediction for any given species
    *  The individual species predictions are combined into a pivot table of counts for false positives, true positives, false negatives and true negatives.
    *  A sample level decision is made to predict if a sample is sterile or contaminated using the pivot table counts for the model predictions.
    <br />
* Evaluation
    *  The process described for training and testing datasets is replicated for the evaluation dataset, which are unseen by the model.
    *  Run the candidates through the XGBoost Classifier model.
    *  Important features are imported for the model.
    *  The standard scaler and machine learning model are imported to run on the evaluation dataset.
    *  Evaluation predictions are compared against the true mask for these samples to score accuracy.

<br /><br />




# Pipeline inputs
* config files - non-barcoded samples (file.txt); file format, six fields (_ delimited):
   * date
   * nucleic-acid
   * sample-name
   * concentration-if-known
   * batch
   * duration-of-sequencing
```
20210727_DNA_TCCacnes_1000CFU_83_18
20211021_DNA_Klebpneu-zymo-16S_10CFU_89_18
```
* config files - barcoded samples (file.csv); file format, important fields:
   * date
   * nucleic-acid
   * sample-name
   * concentration-if-known
   * batch
   * duration-of-sequencing
   * barcode
   * kit
   * identifier (should be the same format as batch e.g. SxBx where batch 818 = S8B18)
```
date	NA	strain	concentration_CFU	batch	duration_h	Flowcell_no	Sample_original	Identifier	Barcode	Total_cycles	Spike_species	Kit
20200324	aDNA	PAO1unfiltTCPA	510000000CFU	1	24	4	Unfiltered 1:1 TC/PA	S1B1	barcode01	25	Pseudomonas aeruginosa	RAB204
20200324	aDNA	PAO1unfiltTCPA	10000000CFU	1	24	4	Unfiltered 1:50 TC/PA	S1B1	barcode02	25	Pseudomonas aeruginosa	RAB204
20220221	DNA	CellfreeMedia-1	0CFU	816	18	33	CellfreeMedia-1	S8B16	barcode05	0	Unknown	RBK004

```
* NanoStats.txt
```
General summary:        
Active channels:                35.0
Mean read length:              324.9
Mean read quality:               9.9
Median read length:            245.0
Median read quality:            10.1
Number of reads:               822.0
Read length N50:               304.0
Total bases:               267,083.0
```
* sequencing_summary.txt
```
filename_fastq	filename_fast5	parent_read_id	read_id	run_id	channel	mux	start_time	duration	num_events	passes_filtering	template_start	num_events_template	template_duration	sequence_length_template	mean_qscore_template	strand_score_template	median_template	mad_template	pore_type	experiment_id	sample_id	end_reason
FAP19063_pass_4f5c1a67_0.fastq.gz	FAP19063_pass_4f5c1a67_0.fast5	163da6a1-c643-4648-b422-3e1770f09a91	163da6a1-c643-4648-b422-3e1770f09a91	4f5c1a6754dc1643a5219daf972a9d8caa44220e	362	1	55.725250	0.780750	624	TRUE	55.840250	532	0.665750	274	10.104866	2.963820	88.026276	9.905292	not_set	ws_bd1r3-1_rk_fc_sampled8march	bd1r_3-1_rk	signal_positive
FAP19063_fail_4f5c1a67_0.fastq.gz	FAP19063_fail_4f5c1a67_0.fast5	d5f574d4-bcd6-41b3-93f8-3b2e23caddff	d5f574d4-bcd6-41b3-93f8-3b2e23caddff	4f5c1a6754dc1643a5219daf972a9d8caa44220e	492	1	55.838250	0.638500	510	FALSE	55.858250	494	0.618500	345	6.720924	3.103649	84.288429	9.905292	not_set	ws_bd1r3-1_rk_fc_sampled8march	bd1r_3-1_rk	signal_positive
```

* species dictionary
```
species = {"Pseudomonas aeruginosa":["PA"], "Cutibacterium acnes":["Cacnes","Pacnes"], \
            "Escherichia coli":"EC", "Klebsiella pneumoniae":"Klebpneu", \
             "Candida albicans":"Calbicans", "Staphylococcus aureus":"Saureus", \
              "Bacillus subtilis": "Bsubtilis", "Minute virus":"MVM", \
               "Feline leukemia":"FELV", "Porcine circovirus": "PCV"}
```

<br />
