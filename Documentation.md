# Documentation

The following is a package containing both scripts and test run data, to predict CG methylation levels from H3K36me3 FPKM values in wildtype mouse oocytes. Please try the test run first and then modify the scripts for your own research/usage.

## Structure of the script-test-package

```
-- Parent folder
 |-- epiDeep_model.py (epiDeep model)
 |-- **data_processing** (Folder to process raw input and output datasets)
   |-- data_process.py
 |-- **training** (Folder for epiDeep training)
   |-- epiDeep_training.py 
 |-- **prediction** (Folder for epiDeep prediction)
   |-- epiDeep_prediction.py
 |-- **raw_data_files** (Folder with all necessary raw files for test run)
   |-- mm10_50kb_nostep.bed
   |-- mm10_bin_start_50000bp_5step.txt
   |-- WT_CG_FGO_Shirane_50000.bedGraph
   |-- WT_H3K36me3_FGO_Xu_50000.bedGraph
```

## Step 1: Data processing

## Step 2: Model training

## Step 3: Output prediction

