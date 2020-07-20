# Documentation

The following is a package containing both scripts and test run data, to predict CG methylation levels from H3K36me3 FPKM values in wildtype mouse oocytes at 50-kb resolution. Run the scripts according the documentation and modify for your own research/usage.

## Structure of the script-test-package

```
-- Parent folder
 |-- epiDeep_model.py (epiDeep model)
 |-- data_processing (Folder to process raw input and output datasets)
   |-- data_process.py
 |-- training (Folder for epiDeep training)
   |-- epiDeep_training.py 
 |-- prediction (Folder for epiDeep prediction)
   |-- epiDeep_prediction.py
 |-- raw_data_files (Folder with all necessary raw files for test run)
   |-- mm10_50kb_nostep.bed
   |-- mm10_bin_start_50000bp_5step.txt
   |-- WT_CG_FGO_Shirane_50000.bedGraph
   |-- WT_H3K36me3_FGO_Xu_50000.bedGraph
```

## Step 1: Data processing

Please modify the header of data_process.py for your own data (Default are the values of the test run).
```
bedGraph_list = ['../raw_data_files/WT_CG_FGO_Shirane_50000.bedGraph', '../raw_data_files/WT_H3K36me3_FGO_Xu_50000.bedGraph'] # BedGraph files for output and input features (name ouput bedGraph first)
feature_name = ['WT_CG', 'WT_K36me3'] # Name of output and input features

feature_data_scaling = [1, 1] # Data scaling before upperlimit cutoff in the form of (WT/KO) e.g. 2 for doubling the data
feature_upperlimit = [100, 0] # Upperlimit N. Each feature with range (0, N) will be scaled to (0, 1). If N = 0, N will be determined during processing (N assigned as 95th percentile value).
feature_max = [100, 1] # Maximum values of each feature in the outpur bedGraph file(s) (eventhough currently input files are not processed)

all_bin_bedFile = '../raw_data_files/mm10_50kb_nostep.bed' # Bed file containing all genomic bins
nearby = '../raw_data_files/mm10_bin_start_50000bp_5step.txt' # Tab-delimited table containing start position of all nearby bins
```
`data_process.py` takes up both input (WT_H3K36me3_FGO_Xu_50000.bedGraph) and output (WT_CG_FGO_Shirane_50000.bedGraph) bedGraph files from the raw_data_files and converts them into numpy arrays.
All available bin coordinates of a genome (mm10_50kb_nostep.bed) can beeasily created using `bedtools makewindow` command.
In addition, a tab-delimited file containing all zero-based start positions of nearby bins is required (mm10_bin_start_50000bp_5step.txt). For example, `chr1    250000  0       50000   100000  150000  200000  300000  350000  400000  450000  500000` means bins surrounding chr1:250,000-300,000 starts from 50,000 to 500,000.

Furthermore, for each input and output feature, data will be scaled to (0,1).
For output CG methylation levels, the original data is in range 0-100. So upperlimit is known.
For input H3K36me3 ChIP-seq FPKM values, this script scales the data with 0.95 equals to 95% percentile of all FPKM values.
Feeding `feature_upperlimit` with zero value will initiate determination of upperlimit. You can resuse from a previous processing output (saved_feature_upperlimit.npy).
If you have different samples for training and prediction (like widltype and mutant/knockout/...), you can apply data scaling before upperlimit determination.

When everything is ready, run the data processing using the script:
```
python data_process.py
```


## Step 2: Model training



## Step 3: Output prediction


