# Documentation of the script-test-package

The following is the details for the script-test-package, which contains both scripts and sample data for test run to predict CG methylation levels from H3K36me3 FPKM values in wildtype mouse oocytes at 50-kb resolution. Perform the test run according the documentation below and modify for your own research/usage.

## Folder structure

```
script-test-package
├── epiNet_model.py (epiNet model)
├── data_processing (Folder to process raw input and output datasets)
│   └── data_process.py
├── training (Folder for epiNet training)
│   └── epiNet_training.py
├── prediction (Folder for epiNet prediction)
│   └── epiNet_prediction.py
└── raw_data_files (Folder with all necessary raw files for test run)
    ├── mm10_50kb_nostep.bed
    ├── mm10_bin_start_50000bp_5step.txt
    ├── WT_CG_FGO_Shirane_50000.bedGraph
    └── WT_H3K36me3_FGO_Xu_50000.bedGraph
```

## Step 1:  Data processing

Please modify the header of **data_process.py** for your own data (Default are the values of the test run).
```
bedGraph_list = ['../raw_data_files/WT_CG_FGO_Shirane_50000.bedGraph', '../raw_data_files/WT_H3K36me3_FGO_Xu_50000.bedGraph'] # BedGraph files for output and input features (name ouput bedGraph first)
feature_name = ['WT_CG', 'WT_K36me3'] # Name of output and input features

feature_data_scaling = [1, 1] # Data scaling before upperlimit cutoff in the form of (WT/KO) e.g. 2 for doubling the data
feature_upperlimit = [100, 0] # Upperlimit N. Each feature with range (0, N) will be scaled to (0, 1). If N = 0, N will be determined during processing (N assigned as 95th percentile value).
feature_max = [100, 1] # Maximum values of each feature in the outpur bedGraph file(s) (eventhough currently input files are not processed)

all_bin_bedFile = '../raw_data_files/mm10_50kb_nostep.bed' # Bed file containing all genomic bins
nearby = '../raw_data_files/mm10_bin_start_50000bp_5step.txt' # Tab-delimited table containing start position of all nearby bins
```
**data_process.py** takes up both input (**WT_H3K36me3_FGO_Xu_50000.bedGraph**) and output (**WT_CG_FGO_Shirane_50000.bedGraph**) bedGraph files from the raw_data_files folder and converts them into a single numpy array.
All available bin coordinates of a genome (**mm10_50kb_nostep.bed**) can be easily created using [bedtools makewindow](https://bedtools.readthedocs.io/en/latest/index.html#) command.
In addition, a tab-delimited file containing all zero-based start positions of nearby bins is required (**mm10_bin_start_50000bp_5step.txt**). For example, `chr1    250000  0       50000   100000  150000  200000  300000  350000  400000  450000  500000` means bins surrounding chr1:250,000-300,000 starts from 50,000 to 500,000.

Furthermore, for each input and output feature, data will be scaled between 0 and 1.
For output CG methylation levels, the original data is in range 0 to 100. So upperlimit is known.
For input H3K36me3 ChIP-seq FPKM values, this script scales the data with 0.95 equals to 95% percentile of all FPKM values.
Feeding `feature_upperlimit` with zero value will initiate determination of upperlimit. You can resuse from a previous processing output (**saved_feature_upperlimit.npy**).
If you have different samples for training and prediction (like widltype and mutant/knockout/...), you can apply data scaling before upperlimit determination.

When everything is ready, run the data processing using the script:
```
python data_process.py
```

After data processing, there will be a data array (**saved_input_in_np.npy**), feature properties (**saved_feature_name.npy**,**saved_feature_upperlimit.npy**,**saved_feature_max.npy**), scaled bedGraph files (**WT_K36me3_50000_ori.bedGraph**,**WT_CG_50000_ori.bedGraph**) and their respective tdf files for IGV visualization.

## Step 2:  Model training

Please modify the header of **epiNet_training.py** for your own data (Default are the values of the test run).
```
no_of_threads = 8 # No. of threads available (not CPU cores)
learning_rate = 0.0001 # Learning rate for Nadam
max_epochs = 9999 # Maximum no. of epochs to try even not pleateau yet
stop_point = 20 # No. of cycles to stop when no more reduction in loss seen (Early stopping)
batch_size = 100 #  Batch size
filter_size = 64 # Filter size of the third layer, filter sizes of other layers will be scaled accordingly

model_file = '../epiNet_model.py' # The epiNet model
all_bin_bedFile = '../raw_data_files/mm10_50kb_nostep.bed' # Bed file containing all genomic bins
input_folder = '../data_processing' # Folder containing of processed data files
feature_name = ['WT_CG', 'WT_K36me3'] # Name of output and input features for training (Exact name used during data processing)
feature_of_input = [2] # List of input fearure(s) (1-based order of feature during data processing)
feature_of_output = [1] # List of output fearure (1-based order of feature during data processing)
chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX'] # Desired order of chromosomes, the 1st entry will be used for testing, the 2nd and 3rd will be used for validation, and others will be used for training
```

The lower half specifies the path of epiNet model files (**epiNet_model.py**) and the data processing folder. The training will split the dataset based on chromosomes for training, validation and testing. The first entry in the `chr_list` will be used for testing, the second and third entries will be used for validation, and others will be used for training. Copy all other parameters **as EXACTLY as stated in data_process.py**.

The upper half specifies the training parameters. We recommend the above learning rate and batch size at a starting point. Do yourself a grid search if you want to further optimize. These parameters work on own hands for the mouse genome at 1-kb, 10-kb and 50-kb resolutions.

When everything is ready, start the training using the script:
```
python epiNet_training.py
```

After training, there will be a **all_corrcoef_table.txt** to summarize the results from all combinations of input features (In the test run, there is only one though). You can find individual training results (numpy arrays and pdf graphs). For example, the only input feature of test run (prefix training_1fea_2), where 1fea and 2 mean the number of input features used for this combination and the features (specified in `feature_for_training`). The output bedGraph of each feature from each input combination is also available (e.g. WT_CG_50000_training_1fea_2_epoch00028.bedGraph), with its converted tdf files.

After running, there will be a data array (**saved_input_in_np.npy**), feature properties (**saved_feature_name.npy**,**saved_feature_upperlimit.npy**,**saved_feature_max.npy**), scaled bedGraph files (**WT_K36me3_50000_ori.bedGraph**,**WT_CG_50000_ori.bedGraph**) and their respective tdf files for IGV visualization.


## Step 3:  Output prediction

Please modify the header of **epiNet_prediction.py** for your own data (Default are the values of the test run).
```
no_of_threads = 8 # No. of threads available (not CPU cores)

training_folder = '../training' # Folder containing training results
input_folder = '../data_processing' # Folder containing of processed data files for prediction
feature_name = ['WT_CG', 'WT_K36me3'] # Name of output and input features for prediction (Exact name used during data processing)
```

This part specifies the path of the training folder (**training**). If you uses a different dataset for prediction, specify the correct `input_folder` and `feature_name`.

```
learning_rate = 0.0001 # Learning rate for Nadam
max_epochs = 9999 # Maximum no. of epochs to try even not pleateau yet
stop_point = 20 # No. of cycles to stop when no more reduction in loss seen (Early stopping)
batch_size = 100 #  Batch size
filter_size = 64 # Filter size of the third layer, filter sizes of other layers will be scaled accordingly

model_file = '../epiNet_model.py' # The epiNet model
all_bin_bedFile = '../raw_data_files/mm10_50kb_nostep.bed' # Bed file containing all genomic bins
feature_of_input = [2] # List of input fearure(s) (1-based order of feature during data processing)
feature_of_output = [1] # List of output fearure (1-based order of feature during data processing)
chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX'] # Desired order of chromosomes, the 1st entry will be used for testing, the 2nd and 3rd will be used for validation, and others will be used for training
```
Copy the above parameters **as EXACTLY as stated in epiNet_training.py**.

When everything is ready, start the prediction using the script:
```
python epiNet_prediction.py
```

After training, there will be a **predict_corrcoef_table.txt** to summarize the results from all combinations of input features. Like training, you can find individual training results (numpy arrays, pdf graphs and output bedGraph/tdf files).

