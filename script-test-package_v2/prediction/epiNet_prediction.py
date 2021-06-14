#!/usr/bin/python3.6
# Use trained epiNet model for prediction, inlcuding iterated combinations of input features

## Global parameters

no_of_threads = 8 # No. of threads available (not CPU cores)
no_of_gpu = 1 # Number of GPU installed (default 1 or 0)

training_folder = '../training' # Folder containing training results
input_folder = '../data_processing' # Folder containing of processed data files for prediction
feature_name = ['WT_CG', 'WT_K36me3'] # Name of output and input features for prediction (Exact name used during data processing)


## Global parameters directly copied from epiNet_training.py

learning_rate = 0.0001 # Learning rate for Nadam
max_epochs = 9999 # Maximum no. of epochs to try even not pleateau yet
stop_point = 20 # No. of cycles to stop when no more reduction in loss seen (Early stopping)
batch_size = 100 #  Batch size
filter_size = 64 # Filter size of the third layer, filter sizes of other layers will be scaled accordingly

model_file = '../epiNet_model.py' # The epiNet model
species = 'mm10' # Genome assembly of the species, or .chrom.size for unsupported genomes
all_bin_bedFile = '../raw_data_files/mm10_50kb_nostep.bed' # Bed file containing all genomic bins
feature_of_input = [2] # List of input fearure(s) (1-based order of feature during data processing)
feature_of_output = [1] # List of output fearure (1-based order of feature during data processing)
chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX'] # Desired order of chromosomes, the 1st entry will be used for testing, the 2nd and 3rd will be used for validation, and others will be used for training


## Copy necessary training data from input folder

import subprocess # Run linux command inside python

nullseparator = ''
command_to_run = nullseparator.join(['cp ', input_folder, '/*.json ', input_folder, '/*.npy ./'])
subprocess.call(command_to_run, shell=True)
command_to_run = nullseparator.join(['cp ', training_folder, '/*_total_training_cycle.txt ', training_folder, '/*.h5 ./'])
subprocess.call(command_to_run, shell=True)


## Generate all combinations out of features for training
#
import itertools
full_feature_of_input = feature_of_input # Save a copy
feature_count = len(full_feature_of_input)
while feature_count > 0:
    if feature_count < len(full_feature_of_input):
        iList_feature_of_input += list(itertools.combinations(feature_of_input, feature_count))
    else:
        iList_feature_of_input = list(itertools.combinations(feature_of_input, feature_count))
    feature_count -=1
#
print("Total combinations", len(iList_feature_of_input))
print("")
#
trial_for_this_combination = 0 # Check try how many times before confirming this combination is not trainable
#
iList_count = 0
while iList_count < len(iList_feature_of_input):
    feature_of_input = iList_feature_of_input[iList_count]
#    
###### Load the model
#
    exec(open(model_file).read())
    dict_row_2_nearbyrow = {} # Manual garbage collection
    gc.collect()
    train_name = nullseparator.join([(training_folder.split("/")[-1]), '_', str(len(feature_of_input)), 'fea_'] + [str(i) for i in feature_of_input])
    run_name = nullseparator.join([run_name, '_', str(len(feature_of_input)), 'fea_'] + [str(i) for i in feature_of_input])
#
#
###### Load info of training result
#
    with open(nullseparator.join([train_name, '_total_training_cycle.txt']), 'r') as file:
        cycle_trained = int(file.read().replace('\n', ''))
#
###### Output preidct epigenetic modifications to bedGraph files (all data)
#
    def cal_all_corrcoef(data_input, data_decoded, feature_no, epoch): # Calculate Pearson product-moment correlation coefficient for each feature
#        print("Pearson product-moment correlation coefficient calculation started") # Debugger
        output = []
        output.append(epoch)
        feature_count = 0
        while feature_count < feature_no:
            this_corrcoef = np.corrcoef(np.moveaxis(data_input, 0, 1)[feature_count].flatten(), np.moveaxis(data_decoded, 0, 1)[feature_count].flatten())
            output.append(this_corrcoef[0][1])
            feature_count += 1
        return output
#
    ulseparator = '_'
    all_bin = [row for row in (csv.reader(open(all_bin_bedFile, 'rt', newline="\n"), delimiter='\t'))]
    bin_size = str(int(all_bin[0][2]) - int(all_bin[0][1]))
    all_bin_dict_nearby = {} # Initialize dictionary for converting genome position to bed bin data
    count = 0
    while count < len(all_bin):
        all_bin_dict_nearby[ulseparator.join([all_bin[count][0], all_bin[count][1]])] = all_bin[count]
        count += 1
#
    predict_corrcoef = [] # Correlation matrix
#
    no_of_feature_to_predict = len(feature_of_output)
#
    cnn = Model(input_epi, x)
    cnn.load_weights(nullseparator.join([train_name, '_final_model.h5'])) # Keras can only load weights but not model (probably bug?)
    cnn.compile(optimizer=my_optimizer, loss='mse') # Tensorflow deprecated warning
#    cnn.summary()
    out_bedGraph_temp = [[] for x in range(no_of_feature_to_predict)]
    array_output_predicted = cnn.predict(array_input, batch_size=batch_size)
    predict_corrcoef.append(cal_all_corrcoef(array_output, array_output_predicted, no_of_feature_to_predict, cycle_trained))
    row_count_out = 0 # Writing test data (all chr)
#
    while row_count_out < len(array_output):
        curr_f_count = 0
        while curr_f_count < no_of_feature_to_predict:
            out_bedGraph_temp[curr_f_count].append(all_bin_dict_nearby.get(dict_row_2_genome.get(str(row_count_out))) + [str(feature_max[curr_f_count]*float(array_output_predicted[row_count_out][curr_f_count]))]) # Generate bin table for each feature
            curr_f_count += 1
        row_count_out += 1
    curr_f_count = 0
    while curr_f_count < no_of_feature_to_predict:
        csv.writer(open(nullseparator.join([feature_name[curr_f_count], '_', bin_size, '_', run_name, '_epoch', format(cycle_trained, '05d'), '.bedGraph']), 'w', newline="\n"), delimiter='\t').writerows(out_bedGraph_temp[curr_f_count])
        curr_f_count += 1
#
    predict_corrcoef = np.array(predict_corrcoef) # Save the arrays
    np.save(nullseparator.join([run_name, '_predict_corrcoef']), predict_corrcoef)
#
    iList_count += 1 # End of loop


## IGVtools converion for converting all bedGraph files to TDF files

all_bedGraph_list = subprocess.check_output(['ls -1 *.bedGraph'], shell=True).decode('utf-8').splitlines()

iterate_count = 0
while iterate_count < len(all_bedGraph_list):
    command_to_run1 = nullseparator.join(['sort -k 1,1 -k 2,2n -o temp.bedGraph ', all_bedGraph_list[iterate_count]])
    out_name = all_bedGraph_list[iterate_count].split('.bedGraph')[0]
    command_to_run2 = nullseparator.join(['igvtools toTDF -z 10 temp.bedGraph ', out_name, '_z10.tdf ', species])
    command_to_run3 = 'rm temp.bedGraph'
    subprocess.call(command_to_run1, shell=True)
    subprocess.call(command_to_run2, shell=True)
    subprocess.call(command_to_run3, shell=True)
    iterate_count += 1


## Generate all combinations out of features and save file header

feature_of_input = full_feature_of_input # Reset again
feature_count = len(full_feature_of_input)
while feature_count > 0:
    if feature_count < len(full_feature_of_input):
        iList_feature_of_input += list(itertools.combinations(feature_of_input, feature_count))
    else:
        iList_feature_of_input = list(itertools.combinations(feature_of_input, feature_count))
    feature_count -=1

dict_genome_2_row = json.load(open('dict_genome_2_row.json'))
no_predict = len(dict_genome_2_row)

del dict_genome_2_row # Manual garbage collection
gc.collect()

output_table = []
output_table.append(["Run :", os.path.basename(os.getcwd()), "", "No. of combinations :", len(iList_feature_of_input), "", "No. of prediction data :", no_predict])
output_table.append([])
output_feature_name = []
fname_count2 = 0
while fname_count2 < len(feature_of_input):
    output_feature_name.append(feature_name[int(feature_of_input[fname_count2]-1)])
    fname_count2 += 1
output_table.append(['No. of features'] + output_feature_name + ['Total epoch', 'Prediction R'])


## Load and summarize corrcoef values
#
iList_count = 0
while iList_count < len(iList_feature_of_input):
    feature_of_input = iList_feature_of_input[iList_count]
#    
    run_name = nullseparator.join([os.path.basename(os.getcwd()), '_', str(len(feature_of_input)), 'fea_'] + [str(i) for i in feature_of_input])
    predict_corrcoef = np.load(nullseparator.join([run_name, '_predict_corrcoef.npy']))
#
    this_feature_name = np.full(len(full_feature_of_input), '')
    fname_counter = 0
    while fname_counter < len(feature_of_input):
        this_feature_name[int(feature_of_input[fname_counter]-len(feature_of_output)-1)] = 'O'
        fname_counter += 1
#
    output_table.append([int(len(feature_of_input))] + this_feature_name.tolist() + [int(predict_corrcoef[-1][0]), predict_corrcoef[-1][1]])
    del predict_corrcoef # Manual garbage collection
    gc.collect()
#
    iList_count += 1


## Save output table

csv.writer(open(nullseparator.join([os.path.basename(os.getcwd()), '_predict_corrcoef_table.txt']), 'w', newline="\n"), delimiter='\t').writerows(output_table)


