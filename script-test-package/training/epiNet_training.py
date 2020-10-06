#!/usr/bin/python3.6
# Train the epiNet model and iterate for combinations of input features

## Global parameters

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


## Copy necessary training data from input folder

import subprocess # Run linux command inside python

nullseparator = ''
command_to_run = nullseparator.join(['cp ', input_folder, '/*.json ', input_folder, '/*.npy ./'])
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
    run_name = nullseparator.join([run_name, '_', str(len(feature_of_input)), 'fea_'] + [str(i) for i in feature_of_input])
#
#
###### Process input data and model fitting
#
    chr_count = 3 # Use chr other than 1st, 2nd and 3rd for training
    while chr_count < len(chr_list):
        if chr_count == 3:
            x_train = array_input[chr_start_row[chr_count]:(chr_end_row[chr_count]+1)]
            if iList_count == 0:
                y_train = array_output[chr_start_row[chr_count]:(chr_end_row[chr_count]+1)]
        else:
            x_train = np.concatenate((x_train, array_input[chr_start_row[chr_count]:(chr_end_row[chr_count]+1)]))
            if iList_count == 0:
                y_train = np.concatenate((y_train, array_output[chr_start_row[chr_count]:(chr_end_row[chr_count]+1)]))
        chr_count += 1
    x_val = np.concatenate((array_input[chr_start_row[1]:(chr_end_row[1]+1)], array_input[chr_start_row[2]:(chr_end_row[2]+1)]))
    if iList_count == 0:
        y_val = np.concatenate((array_output[chr_start_row[1]:(chr_end_row[1]+1)], array_output[chr_start_row[2]:(chr_end_row[2]+1)]))
    x_test = array_input[chr_start_row[0]:(chr_end_row[0]+1)]
    if iList_count == 0:
        y_test = array_output[chr_start_row[0]:(chr_end_row[0]+1)]
#
    history = cnn.fit(x=x_train, y=y_train, shuffle=True, epochs=max_epochs, batch_size=batch_size, callbacks=callbacks_list, validation_data=(x_val, y_val))
#
#
###### Save final model and display curve of loss and accuracy during training
#
    cycle_trained = len(history.history.get('loss'))
#
    if (int(cycle_trained) == int(stop_point+1)) and (trial_for_this_combination < 10): # If no learning at all (probably gradient vanished), restart this training
        subprocess.call('rm saved-model-epoch*.h5', shell=True) # Remove all saaved models
        trial_for_this_combination += 1
        continue
    else:
        trial_for_this_combination = 0
#
    open(nullseparator.join([run_name, '_total_training_cycle.txt']), 'w').write(str(cycle_trained)) # Save the number of training cycle before no improvement (Total - patience)
    json.dump(list(history.history.get('loss')), open(nullseparator.join([run_name, '_list_loss.json']), 'w'))
    json.dump(list(history.history.get('val_loss')), open(nullseparator.join([run_name, '_list_val_loss.json']), 'w'))
    cnn.save_weights(nullseparator.join([run_name, '_final_model.h5']))
    loss = history.history['loss']
    val_loss = history.history['val_loss']
#
    import matplotlib.pyplot as plt
    nullseparator = '' # Print correlation graphs
    pyplot_color = ['', 'b', 'r', 'm', 'g', 'k', 'c', 'b0', 'g0', 'r0', 'c0', 'm0', 'k0', 'bx', 'gx', 'rx', 'cx', 'mx', 'kx', 'b^', 'g^', 'r^', 'c^', 'm^', 'k^']
    cycle = range(1, (cycle_trained + 1))
#
    a = plt.figure()
    plt.plot(cycle, loss[:cycle_trained], 'bo', label='Training loss')
    plt.plot(cycle, val_loss[:cycle_trained], 'b', label='Validation loss')
    plt.title('Training and validation loss')
    plt.legend()
    plt.show()
    nullseparator = ''
    loss_name_list = [run_name, '-epiNet-training-loss.pdf']
    a.savefig(nullseparator.join(loss_name_list), bbox_inches='tight')
    plt.cla()
#
#
###### Output preidct epigenetic modifications to bedGraph files (training data, validation data and test data)
#
    def cal_all_corrcoef(data_input, data_decoded, feature_no, epoch): # Calculate Pearson product-moment correlation coefficient for each feature
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
    train_corrcoef = [] # Correlation matrix
    val_corrcoef = []
    test_corrcoef = []
#
    no_of_feature_to_predict = len(feature_of_output)
#
    iterate_count = 1
    while iterate_count <= cycle_trained:
        cnn = Model(input_epi, x)
        model_name_list = ['saved-model-epoch', format(iterate_count, '05d'), '.h5']
        cnn.load_weights(nullseparator.join(model_name_list)) # Keras can only load weights but not model (probably bug?)
        cnn.compile(optimizer=my_optimizer, loss='mse') # Tensorflow deprecated warning
        out_bedGraph_temp = [[] for x in range(no_of_feature_to_predict)]
        array_output_predicted = cnn.predict(array_input, batch_size=batch_size)
        y_train_predicted = cnn.predict(x_train, batch_size=batch_size)
        y_val_predicted = cnn.predict(x_val, batch_size=batch_size)
        y_test_predicted = cnn.predict(x_test, batch_size=batch_size)
        train_corrcoef.append(cal_all_corrcoef(y_train, y_train_predicted, no_of_feature_to_predict, iterate_count))
        val_corrcoef.append(cal_all_corrcoef(y_val, y_val_predicted, no_of_feature_to_predict, iterate_count))
        test_corrcoef.append(cal_all_corrcoef(y_test, y_test_predicted, no_of_feature_to_predict, iterate_count))
        row_count_out = 0 # Writing test data (all chr)
        if iterate_count == cycle_trained: # Only output the last bedGraph
            while row_count_out < len(array_output):
                curr_f_count = 0
                while curr_f_count < no_of_feature_to_predict:
                    out_bedGraph_temp[curr_f_count].append(all_bin_dict_nearby.get(dict_row_2_genome.get(str(row_count_out))) + [str(feature_max[curr_f_count]*float(array_output_predicted[row_count_out][curr_f_count]))]) # Generate bin table for each feature
                    curr_f_count += 1
                row_count_out += 1
            curr_f_count = 0
            while curr_f_count < no_of_feature_to_predict:
                csv.writer(open(nullseparator.join([feature_name[curr_f_count], '_', bin_size, '_', run_name, '_epoch', format(iterate_count, '05d'), '.bedGraph']), 'w', newline="\n"), delimiter='\t').writerows(out_bedGraph_temp[curr_f_count])
                curr_f_count += 1
        iterate_count += 1
#
    train_corrcoef = np.array(train_corrcoef) # Save the arrays
    val_corrcoef = np.array(val_corrcoef)
    test_corrcoef = np.array(test_corrcoef)
    np.save(nullseparator.join([run_name, '_train_corrcoef']), train_corrcoef)
    np.save(nullseparator.join([run_name, '_val_corrcoef']), val_corrcoef)
    np.save(nullseparator.join([run_name, '_test_corrcoef']), test_corrcoef)
#
    subprocess.call('rm saved-model-epoch*.h5', shell=True) # Remove all saaved models and only leave the last model (saved separately)
#
    acc1 = plt.figure()
    f_count = 1
    while f_count < len(train_corrcoef[0]):
        plt.plot(cycle, np.moveaxis(train_corrcoef, 0, 1)[f_count], pyplot_color[f_count], label=feature_name[int(feature_of_output[(f_count-1)]-1)])
        f_count += 1
    plt.title('Training Pearson correlation coefficient')
    plt.legend()
    plt.show()
    loss_name_list = [run_name, '-epiNet-train_corrcoef.pdf']
    acc1.savefig(nullseparator.join(loss_name_list), bbox_inches='tight')
    plt.cla()
#
    acc1 = plt.figure()
    f_count = 1
    while f_count < len(val_corrcoef[0]):
        plt.plot(cycle, np.moveaxis(val_corrcoef, 0, 1)[f_count], pyplot_color[f_count], label=feature_name[int(feature_of_output[(f_count-1)]-1)])
        f_count += 1
    plt.title('Validation Pearson correlation coefficient')
    plt.legend()
    plt.show()
    loss_name_list = [run_name, '-epiNet-val_corrcoef.pdf']
    acc1.savefig(nullseparator.join(loss_name_list), bbox_inches='tight')
    plt.cla()
#
    acc1 = plt.figure()
    f_count = 1
    while f_count < len(test_corrcoef[0]):
        plt.plot(cycle, np.moveaxis(test_corrcoef, 0, 1)[f_count], pyplot_color[f_count], label=feature_name[int(feature_of_output[(f_count-1)]-1)])
        f_count += 1
    plt.title('Test Pearson correlation coefficient')
    plt.legend()
    plt.show()
    loss_name_list = [run_name, '-epiNet-test_corrcoef.pdf']
    acc1.savefig(nullseparator.join(loss_name_list), bbox_inches='tight')
    plt.cla()
#
    iList_count += 1 # End of loop


## IGVtools converion for converting all bedGraph files to TDF files

all_bedGraph_list = subprocess.check_output(['ls -1 *.bedGraph'], shell=True).decode('utf-8').splitlines()

iterate_count = 0
while iterate_count < len(all_bedGraph_list):
    command_to_run1 = nullseparator.join(['sort -k 1,1 -k 2,2n -o temp.bedGraph ', all_bedGraph_list[iterate_count]])
    out_name = all_bedGraph_list[iterate_count].split('.bedGraph')[0]
    command_to_run2 = nullseparator.join(['igvtools toTDF -z 10 temp.bedGraph ', out_name, '_z10.tdf mm10'])
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
chr_count = 0
chr_start_row = []
chr_end_row = []
while chr_count < len(chr_list):
    nullseparator = ''
    name = nullseparator.join([chr_list[chr_count], '_'])
    tmp_list = []
    for key in dict_genome_2_row:
        if key.startswith(name):
            tmp_list.append(key)
    chr_start_row.append(dict_genome_2_row.get(tmp_list[0]))
    chr_end_row.append(dict_genome_2_row.get(tmp_list[-1]))
    chr_count += 1
chr_count = 3 # Use chr other than 1st, 2nd and 3rd for training
while chr_count < len(chr_list):
    if chr_count == 3:
        no_train = (chr_end_row[chr_count] + 1) - chr_start_row[chr_count]
    else:
        no_train += (chr_end_row[chr_count] + 1) - chr_start_row[chr_count]
    chr_count += 1
no_val = (chr_end_row[2] + 1) - chr_start_row[2] + (chr_end_row[1] + 1) - chr_start_row[1]
no_test = (chr_end_row[0] + 1) - chr_start_row[0]

del dict_genome_2_row # Manual garbage collection
del chr_start_row
del chr_end_row
gc.collect()

output_table = []
output_table.append(["Run :", os.path.basename(os.getcwd()), "", "No. of combinations :", len(iList_feature_of_input), "", "No. of training/validation/testing data :", no_train, no_val, no_test])
output_table.append([])
output_feature_name = []
fname_count2 = 0
while fname_count2 < len(feature_of_input):
    output_feature_name.append(feature_name[int(feature_of_input[fname_count2]-1)])
    fname_count2 += 1
output_table.append(['No. of features'] + output_feature_name + ['Total epoch', 'Training R', 'Validation R', 'Testing R'])


## Load and summarize corrcoef values
#
iList_count = 0
while iList_count < len(iList_feature_of_input):
    feature_of_input = iList_feature_of_input[iList_count]
#    
    run_name = nullseparator.join([os.path.basename(os.getcwd()), '_', str(len(feature_of_input)), 'fea_'] + [str(i) for i in feature_of_input])
    train_corrcoef = np.load(nullseparator.join([run_name, '_train_corrcoef.npy']))
    val_corrcoef = np.load(nullseparator.join([run_name, '_val_corrcoef.npy']))
    test_corrcoef = np.load(nullseparator.join([run_name, '_test_corrcoef.npy']))
#
    this_feature_name = np.full(len(full_feature_of_input), '')
    fname_counter = 0
    while fname_counter < len(feature_of_input):
        this_feature_name[int(feature_of_input[fname_counter]-len(feature_of_output)-1)] = 'O'
        fname_counter += 1
#
    output_table.append([int(len(feature_of_input))] + this_feature_name.tolist() + [int(train_corrcoef[-1][0]), train_corrcoef[-1][1], val_corrcoef[-1][1], test_corrcoef[-1][1]])
    del train_corrcoef # Manual garbage collection
    del val_corrcoef
    del test_corrcoef
    gc.collect()
#
    iList_count += 1


## Save output table

csv.writer(open(nullseparator.join([os.path.basename(os.getcwd()), '_all_corrcoef_table.txt']), 'w', newline="\n"), delimiter='\t').writerows(output_table)


