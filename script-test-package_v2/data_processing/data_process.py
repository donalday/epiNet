#!/usr/bin/python3.6
# Load bedGraph files and process for epiNet training/prediction

## Input parameters

bedGraph_list = ['../raw_data_files/WT_CG_FGO_Shirane_50000.bedGraph', '../raw_data_files/WT_H3K36me3_FGO_Xu_50000.bedGraph'] # BedGraph files for output and input features (name ouput bedGraph first)
feature_name = ['WT_CG', 'WT_K36me3'] # Name of output and input features

feature_data_scaling = [1, 1] # Data scaling before upperlimit cutoff in the form of (WT/KO) e.g. 2 for doubling the data
feature_upperlimit = [100, 0] # Upperlimit N. Each feature with range (0, N) will be scaled to (0, 1). If N = 0, N will be determined during processing (N assigned as 95th percentile value).
feature_max = [100, 1] # Maximum values of each feature in the outpur bedGraph file(s) (eventhough currently input files are not processed)

species = 'mm10' # Genome assembly of the species, or .chrom.size for unsupported genomes
all_bin_bedFile = '../raw_data_files/mm10_50kb_nostep.bed' # Bed file containing all genomic bins
nearby = '../raw_data_files/mm10_bin_start_50000bp_5step.txt' # Tab-delimited table containing start position of all nearby bins


## Common modules, functions, variables

if len(bedGraph_list) != len(feature_upperlimit):
    print('No. of feature list does not match the no. of feature upperlimit values in the input data. Please check. Terminate.')
    exit()
if len(feature_name) != len(feature_max):
    print('No. of feature name does not match the no. of feature maximum values in the input data. Please check. Terminate.')
    exit()

import csv
import numpy as np
import json
import pandas as pd
import subprocess

ulseparator = '_'

def create_data_w_nearby(input):
    output = []
    for line in input:
        output.append([ulseparator.join([line[0], line[1]]), line[3], '-'])
    return output

def cal_data_w_nearby(input, pos): # Position of the value of current bin (zero-based)
    input_dict_nearby = {} # Initialize dictionary for current genome position to current table value
    for line in input:
        input_dict_nearby[line[0]] = line[pos]
    output = []
    for line in input:
        if line[0] in all_nearby_dict_nearby:
            nearby_list1 = all_nearby_dict_nearby.get(line[0])
            total1 = 0.0
            eff_count1 = 0.0 # float
            average1 = 0.0
            for element in nearby_list1:
                if element in input_dict_nearby:
                    total1 += float(input_dict_nearby.get(element))
                    eff_count1 += 1.0
            if eff_count1 > 0.0:
                average1 = total1 / eff_count1
                output_line = []
                output_count = 0
                while output_count < len(line):
                    if output_count == pos:
                        output_line.append(float(line[pos]))
                        output_line.append(average1)
                        output_count += 2
                    else:
                        output_line.append(line[output_count])
                        output_count += 1
                output.append(output_line)
    return output

def std_single_data(input_data, upperlimit, noise_level_inner, pos): # Standardize input data using knwon upperlimit if available
    temp = input_data
    for line in temp:
        line[pos] = float(line[pos]) * noise_level_inner
    if upperlimit == 0:
        pos_cal = []
        for line in temp:
            pos_cal.append(float(line[pos]))
        upperlimit = (1.00/0.95) * np.percentile(pos_cal, 95) # Set 95% Percentile as 0.95 (max is 1)
    for line in temp:
        if float(line[pos]) < upperlimit:
            line[pos] = float(line[pos]) / upperlimit
        else:
            line[pos] = 1
    return [temp, upperlimit]

def join_2D_list(list1, list2):
    a = pd.DataFrame(list1)
    b = pd.DataFrame(list2)
    joined = np.array(pd.merge(a, b, on=0, how='inner')).tolist()
    return joined

def join_multi_2D_list(list_of_list):
    joined = list_of_list[0]
    count = 1
    while count < len(list_of_list):
        joined = join_2D_list(joined, list_of_list[count])
        count += 1
    return joined


## Import bin table (https://stackoverflow.com/questions/10937918/loading-a-file-into-a-numpy-array-with-python)

all_nearby = [row for row in (csv.reader(open(nearby, 'rt', newline="\n"), delimiter='\t'))]
all_nearby_dict_nearby = {}
count = 0
while count < len(all_nearby):
    nearby_ele = all_nearby[count][2:]
    ele_count = 0
    while ele_count < len(all_nearby[count][2:]):
        nearby_ele[ele_count] = ulseparator.join([all_nearby[count][0], all_nearby[count][ele_count+2]])
        ele_count += 1
    all_nearby_dict_nearby[ulseparator.join([all_nearby[count][0], all_nearby[count][1]])] = nearby_ele
    count += 1


## Load input data, standardize to range(0,1) and merge as a combined table in the form of a numpy array

loaded_input = [[] for x in range(len(bedGraph_list))]
feature_count = 0
while feature_count < len(bedGraph_list):
    [loaded_input[feature_count], feature_upperlimit[feature_count]] = std_single_data(create_data_w_nearby([row for row in (csv.reader(open(bedGraph_list[feature_count], 'rt', newline="\n"), delimiter='\t'))]), feature_upperlimit[feature_count], feature_data_scaling[feature_count], 1)
    loaded_input[feature_count] = cal_data_w_nearby(loaded_input[feature_count], 1)
    feature_count += 1

joined = join_multi_2D_list(loaded_input)

dict_genome_2_row = {} # Initialize dictionaries for interconverting genome position and row number
dict_row_2_genome = {}
input_in_np = [] # Actual input in the form of numpy array
con_count = 0
while con_count < len(joined):
    dict_genome_2_row[joined[con_count][0]] = con_count
    dict_row_2_genome[int(con_count)] = joined[con_count][0]
    input_in_np.append(joined[con_count][1:])
    con_count += 1
input_in_np = np.array(input_in_np)

dict_row_2_nearbyrow = {} # Initialize dictionary for getting nearby genome positions (in the form of row no.)
near_count = 0
while near_count < len(input_in_np):
    p_pos = dict_row_2_genome.get(near_count)
    near_g_pos1 = all_nearby_dict_nearby.get(p_pos)
    near_row = []
    for pos in near_g_pos1:
        if pos in dict_genome_2_row:
            near_row.append(dict_genome_2_row.get(pos))
    dict_row_2_nearbyrow[near_count] = near_row
    near_count += 1

np.save('saved_input_in_np', input_in_np) # save the arrays
np.save('saved_feature_upperlimit', feature_upperlimit)
np.save('saved_feature_name', feature_name)
np.save('saved_feature_max', feature_max)
json.dump(dict_genome_2_row, open('dict_genome_2_row.json', 'w'))
json.dump(dict_row_2_genome, open('dict_row_2_genome.json', 'w'))
json.dump(dict_row_2_nearbyrow, open('dict_row_2_nearbyrow.json', 'w'))


## Plot original data

no_of_data_per_feature = int(input_in_np.shape[1]/len(bedGraph_list))

ulseparator = '_'
all_bin = [row for row in (csv.reader(open(all_bin_bedFile, 'rt', newline="\n"), delimiter='\t'))]
bin_size = str(int(all_bin[0][2]) - int(all_bin[0][1]))
all_bin_dict_nearby = {} # Initialize dictionary for converting genome position to bed bin data
count = 0
while count < len(all_bin):
    all_bin_dict_nearby[ulseparator.join([all_bin[count][0], all_bin[count][1]])] = all_bin[count]
    count += 1

epi_shape = (int(input_in_np.shape[1]/no_of_data_per_feature), no_of_data_per_feature, 1)
if (epi_shape[0]) != len(feature_name):
    print('No. of feature name does not match the no. of features in the input data. Please check. Terminate.')
    exit()
input_in_np = input_in_np.astype('float32')
input_in_np = input_in_np.reshape(input_in_np.shape[0], int(input_in_np.shape[1]/no_of_data_per_feature), no_of_data_per_feature, 1) # Reshape to 2D arrray with one extra layer

out_bedGraph_temp = [[] for x in range(len(feature_name))]

row_count_out = 0
while row_count_out < input_in_np.shape[0]: # All input data
    curr_f_count = 0
    while curr_f_count < len(feature_name):
        out_bedGraph_temp[curr_f_count].append(all_bin_dict_nearby.get(dict_row_2_genome.get(row_count_out)) + [str(feature_max[curr_f_count]*float(input_in_np[row_count_out][curr_f_count][0]))]) # Generate bin table for each feature
        curr_f_count += 1
    row_count_out += 1

curr_f_count = 0
nullseparator = ''
while curr_f_count < len(feature_name):
    csv.writer(open(nullseparator.join([feature_name[curr_f_count], '_', bin_size, '_ori.bedGraph']), 'w', newline="\n"), delimiter='\t').writerows(out_bedGraph_temp[curr_f_count])
    curr_f_count += 1

curr_f_count = 0
nullseparator = ''
while curr_f_count < len(feature_name):
    command_to_run1 = nullseparator.join(['sort -k 1,1 -k 2,2n -o temp.bedGraph ', feature_name[curr_f_count], '_', bin_size, '_ori.bedGraph'])
    command_to_run2 = nullseparator.join(['igvtools toTDF -z 10 temp.bedGraph ', feature_name[curr_f_count], '_', bin_size, '_ori_z10.tdf', ' ', species])
    command_to_run3 = 'rm temp.bedGraph'
    subprocess.call(command_to_run1.split())
    subprocess.call(command_to_run2.split())
    subprocess.call(command_to_run3.split())
    curr_f_count += 1

