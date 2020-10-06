#!/usr/bin/python3.6
# Core epiNet model and data loading scripts

## Get folder name as project name

import os
run_name = os.path.basename(os.getcwd())


## Limit CPU thread use
## inter_op_parallelism_threads: controls the maximum parallel speedup for a single operation (e.g. matrix duplication)
## intra_op_parallelism_threads: multithreaded implementation (e.g. independent in your TensorFlow graph)

import tensorflow
config = tensorflow.compat.v1.ConfigProto(intra_op_parallelism_threads=no_of_threads,inter_op_parallelism_threads=no_of_threads)
session = tensorflow.compat.v1.Session(config=config)


## Common modules, functions, variables

import csv
import json
import numpy as np
import keras
from keras import layers
from keras import backend as K
from keras.models import Model
from sklearn.linear_model import LinearRegression
import gc


## Load processed input data

input_in_np = np.load('saved_input_in_np.npy')
feature_max = np.load('saved_feature_max.npy')
processed_feature_name = np.load('saved_feature_name.npy')
dict_genome_2_row = json.load(open('dict_genome_2_row.json'))
dict_row_2_genome = json.load(open('dict_row_2_genome.json'))
dict_row_2_nearbyrow = json.load(open('dict_row_2_nearbyrow.json'))


## Import and initialize CNN network

if len(feature_name) != len(feature_max):
    print('No. of feature name does not match the no. of feature maximum values in the input data. Please check. Terminate.')
    exit()

if len(feature_name) != len(processed_feature_name):
    print('No. of feature name does not match the no. during data processing. Please check. Terminate.')
    exit()

feature_counter = 0
while feature_counter < len(feature_name):
    if feature_name[feature_counter] != processed_feature_name[feature_counter]:
        print('Feature names do not match the names during data processing. Please check. Terminate.')
        exit()
    feature_counter += 1

def relu_output_range(x):
    return K.relu(x, max_value=1.0)

no_of_data_per_feature = int(input_in_np.shape[1]/len(feature_name))
epi_shape = (len(feature_of_input), no_of_data_per_feature, 1)

input_epi = keras.Input(shape=epi_shape, name='Input_array') # Tensorflow deprecated warning

x = layers.Conv2D(filter_size, (1,no_of_data_per_feature), padding='valid', activation='relu', name='Conv2D')(input_epi)
x = layers.Flatten()(x)
x = layers.Dense((4*filter_size), activation='relu', name='Dense_fin')(x)
x = layers.Dense(len(feature_of_output), activation=relu_output_range, name='Output_relu')(x)

cnn = Model(input_epi, x)
my_optimizer = keras.optimizers.Nadam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=None, schedule_decay=0.004)
cnn.compile(optimizer=my_optimizer, loss='mean_squared_error') # Tensorflow deprecated warning
cnn.summary()


## Process input data and model fitting

chr_count = 0
chr_start_row = []
chr_end_row = []

dict_chr_start_row = {}
dict_chr_end_row = {}

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

input_hsplit = np.hsplit(input_in_np, len(feature_name))
del input_in_np # Manual garbage collection
gc.collect()
featureSel_counter = 0
first_feature = 1
while featureSel_counter < len(feature_of_input):
    if first_feature == 0:
        array_input_h = np.concatenate((array_input_h, input_hsplit[(feature_of_input[featureSel_counter]-1):feature_of_input[featureSel_counter]]), 0)
    else:
        array_input_h = input_hsplit[(feature_of_input[featureSel_counter]-1):feature_of_input[featureSel_counter]]
        first_feature = 0
    featureSel_counter += 1
array_input = np.hstack(array_input_h)
del array_input_h # Manual garbage collection
gc.collect()
featureSel_counter = 0
first_feature = 1
while featureSel_counter < len(feature_of_output):
    if first_feature == 0:
        array_output_h = np.concatenate((array_output_h, input_hsplit[(feature_of_output[featureSel_counter]-1):feature_of_output[featureSel_counter]]), 0)
    else:
        array_output_h = input_hsplit[(feature_of_output[featureSel_counter]-1):feature_of_output[featureSel_counter]]
        first_feature = 0
    featureSel_counter += 1
array_output = np.hstack(array_output_h)
del input_hsplit # Manual garbage collection
del array_output_h
gc.collect()

array_input = array_input.astype('float32')
array_input = array_input.reshape(array_input.shape[0], int(array_input.shape[1]/no_of_data_per_feature), no_of_data_per_feature, 1)

array_output = array_output.astype('float32')
array_output = array_output.reshape(array_output.shape[0], int(array_output.shape[1]/no_of_data_per_feature), no_of_data_per_feature, 1)
array_output = np.split(array_output, no_of_data_per_feature, axis=2)[0]
array_output = array_output.reshape(array_output.shape[0], array_output.shape[1])

filepath = "saved-model-epoch{epoch:05d}.h5" # Change model name every save
callbacks_list = [keras.callbacks.ModelCheckpoint(filepath, monitor='val_loss', save_best_only=False, save_weights_only=True, mode='auto', period=1), keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=2), keras.callbacks.EarlyStopping(monitor='val_loss', mode='auto', verbose=1, patience=stop_point)] # With EarlyStopping


