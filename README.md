# epiDeep

epiDeep is a convolutional neural network (CNN)-based regression model to predict CG methylation from epigenetic modifications including histone modifications. This package includes scripts to process raw input and output datasets, and automate the training and prediction steps. It is strict forward and bioinformatics-friendly. It's main features are:
. Datasets implemented as standard bedGraph formats, handling mean CG methylation levels and FPKM values of epigenetic modifications with a defined bin size (for example, 50 kb)
. Single script to handle all raw input and output data processing
. Felxibility in number of input and output features
. Outputs of training and prediction are visualization-friendly with genome browser like IGV
. Automatic iteration of all combinations of multiple input features
