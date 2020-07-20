# epiDeep

epiDeep is a convolutional neural network (CNN)-based regression model to predict CG methylation from epigenetic modifications including histone modifications. This package includes scripts to process raw input and output datasets, and automate the training and prediction steps. It is strict forward and bioinformatics-friendly. 

## Main features of epiDeep

* Datasets implemented as standard [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) formats, handling mean CG methylation levels and FPKM values of epigenetic modifications with a defined bin size (for example, 50 kb)
* Single script to handle all raw input and output data processing
* Felxibility in number of input and output features
* Outputs of training and prediction are visualization-friendly with genome browser like [IGV](http://software.broadinstitute.org/software/igv/)
* Automatic iteration of all combinations of multiple input features

## Installation

epiDeep is written in python and is executed from the command line. The package is installation-free. Just download the package and unzip. The package is the form of a test script scenario. Run the scripts according the documentation and modify for your own research/usage.

epiDeep requires the following tools to be available in the `PATH`:
* [python3](https://www.python.org/downloads/release/python-368/) (developed on 3.6.8)
* [Tensorflow v1](https://www.tensorflow.org) (developed on 1.14.0)
* [Keras (not tf.keras)](https://keras.io) (developed on 2.2.4)
* [IGVTools](http://software.broadinstitute.org/software/igv/igvtools_commandline) (developed on 2.3.67)

For python3, make sure you have the following packages installed: csv, numpy, json, subprocess, pandas, itertools, matplotlib

## Credits

epiDeep was written by Donald, Wan Kin Au Yeung at [Division of Epigenomics and Development, Medical Institute of Bioregulation, Kyushu University, Japan](http://www.bioreg.kyushu-u.ac.jp/labo/epigenome/index_e.html). This work was supported by JSPS KAKENHI grant to Hiroyuki Sasaki ([JP18H05214](https://kaken.nii.ac.jp/ja/grant/KAKENHI-PROJECT-18H05214/)).

## Citation

Please cite 'Deep learning captures epigenetic crosstalk associated with CG methylation' (Au Yeung et al., in preparation).

## Licences

epiDeep is free, licenced under GNU General Public License v3.0.
