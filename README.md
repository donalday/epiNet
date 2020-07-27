# epiDeep

epiDeep is a convolutional neural network (CNN)-based regression model to predict CG methylation from epigenetic modifications including histone modifications. This tool includes scripts to process raw input and output datasets, and automate the training and prediction steps. It is strict forward and compatible with standard bioinformatics files.

## Main features of epiDeep

* Standard [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) input formats, handling mean CG methylation levels and FPKM values of epigenetic modifications with a defined bin size (for example, 50-kb)
* Single script to handle all raw input and output data processing
* Felxibility in number of input and output features
* Outputs of training and prediction are visualization-friendly with genome browser like [IGV](http://software.broadinstitute.org/software/igv/)
* Automatic iteration of all combinations of multiple input features

## Installation

epiDeep is written in python and is executed from the command line. The package is installation-free. Just download the script-test-package folder. The folder contains both scripts and test run data. Perform the test run according the documentation and modify for your own research/usage.

epiDeep requires mainly the following tools to be installed and be available in the `PATH`:
* [python3](https://www.python.org/downloads/release/python-368/) (developed on 3.6.8, 3.6.10)
* [Tensorflow v1](https://www.tensorflow.org) (developed on 1.14.0, 1.15.0)
* [Keras (not tf.keras)](https://keras.io) (developed on 2.2.4, 2.2.5)
* [IGVTools](http://software.broadinstitute.org/software/igv/igvtools_commandline) (developed on 2.3.67)

For python3, make sure you have the following packages installed: csv, numpy, json, subprocess, pandas, itertools, and matplotlib.

## Documentation

Please follows the beginner's guide [here](https://github.com/donalday/epiDeep/blob/master/Documentation.md).

## Credits

epiDeep was written by Donald, Wan Kin Au Yeung at [Division of Epigenomics and Development, Medical Institute of Bioregulation, Kyushu University, Japan](http://www.bioreg.kyushu-u.ac.jp/labo/epigenome/index_e.html). This work was supported by JSPS KAKENHI grant to Hiroyuki Sasaki ([JP18H05214](https://kaken.nii.ac.jp/en/grant/KAKENHI-PROJECT-18H05214/)).

## Citation

Au Yeung et al., in preparation

## Licences

epiDeep is free, licensed under GNU General Public License v3.0.
