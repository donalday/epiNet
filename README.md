# epiNet: a convolutional neural network-based regression model to infer the epigenetic crosstalk responsible for CG methylation patterns

epiNet is a convolutional neural network-based regression model to infer the epigenetic crosstalk responsible for CG methylation patterns (including histone modifications). This tool includes scripts to process raw input and output datasets, and automate the training and prediction steps. It is strict forward and compatible with standard bioinformatics files.

![alt text](https://github.com/donalday/epiNet/blob/master/epiNet_model_schematic.jpg?raw=true)

## Main features of epiNet

* Standard [bedGraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) input formats, handling mean CG methylation levels and FPKM values of epigenetic modifications with a defined bin size (for example, 50-kb)
* Simple script to handle all raw input and output data processing
* Flexibile numbers of input features
* Outputs of training and prediction are visualization-friendly with genome browser like [IGV](http://software.broadinstitute.org/software/igv/)
* Automatic iteration of all combinations of multiple input features

## Version
2021/06/14: version 2 available
* Added support of GPGPU (tested on Titan RTX and CUDA 10.2)
* Added support of non-mouse genomes for IGV visualization
* Remove temporary files of unsuccessful learning trials (those with the number of training cycle equals to stop point + 1)

2020/10/20: version 1 available

## Installation

epiNet is written in python and is executed from the command line. The package is installation-free. Just download the script-test-package folder. The folder contains both scripts and test run data. Perform the test run according the documentation and modify for your own research/usage.

epiNet requires mainly the following tools to be installed and be available in the `PATH`:
* [python3](https://www.python.org/downloads/release/python-368/) (developed on 3.6.8, 3.6.10)
* [Tensorflow v1](https://www.tensorflow.org) (developed on 1.14.0, 1.15.0)
* [Keras (not tf.keras)](https://keras.io) (developed on 2.2.4, 2.2.5)
* [IGVTools](http://software.broadinstitute.org/software/igv/igvtools_commandline) (developed on 2.3.67)

For python3, make sure you have the following packages installed: csv, numpy, json, subprocess, pandas, itertools, and matplotlib.

## Documentation

Please follows the beginner's guide [here](https://github.com/donalday/epiNet/blob/master/Documentation.md).

## Credits

epiNet was written by Donald, Wan Kin Au Yeung at [Division of Epigenomics and Development, Medical Institute of Bioregulation, Kyushu University, Japan](http://www.bioreg.kyushu-u.ac.jp/labo/epigenome/index_e.html). This work was supported by JSPS KAKENHI grant to Hiroyuki Sasaki and Osamu Maruyama ([JP18H05214](https://kaken.nii.ac.jp/en/grant/KAKENHI-PROJECT-18H05214/)).

## Citation

Au Yeung, W.K., Maruyama, O. & Sasaki, H. A convolutional neural network-based regression model to infer the epigenetic crosstalk responsible for CG methylation patterns. BMC Bioinformatics 22, 341 (2021).
https://pubmed.ncbi.nlm.nih.gov/34162326/

## Licences

epiNet is free, licensed under GNU General Public License v3.0.
