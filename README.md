# Systematic Evaluation of Normalization Approaches in Tandem Mass Tag and Label-Free Protein Quantification Data Using PRONE

This repository contains all scripts that were produced when analyzing and evaluating the performance of different normalization methods on spike-in and real-world proteomics data sets.


## Load data

First of all, you need to download the raw data from [here](https://drive.google.com/drive/folders/1MTuCDTsHWcXuAkJwQwxWRkK8At_4GulY?usp=drive_link). Put this data directory in the main directory of the repository.

## Spike-in Data Sets

To assess the performance of the normalization methods on the spike-in data sets, you need to run the [Spike-in_Pipeline.R](spike-in/Spike-in_Pipeline.R) script. This script produces different folders of .rds objects that are required for the final evaluation of all spike-in data sets. After the execution of the [Spike-in_Pipeline.R](spike-in/Spike-in_Pipeline.R), you can run [Spike-in_Evaluation.R](spike-in/Spike-in_Evaluation.R) to evaluate the performance of the normalization methods and produce different types of plots.

## Real-world Data Sets

.. more to come ..
