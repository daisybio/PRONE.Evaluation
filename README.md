# Systematic Evaluation of Normalization Approaches in Tandem Mass Tag and Label-Free Protein Quantification Data Using PRONE

This repository contains all scripts that were produced when analyzing and evaluating the performance of different normalization methods on spike-in and real-world proteomics data sets.


## Load data

First of all, you need to download the raw data from [here](https://drive.google.com/drive/folders/1MTuCDTsHWcXuAkJwQwxWRkK8At_4GulY?usp=drive_link). Put this data directory in the main directory of the repository.

## Spike-in Data Sets

To assess the performance of the normalization methods on the spike-in data sets, you need to run the [Spike-in_Pipeline.R](spike-in/Spike-in_Pipeline.R) script. This script produces different folders of .rds objects that are required for the final evaluation of all spike-in data sets. 

After the execution of the [Spike-in_Pipeline.R](spike-in/Spike-in_Pipeline.R), you can run [Spike-in_Evaluation.R](spike-in/Spike-in_Evaluation.R) to evaluate the performance of the normalization methods and produce different types of plots.

To reproduce the results comparing limma to ROTS, use the script [DE_Methods_Evaluation.R](spike-in/DE_Methods_Evaluation.R).

## Real-world Data Sets

The scripts for the real-world datasets are available in the [real-world](real-world/) directory. The scripts are named according to the dataset they are analyzing. To generate the plots of the paper, use the scripts [Batch_Effect_Correction_Paper.R](real-world/Batch_Effect_Correction_Paper.R) and [DE_Analysis_Paper.R](real-world/DE_Analysis_Paper.R).

## PRONE Package and Shiny App

If you want to perform similar analyses on your proteomics dataset, have a look at the PRONE package available on [Biocondcutor](https://bioconductor.org/packages/release/bioc/html/PRONE.html) and [GitHub](https://github.com/daisybio/PRONE). A web interface is accessible at [https://exbio.wzw.tum.de/prone](https://exbio.wzw.tum.de/prone).

## Citation

The preprint is available on BioRxiv: [https://doi.org/10.1101/2025.01.27.634993](https://doi.org/10.1101/2025.01.27.634993).
