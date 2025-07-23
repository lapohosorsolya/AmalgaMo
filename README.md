# AmalgaMo: DNA Motif Amalgamator

This repository houses AmalgaMo, a tool developed for flexible merging of DNA-binding motifs of any organism, as well as our optimized merged HOCOMOCO human motif set in multiple formats.

**Contents**

- [Overview](#overview)
- [Download](#download)
- [Dependencies](#dependencies)
- [AmalgaMo Usage](#amalgamo-usage)
  - [Input](#input)
  - [Output](#output)
- [Running monaLisa regression using the output of AmalgaMo](#running-monalisa-regression-using-the-output-of-amalgamo)
  - [Converting AmalgaMo output to JASPAR format](#converting-amalgamo-output-to-jaspar-format)
  - [Environment for monaLisa pipeline](#environment-for-monalisa-pipeline)
  - [Running monaLisa's randomized Lasso stability selection](#running-monalisas-randomized-lasso-stability-selection)


## Overview

Given any set of input motifs in HOCOMOCO, JASPAR, or MEME format, AmalgaMo iteratively merges the most similar pairs according to the user-specified parameters, generating position-probability matrices, logos, and merging information for each resultant motif (see [example_output](/example_output)). The default settings have been tuned specifically using the HOCOMOCO v12 human core motif collection for downstream regression-based motif enrichment analysis. However, the five parameters provide ample flexibility for other applications.

AmalgaMo is described in detail in the associated manuscript titled *AmalgaMo: a flexible tool for DNA motif merging* by Orsolya Lapohos and Gregory Fonseca (under review).

## Download

This repository may be downloaded using:

    git clone https://github.com/lapohosorsolya/AmalgaMo

## Dependencies

This repository was written and tested with Python 3.6. All dependencies are listed in the [requirements.txt](/requirements.txt) file. To create a conda environment with the same Python version and dependencies, please enter the repository directory and run: 

    conda create --name amalgamo_env --file requirements.txt

Then, to run the Python programs described below in this environment, you may use the following command line template:

    conda run -n amalgamo_env python <program_name>.py ...

Alternatively,

    <path to amalgamo_env>/bin/python <program_name>.py ...


## AmalgaMo Usage

### Input

The following parameters are **required**:

- `-i` (input directory): the full path to the directory containing the input motif matrices (individual motif files; must be counts or probabilities)
- `-o` (output directory): the directory where the output should be written (will make this directory if it does not exist yet)
- `-f` (motif format code): the file format of the input motifs (lowercase); HOCOMOCO, JASPAR, and MEME format supported

The following parameters are *optional*:

- `-s` (float between 0.5 and 1): the similarity cutoff for merging; default is 0.9
- `-a` (float between 0.5 and 1): the minimum alignment overlap; default is 0.9
- `-m` (int): the maximum allowed length difference between two merged motifs; default is 2
- `-t` (float between 0 and 1): the total information ratio requirement for merging; default is 0.85
- `-r` (int): the difference in number of positions that can be tolerated when comparing high-information windows; default is 0

Here is a command line template using default settings:

    python AmalgaMo.py 
        -i <path to motif directory>
        -o <path to output directory>
        -f <hocomoco|jaspar|meme>

The default settings have been optimized for downstream regression-based motif enrichment analysis using HOCOMOCO v12 human core motifs as input. For other contexts, users may need to modify certain parameter settings. Please see the manuscript for details.

### Output

A successful run of AmalgaMo produces the following files in the specified output directories:

    .
    ├── AmalgaMo_PPMs.npz                   # NumPy archive file with all output PPMs
    ├── AmalgaMo_results.json               # JSON file mapping merged motif names to original motifs
    ├── AmalgaMo_params.json                # JSON file with run parameters
    ├── AmalgaMo_initial_similarities.npy   # NumPy file containing matrix of pairwise similarity scores before merging
    ├── AmalgaMo_initial_names.npy          # NumPy file with ordered motif names before merging
    ├── AmalgaMo_final_similarities.npy     # NumPy file containing matrix of pairwise similarity scores after merging
    ├── AmalgaMo_final_names.npy            # NumPy file with ordered motif names after merging
    ├── merged_logos
    │   ├── merge_1.pdf                     # merged motif information content logo
    │   ├── merge_2.pdf
    │   ├── merge_3.pdf
    │   └── ...
    └── merged_ppms
        ├── merge_1.pfm                     # merged motif PPM in HOCOMOCO (or other specified) format
        ├── merge_2.pfm
        ├── merge_3.pfm
        └── ...

Note that AmalgaMo merged motifs are *always* **position-probability matrices**, regardless of the selected input and output format (counts or probabilities in HOCOMOCO/JASPAR/MEME format).

Example outputs for the optimized merged motif set from the associated manuscript are available in [example_output](/example_output).

## Running monaLisa regression using the output of AmalgaMo

### Converting AmalgaMo output to JASPAR format

An accessory script is available in [extra](/extra) for converting the NumPy archive file output of AmalgaMo into a single JASPAR-formatted file. This converted file can then be used to run the monaLisa pipeine.

    <environment>/bin/python numpy_to_jaspar.py
        -i <path to AmalgaMo output directory>/AmalgaMo_PPMs.npz 
        -o <path to AmalgaMo output directory>/AmalgaMo_motifs.jaspar

If you wish to skip the previous steps, the AmalgaMo-merged HOCOMOCO motif set is available at [formatted_merged_hocomoco_motifs/AmalgaMo_motifs.jaspar](/formatted_merged_hocomoco_motifs/AmalgaMo_motifs.jaspar).

The MEME equivalents of the conversion script and merged motif set are also available in this repository.

### Environment for monaLisa pipeline

To create a separate R environment for the monaLisa pipeline:

    conda create -n Renv-monalisa -c bioconda r-base=4.4 -y
    conda install r::r-essentials
    conda activate Renv-monalisa

Then, open an R session and install the following:

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("monaLisa")
    BiocManager::install("genomation")
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

### Running monaLisa's randomized Lasso stability selection

Now, the script is ready to be executed as follows:

    Rscript run_monaLisa_regression.R
        <path to the tab-separated file with differential chromatin accessibility>
        <path to output directory>/monalisa_selected_motifs.txt
        <path to AmalgaMo output directory>/AmalgaMo_motifs.jaspar
        <index of the column containing the log2 fold change in accessibility>

Remarks:

- This script is designed for the human genome (hg38), but it may be modified easily for other genomes.
- This script follows the corresponding [Bioconductor vignette](https://bioconductor.org/packages/3.22/bioc/vignettes/monaLisa/inst/doc/selecting_motifs_with_randLassoStabSel.html), counting motif hits scoring above the 85th percentile and setting the Lasso stability selection cutoff to 0.8.
- Since AmalgaMo outputs position-probability matrices, care must be taken when reading them in subsequent analyses. This script, along with the previous conversion script, ensures that they are read appropriately.

# Notes

For any questions please contact orsolya.lapohos@mail.mcgill.ca