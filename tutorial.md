# AmalgaMo tutorial with example

This is an example of how you can download and run AmalgaMo. Full documentation is available in the [README](/README.md) file.

## Step 1: download AmalgaMo and set up its environment

Download this repository using the following command:

    git clone https://github.com/lapohosorsolya/AmalgaMo

Enter the AmalgaMo directory and then set up and activate the environment:

    cd AmalgaMo

    conda create --name amalgamo_env --file requirements.txt

    conda activate amalgamo_env

## Step 2: obtain a motif set in HOCOMOCO/JASPAR/MEME format

AmalgaMo takes PCM/PFM/PPM motif matrices in HOCOMOCO/JASPAR/MEME format. These may be provided as individual files within a directory or they may all be in a single file. Examples of each format (single file for each) are provided in the [example_input](/example_input/) directory. Here, we'll use the HOCOMOCO v12 human core motifs concatenated into one file.

## Step 3: run AmalgaMo with default parameters

We'll run AmalgaMo with default parameters for now, as this setting is likely appropriate for most cases. See the documentation and corresponding manuscript for details on parameters.

    python AmalgaMo.py \
    -i example_input/HOCOMOCO_motifs.pcm \
    -o merged_motifs \
    -f hocomoco

The output of the above command will be in a newly created directory called `merged_motifs`, and its contents will look like those of the [example_output](/example_output/) directory.

## Step 4: convert the merged motifs to JASPAR format

Notice that the output of AmalgaMo gives motifs in the same format as the input motifs. Since we used HOCOMOCO formatted input motifs in the last step, we got HOCOCMOCO motifs as output in [example_output/AmalgaMo_PPMs.pfm](/example_output/AmalgaMo_PPMs.pfm). For downstream applications (such as Step 5 below), other formats may be needed. Here, we'll convert the output of AmalgaMo to JASPAR format. Alternatively, we could have run AmalgaMo using pre-formatted motifs, such as those in [example_input/HOCOMOCO_motifs.jaspar](/example_input/HOCOMOCO_motifs.jaspar).

Accessory scripts are provided for converting the NumPy-formatted output of AmalgaMo into HOCOMOCO/JASPAR/MEME format in the [motif_conversion](/motif_conversion/) directory. Use the following command to convert to the output to JASPAR:

    cd motif_conversion

    python numpy_to_jaspar.py \
    -i ../merged_motifs/AmalgaMo_PPMs.npz \
    -o ../merged_motifs/AmalgaMo_motifs.jaspar

At this point, you are ready to use your merged motifs. When you use this motif set in any downstream application, be sure to use the `AmalgaMo_results.json` file to map the merged motifs to the original motif names.

## [Optional] Step 5: run monaLisa's randomized Lasso stability selection

Let's use the AmalgaMo-merged motifs to find transcription factors that may be responsible for differential accessibility in stimulated versus unstimulated CD4 T cells (dataset from Pahl et al. 2024; see manuscript for details). Processed data are available in the [sample_data_monaLisa](/sample_data_monaLisa/) directory.

First, deactivate the AmalgaMo environment and go up a directory:

    conda deactivate

    cd ..

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

After exiting the R session, make a new directory for the results of monaLisa Lasso stability selection:

    mkdir monaLisa_results

Now, you are ready to run the regression-based motif enrichment analysis:

    Rscript run_monaLisa_regression.R \
    sample_data_monaLisa/differential_accessibility.txt \
    monaLisa_results/monaLisa_selected_motifs.txt \
    merged_motifs/AmalgaMo_motifs.jaspar \
    26

Note that the last argument is the (1-based) index of the log2 fold change column in the differential accessibility file. This may change based on the number of ATAC-seq replicates used (and of course, the software used to generate it).

The output file will look like the one given in the [sample_data_monaLisa](/sample_data_monaLisa/) directory. "fracGC" and "oeCpG" are additional features added to the regression, as recommended by the monaLisa [Bioconductor vignette](https://bioconductor.org/packages/3.22/bioc/vignettes/monaLisa/inst/doc/selecting_motifs_with_randLassoStabSel.html).