# RNA2M
Scripts used in [Chialva et al. (2020)](https://doi.org/10.3390/microorganisms8010038) to analyze tomato meta-transcriptome using host RNA-seq data. If you re-use any or part of this code, please reference and cite our paper.
Raw sequence data to run the script are available at NCBI-SRA under accession PRJNA422090. Intermediate results are here provided (raw counts tables) to run downstream analysis.

## Scripts used to trim, filter and taxonomically annotate raw RNA-seq adapter-free reads (50bp, SE)

[1.meta-transcriptome_reads_filtering_annotation.R](https://github.com/mchialva/RNA2M/blob/master/1.meta-transcriptome_reads_filtering_annotation.R)

[2.reads_counting_statistics_normalization.R](https://github.com/mchialva/RNA2M/blob/master/2.reads_counting_statistics_normalization.R)

[Create_taxoner_NCBI_index.sh](https://github.com/mchialva/RNA2M/blob/master/Create_taxoner_NCBI_index.sh)

## Scripts used for meta-transcriptome diversity analysis

[3.bacteria_diversity.R](https://github.com/mchialva/RNA2M/blob/master/3.bacteria_diversity.R)

[4.fungi_diversity.R](https://github.com/mchialva/RNA2M/blob/master/4.fungi_diversity.R)

## Scripts to functionally annotate reads

[5.eggNOG_functional_annotation.R](https://github.com/mchialva/RNA2M/blob/master/5.eggNOG_functional_annotation.R)

## Variance Partitioning Analyses

[6.host_variance_partitioning.R](https://github.com/mchialva/RNA2M/blob/master/6.host_variance_partitioning.R)

## AM colonization analysis

[7.AM_fungi.R](https://github.com/mchialva/RNA2M/blob/master/7.AM_fungi.R)

## plots and graphics

[8.plots.R](https://github.com/mchialva/RNA2M/blob/master/8.plots.R)

[plotting_styles.R](https://github.com/mchialva/RNA2M/blob/master/plotting_styles.R)

## Intermediate result files

[results.tar.gz](https://github.com/mchialva/RNA2M/blob/master/results.tar.gz)

The archive contains raw meta-transcriptome counts (taxonomy and COGs), sample meta-data mapping file, AM colonization data and plant host transcriptome normalized counts

