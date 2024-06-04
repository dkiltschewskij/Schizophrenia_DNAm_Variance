# Schizophrenia_DNAm_Variance
Repository of scripts utilised in the manuscript entitled "Schizophrenia is Associated with Altered DNA Methylation Variance"

# gse_series_matrix_to_pheno.R
Commands utilised to call all relevant phenotype information and covariates from each Gene Expression Omnibus gse file.

# process_arrays_stage1.R
Script for processing DNA methylation array data using methods similar to those outlined in (PMID: 33646943). Briefly, this script reads in a matrix of intensities and P values, performs initial sample filtration (i.e. retaining individuals kept using the gse_series_matrix_to_pheno.R script), and covnerts this data into a MethylSet object. It then conducts QC steps outlined in the accompanying manuscript, including: probe/sample P value filtration, removal of samples with low median intensity, removal of outlier samples, sex estimation/filtration, removal of sex chromosomes/probes near common SNPs/crosshybridising probes/non CpG probes. The data is then normalised via Dasen method, and the resultant beta maatrix, probe intensities and phenotype data are returned.

# process_arrays_stage2.R
This script generates methylation derived variables for inclusion as covariates. Using normalised intensities from "process_arrays_stage1.R", this script estimates cell-type proportions, smoking scores and age. It also visualises normalised beta distributions. 

# DMP_VMP_EWAS.R
Systematically analyses DMPs and VMPs using the beta matrix from stage 1 and the phenotype data from stage 2. Analyses DMPs (mean effects) using a linear model regressing probe beta values as the outcome versus diagnostic status + covariates. VMps (variance effects) are identified via Levene's Test (primary model), Bartlett's Test and the Fligner-Kileen Test. Returns summary statistics in "chunks" to reduce memory requirements.

# DMP_VMP_EWAS_by_sex.R
As above, except enabling the user to define which sex they would like to analyse.

# adjust_EWAS.R
Collates all EWAS "chunk" files and adjusts all mean and variance effects for inflation via the bacon R package.

# DMP_VMP_metaanalysis.R
Reads in adjusted EWAS results from multiple studies and meta-analyses all mean effects (IVW with fixed and random effects, using bacon-adjusted results as input) and variance effects (Stouffer's method using Bacon-adjusted P values). Meta-analysed results by chromosome to enable parallelisation and reduce memory requirements.
