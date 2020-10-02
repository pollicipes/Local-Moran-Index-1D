# Local-Moran-Index-1D
## This is the code that accompanies the paper Lopez de Maturana *et al*., (2020).

The scripts and the data contain the reference bioinformatic steps followed to get to the set of 510 SNPs that were selected by the LMI. Code is given in IPython notebooks, for clarity. The files are intended to be open with the Jupyter Notebook or any other \*.ipynb interpreter. 

Download or clone this repo and start by following the script `00_run_LMI.ipynb`. It will run the LMI calculation on the SNPs on your sumstats. Some complementary data, like sumstats or 1000 Genomes CEU genotypes used in the paper can be [downloaded from here](https://drive.google.com/file/d/1VQ5MEnoqEu6a8f0w1xSY3nyxiKI3GoxW/view?usp=sharing) so the LMI can be run.

Next, the script `01_GWAS_LMI_integration.ipynb` contains the R reference code (with comments) we used to get the 510 SNPs we discuss in the paper. You can download some complementary sample data from here to use to test the script. 
Then, the SNPs will be integrated with HiC data and significant contacts in the `02_TADBit_HOMER_Overlap_SNPs.ipynb` file. This file contains the bash commands to run HOMER and get significant interactions to be crossed with your SNPs.

Finally, the script `calculate_credible_sets_LMI.R` was used as an adaptation of https://github.com/hailianghuang/FM-summary to calculate our credible sets, qqplots and enrichments from the summary statistics.
