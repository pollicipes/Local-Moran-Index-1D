# Local-Moran-Index-1D
## This is the code that accompanies the paper Lopez de Maturana *et al*., (2020).

Download or clone this repo and start by following the script `00_run_LMI.ipynb`. This is an IPython notebook. The files are intended to be open with the Jupyter Notebook or any other \*.ipynb interpreter. The scripts and the data contain the steps to get to the set of 624 SNPs that were selected by the LMI. 

The script `01_GWAS_LMI_integration.ipynb` contains the R reference code (with comments) we used to get the 510 SNPs we discuss in the paper. You can download some sample data from here to use to test the script.

Then, the SNPs will be integrated with HiC data and significant contacts in the `02_TADBit_HOMER_Overlap_SNPs.ipynb` file. This file contains the bash commands to run HOMER and get significant interactions to be crossed with your SNPs.
