# FGTLDA
Perform four-gamete test and linkage disequilibrium analysis on biallelic SNP loci. <br />
FGTLDA.R is a script to conduct four-gamete test and linkage disequilibrium analysis using SNP calling results.

# Prerequisite
library(data.table)

# Input: 
A tab-delimited text file containing genotypes of a population with the sample IDs as the column names and the SNP IDs as the row names. Both file path and file name need to be defined while running the script.

# Output: 
Two files will be generated. A tab-delimited text file with each row listing the counts for genotype combinations of the tested two loci, SNP IDs, normalized D, $r^2$, and p-value. A text file ending with .FDR, containing adjusted p-values, will also be generated.

# Running the script:
Rscript /path/to/FGTLDA.R /path/to/input.txt /path/to/output.txt <br />

Use the data in the example folder for testing purposes.<br />
