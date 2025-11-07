---
title: "README"
output: github_document
---

# ComBatSuite
## Introduction
ComBatSuite is built upon the original sva package found on Bioconductor.
ComBatSuite expands the ComBat tools available to various data sets other than
sequencing data (including methylation data, single cell data, and microbiome
data). 

The `ComBatSuite` package contains functions for removing batch effects and
other unwanted variation in high-throughput experiments. Specifically, the 
`ComBatSuite` package contains functions for identifying and building surrogate 
variables for high-dimensional data sets. Surrogate variables are covariates
constructed directly from high-dimensional data (like gene expression/RNA
sequencing/methylation/brain imaging data) that can be used in subsequent
analyses to adjust for unknown, unmodeled, or latent sources of noise. 

The `ComBatSuite` package can be used to remove artifacts in two ways:

1. Identifying and estimating surrogate variables for unknown sources of 
variation in high-throughput experiments
2. directly removing known batch effects using ComBat [@johnson:2007aa] and its
various versions. 

## Installation
### Github Version
To install the most up-to-date version of ComBatSuite, please install directly
from github. You will need the devtools package. You can install both of these
with the following commands:

``` r
if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools")   
}
library(devtools)
install_github("wejlab/ComBatSuite")
```

## Citation
The ComBatSuite package includes multiple different methods created by different
faculty and students, including the original sva package. It would really help
them out if you would cite their work when you use this software. 

### Overall Package
To cite the overall ComBatSuite package please cite the original sva package
and the updated ComBatSuite package: 

* Leek JT, Johnson WE, Parker HS, Jaffe AE, and Storey JD. (2012) The sva 
package for removing batch effects and other unwanted variation in
high-throughput experiments. Bioinformatics DOI:10.1093/bioinformatics/bts034

### Specific Methods 
When using specific methods, please cite those methods as well:

#### sva method
For sva please cite: 

* Leek JT and Storey JD. (2008) A general framework for multiple testing
dependence. Proceedings of the National Academy of Sciences , 105: 
18718-18723. 
* Leek JT and Storey JD. (2007) Capturing heterogeneity in gene expression 
studies by `Surrogate Variable Analysis'. PLoS Genetics, 3: e161.

#### ComBat method
For ComBat please cite: 

* Johnson WE, Li C, Rabinovic A (2007) Adjusting batch effects in microarray 
expression data using empirical Bayes methods. Biostatistics,  8 (1), 118-127

#### mean-only or reference-batch ComBat method
For mean-only or reference-batch ComBat please cite:

* Zhang, Y., Jenkins, D. F., Manimaran, S., Johnson, W. E. (2018). Alternative
empirical Bayes models for adjusting for batch effects in genomic studies. BMC
bioinformatics, 19 (1), 262.

#### ComBat-Seq method

* Zhang, Y., Parmigiani, G., Johnson, W. E., ComBat-seq: batch effect adjustment
for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3,
September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078

#### svaseq method
For svaseq please cite: 

* Leek JT (2014) svaseq: removing batch and other artifacts from count-based
sequencing data.  Nucleic Acids Res. doi: 10.1093/nar/gku864.

#### supervised sva method
For supervised sva please cite: 

* Leek JT (2014) svaseq: removing batch and other artifacts from count-based
sequencing data. Nucleic Acids Res. doi: 10.1093/nar/gku864.
* Gagnon-Bartsch JA, Speed TP (2012) Using control genes to correct for unwanted
variation in microarray data. Biostatistics 13:539-52. 

#### fsva method
For fsva please cite:

* Parker HS, Bravo HC, Leek JT (2013) Removing batch effects for prediction
problems with frozen surrogate variable analysis arXiv:1301.3947

#### psva method
For psva please cite: 

* Parker HS, Leek JT, Favorov AV, Considine M, Xia X, Chavan S, Chung CH, 
Fertig EJ (2014) Preserving biological heterogeneity with a permuted surrogate
variable analysis for genomics batch correction Bioinformatics doi:
10.1093/bioinformatics/btu375
