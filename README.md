# Transposable element polymorphisms improve prediction of complex agronomic traits in rice
**Scripts and data material related to the paper "Transposable element polymorphisms improve prediction of complex agronomic traits in rice".**

*Ioanna-Theoni Vourlaki, Raúl Castanera, Sebastián E. Ramos-Onsins, Josep M. Casacuberta, Miguel Pérez-Enciso*

## Files
  ### Scripts 
  * BayesC_PREDICTION.MODEL.R : An R script for running genomic prediction for within and across population analysis, applying method "BayesC".   
  * RKHS_PREDICTION.MODEL.R :   An R script for running genomic prediction for within and across population analysis, applying method "RKHS".
  * RKHS_GENETIC_VARIANCE_INFERENCE.R : An R script for genetic variance inference, within and across population, applying method "RKHS".    
   ### Data-Material
   * all.phenotypes.csv : A file of the 11 phenotypes and their corresponding phenotypic values for 738 assecions originated by 5 groups in csv format. 
   * iris_pedigree.csv :  A file of the 11 phenotypes and their corresponding values of their improved varieties, in csv format.
   * snps. RData : A SNPs matrix in R format.
   * Additive_Matrix.RData : The three additive-relationship matrices for each marker (SNPs, MITE/DTX, RLX/RIX) to be used in RKHS method script. 
   * PCAs_fixed_effect.RData: A matrix in R format generating by the three markes to be used as fixed effect in any prediction script.
   
