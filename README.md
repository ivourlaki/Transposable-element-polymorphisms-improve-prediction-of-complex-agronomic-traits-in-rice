# Transposable element polymorphisms improve prediction of complex agronomic traits in rice

Transposable Elements Polymorphisms (TIPs) are significant sources of genetic variation. Previous work has shown that TIPs can improve detection of causative loci on agronomic traits in rice. Here, we quantify the fraction of variance explained by Single Nucleotide Polymorphisms (SNPs) compared to TIPs, and we explore whether TIPs can improve prediction of phenotypes when compared to using only SNPs. We used eleven traits of agronomic relevance from by five different rice population groups (Aus, Indica, Aromatic, Japonica and Admixed), 738 varieties in total. We assess prediction by applying data split validation in two scenarios. In the within population scenario, we predicted performance of improved Indica varieties using the rest of Indica and additional samples. In the across population scenario, we predicted all Aromatic and Admixed samples using the rest of populations. In each scenario, Bayes C and a Bayesian reproducible kernel Hilbert space regression were compared. We find that TIPs can explain an important fraction of total genetic variance, often more than the fraction explained by SNPs, and that they also improve genomic prediction, especially in the across population prediction scenario, where TIPs outperformed SNPs in nine out of the eleven traits analyzed. In some phenotypes like leaf senescence or grain width, using TIPs increased predictive correlation by 40%. Our results evidence, for the first time, that TIPs genotyping can improve prediction on complex agronomic traits in rice, especially when samples to be predicted are less related to training samples. 

**Citation**:  Ioanna-Theoni Vourlaki, Raúl Castanera, Sebastián E. Ramos-Onsins, Josep M. Casacuberta, Miguel Pérez-Enciso. Transposable element polymorphisms improve prediction of complex agronomic traits in rice. Submitted.



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
   
