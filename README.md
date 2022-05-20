
# Transposable element polymorphisms improve prediction of complex agronomic traits in rice

Transposon Insertion Polymorphisms (TIPs) are a significant source of genetic variation. Previous work ([Castanera et al., 2021](https://onlinelibrary.wiley.com/doi/10.1111/tpj.15277)) has shown that TIPs can improve detection of causative loci on agronomic traits in rice. Here, we quantify the fraction of variance explained by Single Nucleotide Polymorphisms (SNPs) compared to TIPs, and we explore whether TIPs can improve prediction of phenotypes when compared to using only SNPs. We used eleven traits of agronomic relevance from by five different rice population groups (Aus, Indica, Aromatic, Japonica and Admixed), 738 varieties in total. We assess prediction by applying data split validation in two scenarios. In the within population scenario, we predicted performance of improved Indica varieties using the rest of Indica and additional samples. In the across population scenario, we predicted all Aromatic and Admixed samples using the rest of populations. In each scenario, Bayes C and a Bayesian reproducible kernel Hilbert space regression were compared. We find that TIPs can explain an important fraction of total genetic variance, often more than the fraction explained by SNPs, and that they also improve genomic prediction, especially in the across population prediction scenario, where TIPs outperformed SNPs in nine out of the eleven traits analyzed. In some phenotypes like leaf senescence or grain width, using TIPs increased predictive correlation by 40%. Our results evidence, for the first time, that TIPs genotyping can improve prediction on complex agronomic traits in rice, especially when samples to be predicted are less related to training samples. 

**Citation**:  Ioanna-Theoni Vourlaki, Raúl Castanera, Sebastián E. Ramos-Onsins, Josep M. Casacuberta, Miguel Pérez-Enciso. Transposable element polymorphisms improve prediction of complex agronomic traits in rice. Submitted. (https://assets.researchsquare.com/files/rs-1280195/v1_covered.pdf?c=1645565639)



## Files
  ### Scripts 
  * BayesC_PREDICTION.MODEL.R : An R script for genomic prediction within and across population analysis applying "BayesC".   
  * RKHS_PREDICTION.MODEL.R :   An R script for genomic prediction within and across population analysis applying  "RKHS".
  * RKHS_GENETIC_VARIANCE_INFERENCE.R : An R script for genetic variance inference within and across population applying  "RKHS".
     
   ### Data
   * Accessions_Traits.csv : A csv file of the 11 traits and their corresponding phenotypic values for the 738 accessions. Data transformed as described in manuscript.
   * snps. RData : SNPs matrix in R format.
   * mitedtx_matrix.RData: Merged matrix of MITE and DTX TIPs in R format.
   * rlxrix_matrix.RData: Merged matrix of RLX and RIX TIPs in R format.
   * Additive_Matrix.RData : The three additive-relationship matrices for each marker (SNPs, MITE/DTX, RLX/RIX) to be used in RKHS method script. 
   
