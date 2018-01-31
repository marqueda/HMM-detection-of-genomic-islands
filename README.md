# HMM-detection-of-genomic-islands

Hidden Markov model (HMM) approach to detect genomic islands, based on outlier probabilites / quantile-values from an outlier analyses (Hofer et al. 2012, Marques et al. 2016a) or from Fst values (Sorio-Carrasco et al. 2014, Marques et al. 2016b). Please cite Hofer et al. 2012 and the respective references above for the use of either method / scripts. Provided are scripts using either a 3-state HMM or a 2-state HMM.

### Requirements:

R-libraries HiddenMarkov, foreach, doParallel, PLIS

### References:
Hofer T, M Foll & L Excoffier (2012) Evolutionary forces shaping genomic islands of population differentiation in humans. BMC Genomics 13:107.

Marques DA, K Lucek, JI Meier, S Mwaiko, CE Wagner, L Excoffier & O Seehausen (2016a) Genomics of rapid incipient speciation in sympatric threespine stickleback. PLOS Genetics 12: e1005887.

Marques DA, K Lucek, MP Haesler, AF Feller, JI Meier, CE Wagner, L Excoffier & O Seehausen (2016b) Genomic landscape of early ecological speciation initiated by selection on nuptial color. Molecular Ecology 26: 7-24.

Sorio-Carrasco V, Z Gompert, AA Comeault, TE Farkas, TL Parchman, JS Johnston, CA Burkle, JL Feder, J Bast, T Schwander, SP Egan, BJ Crespi & P Nosil (2014) Stick insect genomes reveal natural selection's role in parallel speciation. Science 344: 738-742.
