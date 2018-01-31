# HMM-detection-of-genomic-islands

Hidden Markov model (HMM) approach to detect genomic islands, based on outlier probabilites / quantile-values from an outlier analyses (Hofer et al. 2012, Marques et al. 2016a) or from log10-transformed Fst values (Marques et al. 2016b). Please cite the references above for the use of either method / scripts. Provided are two scripts each using either a 3-state HMM or a 2-state HMM.

## Requirements:

R-libraries HiddenMarkov, foreach, doParallel, PLIS

## References:
Hofer T, M Foll & L Excoffier (2012) Evolutionary forces shaping genomic islands of population differentiation in humans. BMC Genomics 13:107.

Marques DA, K Lucek, JI Meier, S Mwaiko, CE Wagner, L Excoffier & O Seehausen (2016a) Genomics of rapid incipient speciation in sympatric threespine stickleback. PLOS Genetics 12: e1005887.

Marques DA, K Lucek, MP Haesler, AF Feller, JI Meier, CE Wagner, L Excoffier & O Seehausen (2016b) Genomic landscape of early ecological speciation initiated by selection on nuptial color. Molecular Ecology 26: 7-24.
