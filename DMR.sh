module load combined-pvalues
#https://github.com/brentp/combined-pvalues/tree/master/examples

comb-p pipeline -c 4 --dist 500 --seed 1.0e-4 -p DMR/DMR dt_DMR.bed
