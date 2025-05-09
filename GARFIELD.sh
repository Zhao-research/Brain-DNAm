
ml garfield/2
ml R

INPUTNAME=Jansen_AD
DATADIR=/path/GARFIELD_Jansen
PLOTDIR=/path/GARFIELD_MAP/
PRUNETAGSDIR=$DATADIR/tags/r01
CLUMPTAGSDIR=$DATADIR/tags/r08
MAFTSSDDIR=$DATADIR/maftssd
PVALDIR=$DATADIR/pval/$INPUTNAME
ANNOTDIR=$DATADIR/annotation
OUTDIR=$DATADIR/output/$INPUTNAME
mkdir -p $OUTDIR

ANNOTLINKFILE=$ANNOTDIR/link_file.txt
PTHRESH=0.1,0.01,0.001,1e-04,1e-05,1e-06,1e-07,1e-08
BINNING=m5,n5,t5
CONDITION=0
SUBSET="1-1005"

F1=$OUTDIR/garfield.prep.$INPUTNAME.out
F0=$OUTDIR/garfield.Meff.$INPUTNAME.out

echo 'Prune and Clump'
echo -n > $F1
for CHR in `seq 1 22` #X
do
	echo 'CHR'$CHR
	${HPC_GARFIELD_BIN}/garfield-prep-chr -ptags $PRUNETAGSDIR/chr$CHR -ctags $CLUMPTAGSDIR/chr$CHR -maftss $MAFTSSDDIR/chr$CHR -pval $PVALDIR/chr$CHR -ann $ANNOTDIR/chr$CHR -excl 895,975,976,977,978,979,980 -chr $CHR -o $F1 || { echo 'Failure!'; } 
done

echo 'Calculate effective number of annotations'
Rscript ${HPC_GARFIELD_BIN}/garfield-Meff-Padj.R -i $F1 -o $F0
NEA=$(head -1 $F0 |awk '{print $2}')
Padj=$(tail -1 $F0 |awk '{print $2}')

echo 'Calculate Enrichment and Significance'
F2=$OUTDIR/garfield.test.$INPUTNAME.out
Rscript ${HPC_GARFIELD_BIN}/garfield-test.R -i $F1 -o $F2 -l $ANNOTLINKFILE -pt $PTHRESH -b $BINNING -s $SUBSET -c $CONDITION
echo 'GARFIELD single annotation analysis complete'
echo $Padj 
echo 'Create Plots'
Rscript ${PLOTDIR}garfield-plot.R -i $F2 -o $F2 -l $ANNOTLINKFILE -t " " -f 10 -padj $Padj

echo 'GARFIELD Analysis Complete!'
