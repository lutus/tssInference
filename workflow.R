
#root dir
root=paste(getwd(),
            '/_dev/guertin/TSSinference/func/',
            sep=''
          )
# load test data
library(bigWig)
loc=paste(getwd(),'/_dev/guertin/TSSinference/_data/gencode.hg38.firstExon.bed.gz',sep='')
bed=read.table(gzfile(loc),
  col.names=c('chrom', 'start', 'end', 'gene', 'v5', 'strand'))
bwMinusLoc=paste(getwd(),'/_dev/guertin/TSSinference/_data/',
'H9_minus_PE2_merged.bigWig', sep='')
bwPlusLoc=paste(getwd(),'/_dev/guertin/TSSinference/_data/',
'H9_plus_PE2_merged.bigWig', sep='')
bwMinus=load.bigWig(bwMinusLoc)
bwPlus=load.bigWig(bwPlusLoc)

#source functions
source(paste(root,'versions/TSSinference_v0_6.R',sep=''))
#start time trial
t0=Sys.time()
agDf=activeGene(bwPlus,bwMinus,bed,tssWin=100)
agDfSec=agDf[1:4000,]

x=optimizer(agDf,bwPlus, bwMinus,tssWin=100,densityWin=100)
t1=Sys.time()
dt=difftime(t1,t0,units='mins')
dt
dim(x)
x
