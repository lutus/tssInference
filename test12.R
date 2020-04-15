
#root dir
root=paste(getwd(),
            '/_dev/guertin/TSSinference/func/',
            sep=''
          )
#source Data
source(paste(root,'tests/loadData.R',sep=''))

#source functions
source(paste(root,'versions/TSSinference_v0_6.R',sep=''))

t0=Sys.time()
agDf=activeGene(bwPlus,bwMinus,bed,tssWin=100)
#agDf[1,]
#agDf1=agDf[1:4000,]
#write.table(agDf1,'/home/lutus/_projects/_dev/guertin/TSSinference/docs/peaks.txt',sep='\t')
dim(agDf1)
tssWin=100
densityWin=100
x=optimizer(agDf1,bwPlus, bwMinus,tssWin=100,densityWin=100)

t1=Sys.time()
dt=difftime(t1,t0,units='mins')
dt

length(x)
sec1=unlist(x[1])
length(sec1)
sec1
agDf[3457,]
agDf$end[3457]-agDf$start[3457]
#end is not included in peaks
bpLoc=seq(agDf$start[3457]-tssWin-densityWin,agDf$end[3457]-1, by=1)
length(bpLoc)

peaksAnnotated=data.frame(bpLoc,sec1,stringsAsFactors=FALSE)
dim(peaksAnnotated[peaksAnnotated$seq1>0,])
peaksAnnotated
sec1[1]
