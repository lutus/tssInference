#Title: optimizer
#Author: LDV
#Created: 4/20/2020
#Updated: 4/25/2020
#Description: meassures size gene data frame and chunks out parts
            # to ensure not over stressing machine
            # ideally not to exceed 2 gb
optimizer=function(agDf,bwPlus, bwMinus,tssWin, densityWin){
  ##--Arg--##
  #agDf: data frame of active Genes
  #bwPlus: bigWig file for the plus strand
  #bwMinus: bigWig file for the minus strand
  #tssWin: the number of bp to add TSS window to

  #tresult is a data frame of gene and possible start site

  #-------------calculate break points-------------#

  #95,796,160 bp gives an object of roughly 2GB
  #sections is a vector of break points it is the upper end range
  # ie use ranges 1:sections[1], sections[1]+1:sections[2],
  # sections[n]+1:sections[n+1], ... , sections[n]+1:sections[nrow(agDf)]
  agDf$numBP=agDf$end-agDf$start
  sections=c(1)
  bpCount=0
  st=1
  for(i in st:nrow(agDf)){
    bpCount=bpCount+agDf$numBP[i]
    if(bpCount<95796160){next}
    else if(bpCount>=95796160){
      sections=append(sections,i-1)
      bpCount=0
      next
    }
  }
  #create column for tssPeak
  agDf$tssPeak=NA
  #-------------loop through sections-------------#
  #create result df:
  ##  gene, bpLoc,
  resultGene=c()
  resuiltBpLoc=c()
  resultDf=data.frame(resultGene,resuiltBpLoc,
                      stringsAsFactors = FALSE)

  for(i in 1:length(sections)){
    ##-------------determine start stop of section-------------##
    #conditional for first set
    if(i==1){
      agDfIndex=sections[i]
      agDfStop=sections[i+1]
    }
    #condition for the rest
    else if (i>1){
      agDfIndex=sections[i]+1
      if(i<length(sections)){
        agDfStop=sections[i+1]
      }
      else if(i==length(sections)){
        agDfStop=nrow(agDf)
      }
    }
    ##-------------run analysis on section-------------##
    # cat(paste('section ',i, '\n','start: ',agDfIndex,
    #  '\nstop: ',agDfStop, '\n\n', sep=''))

    if(abs(agDfIndex-agDfStop)>0){
      agDfSec=agDf[agDfIndex:agDfStop,]
    }else if(abs(agDfIndex-agDfStop)==0){
      agDfSec=agDf[agDfIndex,]
    }
    cat('input Dim: ', dim(agDfSec)[1],'\n',sep='')
    # cat(as.character(agDfSec$gene[1]),'\n',
    #     agDfSec$chrom[1],'\n',
    #     agDfSec$start[1],'\n',
    #     agDfSec$end[1],'\n',
    #     sep='')
    cat('section: ', i,'\nstart:', agDfIndex,
         '\nend: ',agDfStop,'\n',sep='')
    #clean data frame with geneSort
    agDfSec=geneSort(agDfSec)
    #find peaks
    if(dim(agDfSec)[1]==0){
      cat('geneSort removed from possible calculation',
      '\ndone\n','#----------------#\n',sep='')
      next
    }

    agDfSecPeaks=peakFinder(agDfSec,bwPlus, bwMinus, tssWin, densityWin)
     cat('\ndone\n','#----------------#\n',sep='')
    agDfSec=decisionTree(agDfSec,agDfSecPeaks,bwPlus,bwMinus,tssWin,densityWin)
##-------------append section data to ultimate df-------------##
    for(j in 1:nrow(agDfSec)){
      agDf[agDf$gene==agDfSec$gene[j],]$tssPeak=agDfSec$tssPeak[j]

    }

  }

  #result is a data frame

  return(agDf)
}
