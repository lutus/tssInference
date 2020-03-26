#Title: geneSort
#Author: LDV
#Created: 2/20/2020
#Updated: 4/25/2020
#Description: runs through bpLocator function to find peaks
            #

peakFinder=function(agDf, bwPlus, bwMinus, tssWin, densityWin){
  ##--Arg--##
  #agDf: data frame of active Genes
  #bwPlus: bigWig file for the plus strand
  #bwMinus: bigWig file for the minus strand
  #tssWin: the number of bp to add TSS window to

  #result is a list of peaks for the given gene

  #prepares bed for bed6
  tmpDf=data.frame(agDf$chrom,agDf$start,agDf$end,
                  agDf$gene,agDf$active,agDf$strand,
                  stringsAsFactors=FALSE
  )
  tmpDf$agDf.start=tmpDf$agDf.start-tssWin - densityWin
  peaks=bed6.step.bpQuery.bigWig(bwPlus,bwMinus,tmpDf,step=1,op='sum')

  return(peaks)
}
