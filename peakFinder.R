#Title: peakFinder
#Author: LDV
#Created: 2/20/2020
#Updated: 4/4/2020
#Description: loops through agDf to find peaks
            #

peakFinder=function(agDf, bwPlus, bwMinus, tssWin, densityWin){
  ##--Arg--##
  #input:
        # agDf: data frame of active Genes
        # bwPlus: bigWig file for the plus strand
        # bwMinus: bigWig file for the minus strand
        # tssWin: the number of bp to add TSS window to
  #ouput:
    # peaks: list of peaks for each gene in agDf
    
  #prepares bed for bed6
  tmpDf=data.frame(as.character(agDf$chrom),agDf$start,agDf$end,
                  as.character(agDf$gene),as.character(agDf$active),as.character(agDf$strand),
                  stringsAsFactors = FALSE
  )
  if (is.null(tmpDf) || dim(tmpDf)[1] == 0) {
    cat('Null: ',is.null(tmpDf),'\ndim: ', dim(tmpDf)[1],
    sep='')
  }
  tmpDf$agDf.start=tmpDf$agDf.start-tssWin - densityWin
  peaks=bed6.step.bpQuery.bigWig(bwPlus, bwMinus, tmpDf, step=1,
    op='sum', with.attributes=TRUE)

  return(peaks)
}
