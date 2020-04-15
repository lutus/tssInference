#Title: decisionTree
#Author: LDV
#Created: 3/25/2020
#Updated: 4/4/2020
#Description: determines the peak for a data frame of peaks
            #Tests:
            #     Peak profile
            #     Upstream/downstream density
            #      next nearest neighbor

decisionTree=function(agDf,agDfPeaks,bwPlus,bwMinus,tssWin,densityWin){
  ##--Arg--##
  #input:
        # agDf: data frame  of genes and coresponding info
        # agDfPeaks: list of peak heights for each gene in agDf
  #output:
        # agDf: original agDf with column for peakLoc

  ## loop through agDf
  for(i in 1:nrow(agDf)){
    #-------------Prep data row-------------#
    #gene window start
    gWinStart=agDf$start[i]
    #gene window end
    gWinEnd=agDf$start[i]

    #peaks in gene
    bpHeight=unlist(peaks[i])
    #create vector of bp loc
    bpLoc=seq(agDf$start[i]-tssWin-densityWin,agDf$end[i]-1, by=1)
    #align bpLoc to bpHeight in df
    bpPeaks=data.frame(bpLoc,bpHeight,stringsAsFactors=FALSE)
    bpPeaks=bpPeaks[order(-bpPeaks$bpHeight),]
    bpPeaks=bpPeaks[bpPeaks$bpHeight!=0,]
    top20=head(bpPeaks,20)
    #-------------Upstream Downstream density check-------------#
    #calculate upstream/downstream densities
    if(agDf[i,]$strand=='+'){
      #loop through top20
      for(j in 1:nrow(top20)){
        x=densityCalc(bwPlus,agDf$chrom[i],top20$bpLoc[j],agDf$strand[i],densityWin)
        top20$upDen[j]=as.numeric(x[1])
        top20$dwnDen[j]=as.numeric(x[2])
        loc=top20$bpLoc[i]
        tmp=top20
        tmp$nnn=abs(tmp$bpLoc-loc)
        tmp=tmp[order(tmp$nnn),]
        nnn=tmp$nnn[1]
      }
    }else if(agDf[i,]$strand=='-'){
      #loop through top20
      for(j in 1:nrow(top20)){
        x=densityCalc(bwMinus,agDf$chrom[i],top20$bpLoc[j],agDf$strand[i],densityWin)
        top20$upDen[j]=as.numeric(x[1])
        top20$dwnDen[j]=as.numeric(x[2])
        loc=top20$bpLoc[i]
        tmp=top20
        tmp$nnn=abs(tmp$bpLoc-loc)
        nnn=tmp$nnn[1]
      }
    }


    #determine peak profile
    #peakCnt=dim(bpPeaks[bpPeaks$bpHeight>0,])[1]
    #avgPeakHeight=average()
  }## close loop through agDf
  #-------------Set peak to agDf-------------#

  ## return agDf
  return(agDf)
}## decisionTree close function
