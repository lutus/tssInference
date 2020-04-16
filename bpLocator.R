#Title: optimizer
#Author: LDV
#Created: 2/2/2020
#Updated: 4/25/2020
#Description: uses to calculate the
            # up stream and down stream, bpHeight of each
            # peak

bpLocator=function(bw,chrom,start,end,strand,window){
  ##--Arg--##
  #bw: bigWig file of pileups
  #chrom, start, end, gene, score, strand
  #window: number of base pairs to look in either direction and compare densities

  #start at downstream extreme iterate to upstream
  bpLoc=as.integer(c())
  upDen=as.numeric(c())
  dwnDen=as.numeric(c())
  denDiff=as.numeric(c())
  bpHeight=as.numeric(c())
  j=1
  for(i in end:start) {
    bpLoc[j]=i
    dens=densityCalc(bw,chrom,i,strand,window)
    upDen[j]=dens[1]
    dwnDen[j]=dens[2]
    denDiff[j]=dens[2]-dens[1]
    bpHeight[j]=dens[3]
    #index
    j=j+1
  }
  denDiffDf=data.frame(bpLoc,
                upDen,
                dwnDen,
                denDiff,
                bpHeight,
                stringsAsFactors=FALSE)

  return(denDiffDf)

}
