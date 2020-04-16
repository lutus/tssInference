#Title: activeGene
#Author: LDV
#Created: 2/3/2020
#Updated: 4/31/2020
#Description: takes a bed file of 1st exons and creates a
            # new bed file of 1st exons ov genes that are
            # active and returns a data frame with
            # gene,chrom,start,end,strand,
            #  active,plusCnt, minusCnt, region

activeGene = function(bwPlus, bwMinus, bed, tssWin){
  ##--Arg--##
  #bed: bed as a data frame of 1st exons in a gene
  #bwPlus: bigWig file for the plus strand
  #bwMinus: bigWig file for the minus strand
  #tssWin: integer number of bp to look at beforegene start
  cat('Debugging',sep='')
  # list unique genes in bed file
  gene=unique(bed$gene)
  # initiate empty vectors
  chrom=as.character(c())
  start=as.integer(c())# start position of gene
  end=as.integer(c())# end position of gene
  strand=as.character(c())# strand to use for preak calculation

  plusCnt=as.integer(c())# number of exons on plus strand in gene
  minusCnt=as.integer(c())#  number of exons on minus strand in gene
    # calculated using bigWig functions
  active=as.logical(c())# logical whether gene is active
  region=as.integer(c())# integer of pileups in region

  warning=as.integer(c())# vector of integers describing warning
  warningFlag=as.logical(c())# logical flag for sorting

  # loop through unique genes
  for(i in gene){
    #reset gene warning vector
    gWarning=as.integer(c())

    # create bed of gene isoforms
    gBed=bed[bed$gene==i,]

    ################
    # determine chrom
    gChrom=as.character(unique(gBed$chrom))
    if(length(gChrom)>1){
      chrom=append(chrom,'chrX')# append to chrom vector
      gWarning=append(gWarning,1)
    }
    else if(length(gChrom)==1){
      chrom=append(chrom,as.character(gChrom))
    }

    ################
    # determine strandness
    pC=dim(gBed[gBed$gene==i&gBed$strand=='+',])[1]
    mC=dim(gBed[gBed$gene==i&gBed$strand=='-',])[1]
    if(pC>mC){
      gStrand='+'
    }
    else if(mC>pC){
      gStrand='-'
    }
    #-----------------------#
    else if(mC==pC){
      gbedPlus=gBed[gBed$strand=='+']
      #num of bp
      gbedPlusbp=max(gBedPlus$start)-min(gBedPlus$start)
      gbedMinus=gBed[gBed$strand=='-']
      gbedMinusbp=max(gBedMinus$start)-min(gBedMinus$start)
      if(gbedPlusbp>gbedMinusbp){
        gStrand='+'
      }
      else if(gbedPlusbp<gbedMinusbp){
        gStrand='-'
      }
      else if(gbedPlusbp==gbedMinusbp){
        gStrand=NA
      }
    }# need to test
    strand=append(strand,gStrand)
    #-----------------------#

    #add to vectors
    plusCnt=append(plusCnt,pC)
    minusCnt=append(minusCnt,mC)
    #if(gStrand==NA){gWarnings=append(gWarnings,2}
    ################
    # determine start, end
    gBedStarts=c()
    for(row in 1:nrow(gBed)){
      if(as.character(gBed[row,'strand'])=='+'){
        gBedStarts=append(gBedStarts,as.integer(gBed[row,'start']))
      }
      else if(as.character(gBed[row,'strand'])=='-'){
        gBedStarts=append(gBedStarts,as.integer(gBed[row,'end']))
      }
    }
    # append gene start and end
    gStart=min(gBedStarts)
    start=append(start,min(gBedStarts))
    gEnd=max(gBedStarts)
    end=append(end,max(gBedStarts))

    if(length(gChrom)>1){gChrom=gChrom[1]}
    #use start, end to determine active
    if(gStart!=gEnd){
      if(gStrand=='+'){
        gRegion=region.bpQuery.bigWig(bwPlus,start=gStart-tssWin,
                                  end=gEnd,
                                  chrom=gChrom)
                                }
      if(gStrand=='-'){
        gRegion=region.bpQuery.bigWig(bwMinus,start=gStart-tssWin,
                                    end=gEnd,
                                    chrom=gChrom)
                                }}
    else if(gStart==gEnd){gRegion=0}
    if(gRegion>0){
        active=append(active,TRUE)
        region=append(region, gRegion)
    }
    else if(gRegion==0){
        active=append(active,FALSE)
        region=append(region, gRegion)
    }

  }
  # if(length(gWarnings)>0){
  #   warnings=append(warningFlag,TRUE)
  # }
  # else if(length(gWarnings)==0){
  #   warnings=append(warningFlag,FALSE)
  # }
  return(data.frame(gene,chrom,start,end,strand,
    active,plusCnt, minusCnt, region,stringsAsFactors = FALSE))
}
