#Title: geneSort
#Author: LDV
#Created: 2/20/2020
#Updated: 4/25/2020
#Description: filters the result of activeGene function
            # removes genes that cause errors
geneSort=function(activeGeneDF){
  ##--Arg--##
  #activeGeneDF: data frame prpoduced by activeGene function

  #returns a cleaned data frame of genes

  #remove those with no counts
  activeGeneDF=activeGeneDF[activeGeneDF$active==TRUE,]

  #remove strand errors
  activeGeneDF=activeGeneDF[(activeGeneDF$strand=='+'|activeGeneDF$strand=='-'),]

  #remove genes over 2.5M bp
  activeGeneDF=activeGeneDF[(activeGeneDF$end-activeGeneDF$start)<2500000,]

  return(activeGeneDF)
}#end function
#[c('chrom','start','end',
              #        'gene','active','strand')]
