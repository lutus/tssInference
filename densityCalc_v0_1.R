#Title: densityCalc
#Author: LDV
#Created: 2/3/2020
#Updated: 4/7/2020
#Description: support function for other parts of the code
            # calculates densties upstream and downstream


densityCalc = function(bw, chrom, bpLoc, strand, window){
  ##--Arg--##
  #bw: bigWig file of pileups
  #chrom, start, end, gene, score, strand
  #window: number of base pairs to look in either direction
  #and compare densities. Also include bp height

  upWinStart=bpLoc-window
  upWinEnd=bpLoc-1
  dwnWinStart=bpLoc+1
  dwnWinEnd=dwnWinStart+window
  #upstream density
  upDen=region.bpQuery.bigWig(bw,chrom,upWinStart,upWinEnd, op='avg')
  dwnDen=region.bpQuery.bigWig(bw,chrom,dwnWinStart,dwnWinEnd, op='avg')
  bpHeight=region.bpQuery.bigWig(bw,chrom,bpLoc,bpLoc+1,op='sum')
  return (c(upDen, dwnDen,bpHeight))
}
