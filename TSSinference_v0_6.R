#root dir
root=paste(getwd(),
            '/_dev/guertin/TSSinference/func/',
            sep=''
          )
#source functions
#density
source(paste(root,'densityCalc/densityCalc_v0_1.R',sep=''))
#activeGene
source(paste(root,'activeGene/activeGene_v0_5.R',sep=''))
#bpLocator
source(paste(root,'bpLocator/bpLocator_v0_2.R',sep=''))
#plots
source(paste(root,'plots/plotProfile_v0_1.R',sep=''))
#geneSort
source(paste(root,'geneSort/geneSort_v0_0.R',sep=''))
#peakFinder
source(paste(root,'peakFinder/peakFinder_v0_3.R',sep=''))
#optimizer
source(paste(root,'optimizer/optimizer_v0_1.R',sep=''))
#decisiontree
source(paste(root,'decisionTree/decisionTree_v0_1.R',sep=''))
