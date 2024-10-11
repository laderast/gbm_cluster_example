source("calculate-neighbor-mutations.R")
load("pancan.RData")

res1 <- calcNeighborMutations(intome, nodeset=tcgaNodes, focalBRCAnodeMapped, 
                              BRCAmutsNodeMapped, cellLines=patientsComplete[1:3], 
                              cores=20, prefix="pancanBRCA-", numPermutes=10000)

save.image("pancan-results.RData")