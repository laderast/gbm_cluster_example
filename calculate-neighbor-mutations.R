library(BioNet)
source("permutation-script.R")


calcNeighborMutations <- function(intome, nodeset=NULL, gisticCopyCalls, 
                                  nodeMappedDream, cellLines=NULL, cores=NULL, 
				  prefix="", numPermutes=10000){
  #for permutation script, look at venus-permutation-script.R
  intNames <- nodes(intome)

  #build adjacency matrix
  intAdjmat <- as(intome, "graphAM")
  intAdjmat <- intAdjmat@adjMat
  rownames(intAdjmat) <- colnames(intAdjmat)
  
  #print(intAdjmat[1:10,1:10])
  
  #limit adjacency matrix to nodeset, if provided  
  if(!is.null(nodeset)){
    intAdjmat <- intAdjmat[,nodeset]
  }
  
  #mark nodes as being connected to themselves (necessary in the count)
  for(rp in colnames(intAdjmat)){ intAdjmat[rp,rp] <- 1}

  #cntMuts is function for apply to look for intersection in adjacency matrix and mutations
  cntMuts <- function(x){
    ed <- which(x>0)
    idx <- ed[ed %in% names(freqvec)]
    #print(idx)
    mutcount <- sum(freqvec[as.character(idx)])
    mutcount
  }

  cls <- intersect(colnames(gisticCopyCalls), as.character(as.character(unique(nodeMappedDream$CellLine))))
  
  if(!is.null(cellLines)){cellLines <- intersect(cls, cellLines)}
  else {cellLines <- cls}
  
  cellResults <- list()
  distMatrices <- list()
  mutCopyFrames <- list()

  for(cell in cellLines){
    print(cell)
  
  #cell <- "UACC812"
  
    cellCalls <- gisticCopyCalls[,cell]
    names(cellCalls) <- gisticCopyCalls$NodeName
    cellCalls <- cellCalls[cellCalls != 0]
    cellCopyCalls <- rep(1, length(cellCalls))
    names(cellCopyCalls) <- names(cellCalls)
  
    cellMutCalls <- table(as.character(nodeMappedDream[nodeMappedDream$CellLine == cell, "NodeName"]))
    cellMutCalls[na.omit(names(cellMutCalls))]
  
    #build a vector with all mutations and copy calls
    allMuts <- na.omit(union(names(cellCopyCalls), names(cellMutCalls)))
    allMutCopyCalls <- rep(0,length=length(allMuts))
    names(allMutCopyCalls) <- allMuts
    #allMutCopyCalls <- allMutCopyCalls[na.omit(names(allMutCopyCalls))]
  
    allCopyCalls <- allMutCopyCalls
    allCopyCalls[names(cellCopyCalls)] <- cellCopyCalls
    allCopyCalls <- allCopyCalls[na.omit(names(allCopyCalls))]
    allMutCalls <- allMutCopyCalls
    allMutCalls[names(cellMutCalls)]  <- cellMutCalls
    allMutCalls <- allMutCalls[na.omit(names(allMutCalls))]
  
    allMutCopyCalls <- allCopyCalls + allMutCalls
  
    mutCopyFrame <- data.frame(NodeName=names(allCopyCalls),copy=allCopyCalls, mutations=allMutCalls, combined=allMutCopyCalls)
    mutCopyFrame <- merge(mutCopyFrame, geneIntTable, by="NodeName")
  
    sampleSize <- sum(allMutCopyCalls)
  
    #allMutCopyCalls
    freqframe <- data.frame(allMutCopyCalls, NodeName=names(allMutCopyCalls))
  
    #need node indices
    freqind <- which(nodes(intome) %in% names(allMutCopyCalls))
    freqnames <- nodes(intome)[freqind]
    freqindframe <- data.frame(freqind, NodeName=freqnames)
  
    freqframe <- merge(freqframe, freqindframe, by="NodeName")
  
    freqvec <- freqframe$allMutCopyCalls
    names(freqvec) <- freqframe$freqind
    neighborVec <- apply(intAdjmat, 2, cntMuts)
    #print(neighborVec)
  
    if(!is.null(cores)){
      dist1res <- permuteNetsParallel(sampleSize, cellLine=cell, intAdjmat=intAdjmat, 
                                    numPermutes=numPermutes, cores=cores, prefix=prefix)
    }
    else{
      dist1res <- permuteNets(sampleSize, cellLine=cell, intAdjmat=intAdjmat, 
                              numPermutes=numPermutes, prefix=prefix)
    }
  
    colnames(dist1res) <- names(neighborVec)

    pvaldist1mutcopycount <- vector(length=length(neighborVec),mode="numeric")
    names(pvaldist1mutcopycount) <- names(neighborVec)
  
    for(nd in names(neighborVec)){
      distRP <- table(dist1res[,nd])
      obsvalue <- neighborVec[names(neighborVec) == nd]
      pvaldist1mutcopycount[nd] <- sum(distRP[as.numeric(names(distRP)) >=  obsvalue])/sum(distRP)
    }
  
    isMutated <- as.numeric(names(neighborVec) %in% names(allMutCopyCalls))
    degNode <- degree(intome)[names(neighborVec)]
  
    cellRes <- data.frame(NodeName = names(neighborVec), isMutated, degree=degNode, neighborVec, pvalue=pvaldist1mutcopycount[names(neighborVec)])
    cellRes <- merge(cellRes, geneIntTable, by="NodeName")
    write.table(cellRes, file=paste(prefix,cell, "-pvalue-results.txt", sep=""), sep="\t",quote=F, row.name=F)
    write.table(mutCopyFrame, paste(prefix,cell,"-copyNumberAndMutations.txt",sep=""), quote=F, sep="\t")
    
    distMatrices[[cell]] <- dist1res
    cellResults[[cell]] <- cellRes
    mutCopyFrames[[cell]] <- mutCopyFrame
  
 }
    outresults <- list(cellResults=cellResults, 
                       mutCopyFrames=mutCopyFrames)
    outresults
 
}

