# Input: "pathway" is a list including pathway information, 
#        "data" is a scaled expression data frame (rows: genes, cols: samples)
# Output: Returns a list of matrices: [pathway1, ..., pathwayN], 
#         where pathway1 = matrix (rows: nodes, columns: samples), missing values are NA 

MeasurementsToNodes = function(pathways, data){
  
  samplenumber = ncol(data)
  nodevalues = lapply(pathways, function(p){
    node = p$nodeinfo
    nodedirections = node$Role
    
    # Initialize node value matrix for pathway p
    nodenumber = nrow(node)
    pnodevalues = matrix(NA, nrow=nodenumber, ncol=samplenumber)
    rownames(pnodevalues) = node$Id
    colnames(pnodevalues) = colnames(data)
    
    # Fill the value matrix by rows (nodes)
    for(i in 1:nodenumber){
      entrez = node$Entrez[i]
      index = match(entrez, rownames(data), nomatch=0)
      
      if(sum(index)!=0){
        index = index[index>0]
        
        # A node is type gene -> one Entrez mapped to the node
        if(length(index) == 1){ values = data[index,]
        # A node is type group -> multiple Entrez mapped to the node, node value is the mean
        } else{
          values = colMeans(data[index,], na.rm=T)
          values[is.nan(values)] = NA
        }
        
        pnodevalues[i,] = nodedirections[i]*values
        
      }
    }
    
    return(pnodevalues)
  })
  names(nodevalues) = names(pathways)
  
  return(nodevalues)
}