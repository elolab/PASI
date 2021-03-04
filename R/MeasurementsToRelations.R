# Input: "pathway" is a list including pathway information, 
#        "zscores" is a matrix including z-scores of expression data (rows: genes, cols: samples) 
# Output: Returns a list of matrices: [pathway1, ..., pathwayN], 
#         where pathway1 = matrix (rows: relations, columns: samples), missing values are NA 

MeasurementsToRelations = function(pathways, zscores){
  
  relationvalues = lapply(pathways, function(p){

    if(is.list(p$relationinfo)){
      node = p$nodeinfo
      relation = p$relationinfo
      roles = relation$relationrole
      
      # Initialize relation value matrix for pathway p
      relationnumber = length(relation$startnode)
      prelationvalues = mat.or.vec(relationnumber, ncol(zscores))
      rownames(prelationvalues) = paste("Relation", 1:relationnumber, sep="")
      colnames(prelationvalues) = colnames(zscores)
      
      # Fill the value matrix by rows (relations)
      for(i in 1:relationnumber){
        
        # Detect which Entrez ids are mapped to the parent and child node of the relation
        startnodeid = relation$startnode[i]
        endnodeid = relation$endnode[i]
        interaction = relation$direction[i]
        startentrez = node$nodeentrez[[which(node$nodeid %in% startnodeid)]]
        endentrez = node$nodeentrez[[which(node$nodeid %in% endnodeid)]]
        
        # Detect the scaled value of the parent node of the relation
        startindex = which(rownames(zscores) %in% startentrez)
        if(length(startindex)!=0){
          startvalues = colMeans(matrix(zscores[startindex,], nrow=length(startindex)), na.rm=T)
          startvalues[is.nan(startvalues)] = NA
          
        } else{
          startvalues = rep(NA, ncol(zscores))
        }
        
        # Detect the scaled value of child node of the relation
        endindex = which(rownames(zscores) %in% endentrez)
        if(length(endindex)!=0){
          endvalues = colMeans(matrix(zscores[endindex,], nrow=length(endindex)), na.rm=T)
          endvalues[is.nan(endvalues)] = NA
          
        } else{
          endvalues = rep(NA, ncol(zscores))
        }
        
        # Calculate scaling parameter (high with active parent)
        scalingparam = (0.5*startvalues+1.5)/dnorm(0,0,0.7)
        
        # Calculate relation values (NA for samples with missing parent or child measurement)
        # Note: this could be further developed by considering also other parent nodes of the child node
        if(interaction == "activation"){ prelationvalue = scalingparam*dnorm(endvalues, mean=startvalues, sd=0.7)-1
        } else prelationvalue = scalingparam*dnorm(-endvalues, mean=startvalues, sd=0.7)-1
        
        # Save the values for the relation
        prelationvalues[i,] = roles[i]*prelationvalue
      }  
      
      return(prelationvalues)
      
    } else return(NA)
    
  })
  
  return(relationvalues)
}