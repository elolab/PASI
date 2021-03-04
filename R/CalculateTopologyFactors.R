# Input: argument "pathway"  is a list including pathway information, "neighbours" and "secondneighbours" are lists including
#        argument "statistics" is a list including information about neighbors, second neighbors and occurrences 
# Output: Returns a list of importance values (0.9-1.1), one value per node in each "pathway".

CalculateTopologyFactors = function(pathways, statistics){
  
  importances = list()
  pathwayoccurrences = statistics$occurrences
  pathwayfneighbors = statistics$fneighbors
  pathwaysneighbors = statistics$sneighbors
  for(p in 1:length(pathways)){
    
    # Pick initial info
    pathway = pathways[[p]]
    occurrences = pathwayoccurrences[[p]]
    fneighbors = pathwayfneighbors[[p]]
    sneighbors = pathwaysneighbors[[p]]
    nodes = pathway$nodeinfo$nodeid
    
    # Detect neighborvalues and scale them to interval [0,1] within pathway
    neighboreffect = fneighbors + 0.5*sneighbors
    greatest = max(neighboreffect, na.rm=T)
    if(greatest != 0) neighboreffect = neighboreffect/greatest
    
    # Scale occurrences to interval [0,1] within pathway
    scaledoccurrences = occurrences/max(occurrences) 
    
    # Calculate and scale importances to interval (0.9, 1.1]
    rawimportances = neighboreffect/scaledoccurrences
    block = 0.2/length(fneighbors)
    scaledimportances = 0.9 + rank(rawimportances)*block
    
    importances[[p]] = scaledimportances
  }
 return(importances) 
}
