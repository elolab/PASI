# Input: list of cleaned pathways
# Output: A list of vectors(pathways), including -1 or 1 for each node indicating if the node should
#         be active (1) or inactive (-1) when the pathway is active.

DetermineRoles = function(pathways, score){
  updatedpathways = lapply(pathways, function(x){
    
    # Extract relevant information
    relations = x$relationinfo
    nodes = x$nodeinfo
    id = nodes$nodeid
    if(!is.list(relations)){ # All nodes in pathways without edges are activating
      x$nodeinfo$noderole = rep(1, length(id))
      return(x)
    }  
    ends = relations$endnode
    starts = relations$startnode
    directions = relations$direction
    
    # Initialize node roles within pathway as active and relations according to their directions
    noderoles = rep(1, length(id))
    relationroles = c("activation"=1, "inhibition"=-1)
    relationroles = relationroles[directions]
    
    if(score == "activity"){
      
      # Update node roles based on their majority leaving relation
      noderoles = UpdateNodeRoles(starts, id, relationroles, noderoles)
      
      continue = TRUE
      iterationcounter = 0
      while(continue & (iterationcounter < 20)){
        
        nodeprev = noderoles
        relationprev = relationroles
        iterationcounter = iterationcounter + 1
        
        # Set relations inhibiting inhibitor nodes as activating
        inhinhrelindex = which(directions == "inhibition" & noderoles[match(ends, id, 0)] == -1)
        if(length(inhinhrelindex) > 0){
          relationroles[inhinhrelindex] = 1
        }
        
        # Update node roles based on their majority leaving relation
        noderoles = UpdateNodeRoles(starts, id, relationroles, noderoles)
        
        # Set relations activating inhibitor nodes as inhibiting
        actinhrelindex = which(directions == "activation" & noderoles[match(ends, id, 0)] == -1)
        if(length(inhinhrelindex) > 0){
          relationroles[actinhrelindex] = -1
        }
        
        # Update node roles based on their majority leaving relation
        noderoles = UpdateNodeRoles(starts, id, relationroles, noderoles)
        
        # If no change from previous iteration, stop
        if(all(noderoles == nodeprev) & all(relationroles == relationprev)){ 
          continue=FALSE
        }
      }
    }
    
    
    
    # Add node and relation roles to pathway info
    x$nodeinfo$noderole = noderoles
    x$relationinfo$relationrole = relationroles
    
    return(x)
  })
  
  return(updatedpathways)
}

# Helper script for DetectRoles
# Updates node roles so that each node is the majority relation role leaving from the node

UpdateNodeRoles = function(starts, nodeids, relationroles, noderoles){
  for(i in 1:length(noderoles)){
    nodeid = nodeids[i]
    leavingrelindex = which(starts %in% nodeid)
    if(length(leavingrelindex) > 0){
      leavingrelroles = relationroles[leavingrelindex]
      majorityrole = 1
      if(sum(leavingrelroles) < 0){
        majorityrole = -1
      }
      noderoles[i] = majorityrole
    }
  }
  return(noderoles)
}


