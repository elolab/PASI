# Input: list of cleaned pathways
# Output: A list of vectors(pathways), including -1 or 1 for each node indicating if the node should
#         be active (1) or inactive (-1) when the pathway is active.

DetermineRoles = function(pathways, score){
  updatedpathways = lapply(pathways, function(x){
    
    # Extract relevant information
    relations = x$relationinfo
    nodes = x$nodeinfo
    id = nodes$Id
    if(!is.list(relations)){ # All nodes in pathways without edges are activating NOTE: is.list is fine as data frames return TRUE
      x$nodeinfo$Role = rep(1, length(id))
      return(x)
    }  
    ends = relations$EndId
    starts = relations$StartId
    relationdirs = as.numeric(c("activation"=1, "inhibition"=-1)[relations$Direction])
    # NOTE: relation direction means here inhibiting vs activiting relation, NOT defining which end is parent and which one is child
    
    # Initialise all node and relation roles as active (+1)
    noderoles = rep(1, length(id))
    relationroles = rep(1, length(ends))
    
    if(score == "activity_future"){ # Hopefully we will find a way to do this reliably for the next version of the package
      
      # Update node roles based on their majority leaving relation
      noderoles = UpdateNodeRoles(starts, ends, id, relationroles, noderoles, relationdirs)
      
      continue = TRUE
      iterationcounter = 0
      while(continue & (iterationcounter < 20)){
        
        nodeprev = noderoles
        relationprev = relationroles
        iterationcounter = iterationcounter + 1
        
        # Set relations inhibiting inhibitor nodes as activating
        relationroles = UpdateRelationRoles(ends, id, noderoles, relationdirs)
        
        # Update node roles
        noderoles = UpdateNodeRoles(starts, ends, id, relationroles, noderoles, relationdirs)
        
        # If no change from previous iteration, stop
        if(all(noderoles == nodeprev) & all(relationroles == relationprev)){ 
          continue=FALSE
        }
      }
      
      # If more than half of the nodes have role -1, the algorithm probably went "wrong way" and all the roles should be the opposite
      if(sum(noderoles) < 0){
        noderoles = -1 * noderoles
        relationroles = -1 * relationroles
        islandindex = which(!(id %in% starts) & !(id %in% ends))
        if(length(islandindex) > 0) noderoles[islandindex] = 1
      }
    }
    
    # Add node and relation roles to pathway info
    x$nodeinfo$Role = noderoles
    x$relationinfo$Role = relationroles
    
    return(x)
  })
  
  return(updatedpathways)
}

# Helper script for DetectRoles
# Updates node roles so that each node is the majority relation role leaving from the node

UpdateNodeRoles = function(starts, ends, nodeids, relationroles, noderoles, relationdirs){
  
  # Calculate number of child nodes for each node (named with node ids)
  childnumber = rep(0,length(nodeids))
  names(childnumber) = nodeids
  childnumber[nodeids[nodeids %in% starts]] = table(starts)[as.character(nodeids[nodeids %in% starts])]
  
  # Define indices of nodes with children, island nodes without parents or children, and leaf nodes without children but with parents
  index_parent = which(childnumber > 0)
  index_island = which((childnumber == 0) & !(nodeids %in% ends))
  index_youngest = setdiff(which(childnumber == 0), index_island)
  
  # Initialize all roles as 1
  noderoles = rep(1, length(nodeids))
  
  #index = which((nodeids %in% starts) | (nodeids %in% ends))
  #influence = relationroles * relationdirs
  #if(length(index) > 0){
  #  for(i in index){
  #    arrivingrelindex = which(ends == nodeids[i])
  #    leavingrelindex = which(starts == nodeids[i])
  #    arrivingeffect = sum(influence[arrivingrelindex])
  #    leavingeffect = sum(relationroles[leavingrelindex])
  #    if(arrivingeffect + leavingeffect <= 0) noderoles[i] = -1
  #  }
  #}
  
  # Leaf nodes' roles are defined by majority vote from arriving relations (role*direction)
  parentalinfluence = relationroles * relationdirs
  if(length(index_youngest) > 0){
    for(i in index_youngest){
      arrivingrelindex = which(ends %in% nodeids[i])
      majorityrole = 1
      if(sum(parentalinfluence[arrivingrelindex]) <= 0) noderoles[i] = -1
    }
  }
  
  # Nodes with children are updated according to majority vote of leaving relation roles
  if(length(index_parent) > 0){
    for(j in index_parent){
      leavingrelindex = which(starts %in% nodeids[j])
      leavingrelroles = relationroles[leavingrelindex]
      if(sum(leavingrelroles) <= 0) noderoles[j] = -1
      #if(sum(leavingrelroles*(childnumber[ends[leavingrelindex]]+1)) <= 0) noderoles[j] = -1 
    }
  }

  return(noderoles)
}


#
UpdateRelationRoles = function(ends, nodeids, noderoles, relationdirs){
  
  # Define which relations activate inhibitor nodes or inhibit activator nodes (these relations' roles are -1)
  endnodeindices = match(ends, nodeids, nomatch=0)
  index_inhibitorsactivator = which((noderoles[endnodeindices] == -1) & (relationdirs == 1))
  index_activatorsinhibitor = which((noderoles[endnodeindices] == 1) & (relationdirs == -1))
  
  # Set new relation roles
  relationroles = rep(1, length(ends))
  relationroles[c(index_inhibitorsactivator, index_activatorsinhibitor)] = -1
  
  return(relationroles)
}