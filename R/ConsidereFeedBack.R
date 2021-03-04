# Input: numeric vector "values" include a value for each node (-1 if the node isn't measured), 
#        list "pathway" includes pathway information, vector "occurances" tells how many pathways include each node in 
#        "pathway" and logical vector "measured" tells if a node is measured or not.
# Output: Returns a numeric vector. The values are node values related to pathway "pathway" when feed back   
#         phenomena is taken account.

ConsidereFeedBack = function(nodevalues, pathways, alloccurances){
  
  for(i in 1:length(nodevalues)){ # i is pathway

    # Extract pathway info
    pathway = pathways[[i]]
    
    if(is.list(pathway$relationinfo)){
      
      allend = pathway$relationinfo$endnode
      allstart = pathway$relationinfo$startnode
      
      for(j in 1:ncol(nodevalues[[1]])){ # j is sample
     
        values = nodevalues[[i]][,j]
        measured = which(!is.na(values))
        originalvalues = values
        nodes = pathway$nodeinfo$nodeid[measured]
        values = values[measured]
        occurances = alloccurances[[i]][measured]
        
        keeprelationindex = which((allend %in% nodes) & (allstart %in% nodes))
        end = allend[keeprelationindex]
        start = allstart[keeprelationindex]
        
        structure = FindSCC(start, end, nodes)
        scc = unlist(lapply(structure, length))
        
        if(all(scc == 1)){ # The pathway is acyclic
          inferredvalues = WithoutCycles(values, nodes, start, end, occurances, measured)   
          
        } else{ # The pathway includes cycles
          inferredvalues = WithCycles(structure, values, nodes, start, end, occurances, measured)    
        }
        
        originalvalues[measured] = inferredvalues
        nodevalues[[i]][,j] = originalvalues
      }
    }
  }
  
  return(nodevalues)
}


# Helper script for ConsidereFeedBack
# Returns a list of strongly connected components

FindSCC = function(start, end, nodes){
  scc = list()
  
  #Preprocessing: detect unconnected nodes and remove them
  islands = setdiff(nodes, c(start,end)) 
  if(length(islands)!=0){
    scc = as.list(islands)
    nodes = setdiff(nodes, islands) 
  }
  
  #Preprocessing: detect root and leaf nodes and remove them iteratively
  roots = setdiff(start, end) 
  leafs = setdiff(end, start)
  rootsandleafs = unique(c(roots,leafs))
  while(length(rootsandleafs)!=0){
    scc = as.list(c(rootsandleafs, unlist(scc)))
    nodes = setdiff(nodes, rootsandleafs)
    relationdropindex = unique(which((end %in% rootsandleafs) | (start %in% rootsandleafs)))
    start = start[-relationdropindex] 
    end = end[-relationdropindex]
    roots = setdiff(start, end)
    leafs = setdiff(end, start)
    rootsandleafs = unique(c(roots,leafs))
  }
  
  #Preprocessing: detect unconnected nodes after removing root and leaf nodes
  islands = setdiff(nodes, c(start,end)) 
  if(length(islands)!=0){
    scc = as.list(c(islands, unlist(scc)))
    nodes = setdiff(nodes, islands) 
  }
  
  #Actual algorithm
  index = 1 
  nodestoanalyze = nodes
  while(length(nodestoanalyze)!=0){ 
    dfnumbers = DF(start, end, nodestoanalyze, nodestoanalyze[1], index)
    visitedindex = which(dfnumbers!=0)
    visitednodes = nodestoanalyze[visitedindex]
    ancestornumbers = FindOldestAncestor(dfnumbers[visitedindex], visitednodes, start, end)
    differentancestors = unique(ancestornumbers)
    for(j in 1:length(differentancestors)){
      ancestor = differentancestors[j]
      component = visitednodes[which(ancestornumbers == ancestor)]
      scc[[length(scc)+1]] = component
    }
    
    nodestoanalyze = nodestoanalyze[-visitedindex]
    relationdropindex = which((end %in% visitednodes)|(start %in% visitednodes))
    end = end[-relationdropindex]
    start = start[-relationdropindex]
    index = index+length(visitednodes)
  }
  
  return(scc)
}


# Helper script for FindSCC
# Returns an integer vector indicating the order nodes are visited from node "from" (0 if not visited at all)  

DF = function(start, end, nodes, from, indexfrom){
  history = from
  visited = rep(0, length(nodes))
  index = indexfrom
  while(length(history)!=0){
    current = history[1]
    history = history[-1]
    currentnodeindex = which(nodes %in% current)
    
    if(visited[currentnodeindex]==0){
      visited[currentnodeindex] = index
      index = index+1
      currentstartindex = which(start %in% current)
      nextnodes =  unique(end[currentstartindex])
      history = c(nextnodes, history)
    }
  }
  return(visited)
}


# Helper script for FindSCC
# Returns an integer vector indicating the root most parent id of each node

FindOldestAncestor = function(dfnumbers, nodes, start, end){ #tarkasta
  children = lapply(nodes, function(x){
    relationindex = which(start %in% x)
    kids = end[relationindex]
    return(kids)
  })
  
  reachable = as.list(nodes)
  change = T
  while(change){
    change = F
    for(i in 1:length(nodes)){
      self = reachable[[i]]
      index = which(nodes %in% self)
      kids = unlist(children[index])
      newreachable = unique(c(self,kids))
      if(length(reachable[[i]])!=length(newreachable)){
        reachable[[i]] = newreachable
        change = T
      }
    }
  }
  
  oldest = unlist(lapply(reachable, function(y){
    nodeindex = which(nodes %in% y)
    reachabledf = dfnumbers[nodeindex]
    minimum = min(reachabledf)
    return(minimum)
  }))
  return(oldest)
}


# Helper script for ConsidereFeedBack
# Returns vector of feed back modified values in case the pathway is acyclic

WithoutCycles = function(values, nodes, start, end, occurances, measured){
  inferred = rep(F,length(occurances))
  inferredvalues = values
  
  currentnodeindex = which(nodes %in% setdiff(end, start)) # Initially leaf nodes
  inferred[currentnodeindex] = T
  islandindex = which(nodes %in% setdiff(nodes, c(start,end)))
  inferred[islandindex] = T
  
  while(length(currentnodeindex)!=0){
    
    addnode = NULL
    dropindex = NULL
    
    parentnodes = lapply(currentnodeindex, function(x){
      parentrelationindex = which(end %in% nodes[x])
      return(start[parentrelationindex])
    })
    
    for(i in 1:length(parentnodes)){ 
      parents = parentnodes[[i]]
      votefori = T # Should current node i be dropped 
      
      # Go through parents if i isn't a root node
      if(length(parents)!=0){ 
        for(j in parents){
          
          childrelationindex = which(start %in% j)
          kids = end[childrelationindex]
          childnodeindex = which(nodes %in% kids)
          parentnodeindex = which(nodes %in% j)
          
          if(!inferred[parentnodeindex]){
            if(all(inferred[childnodeindex])){
              childvalues = inferredvalues[childnodeindex]
              if(median(childvalues, na.rm=T) >= 0.8*values[parentnodeindex]){
                value = weighted.mean(childvalues, 1/occurances[childnodeindex])               
                inferredvalues[parentnodeindex] = max(value,values[parentnodeindex],na.rm=T)
              }              
              inferred[parentnodeindex] = T
              
              addnode = c(addnode,j)

            } else votefori = F
          }        
        }
      }
      if(votefori){
        dropindex = c(dropindex,i)
      }
    }
    
    if(length(dropindex) != 0) currentnodeindex = currentnodeindex[-dropindex]
    if(length(addnode) != 0) currentnodeindex = c(currentnodeindex, which(nodes %in% addnode))             
  }
  return(inferredvalues)
}


# Helper script for ConsidereFeedBack
# Returns vector of feed back modified values in case the pathway includes cycles

WithCycles = function(structure, values, nodes, start, end, occurances, measured){
  inferredvalues = values
  inferred = rep(F,length(structure))
  bigstructure = CycleLevelRelations(structure, start, end)
  start = bigstructure$startrelation # Unlike with nodes, these are cycle indexes!
  end = bigstructure$endrelation
  
  currentcycleindex = setdiff(end, start) # Initially leaf cycles
  inferred[currentcycleindex] = T
  islandindex = setdiff(1:length(structure), c(start,end))
  inferred[islandindex] = T
  
  while(length(currentcycleindex)!=0){
    
    addcycle = NULL
    dropindex = NULL
    
    parentcycles = lapply(currentcycleindex, function(x){
      parentrelationindex = which(end %in% x)
      return(start[parentrelationindex])
    })
    
    for(i in 1:length(parentcycles)){ 
      parents = parentcycles[[i]]
      votefori = T # Should current cycle i be dropped 
      
      # Go through parents if i isn't a root cycle
      if(length(parents)!=0){ 
        for(j in parents){
          
          childrelationindex = which(start %in% j)
          childcycleindex = end[childrelationindex] 
          childnodeindex = which(nodes %in% unlist(structure[childcycleindex]))
          parentnodeindex = which(nodes %in% structure[[j]])
          
          if(!inferred[j]){
            if(all(inferred[childcycleindex])){
              childvalues = inferredvalues[childnodeindex]
              if(median(childvalues, na.rm=T) >= 0.8*median(values[parentnodeindex], na.rm=T)){
                value = weighted.mean(childvalues, 1/occurances[childnodeindex])               
                for(k in parentnodeindex){
                  inferredvalues[k] = max(value,values[k],na.rm=T)  
                }
              }
              
              inferred[j] = T
              
              addcycle = c(addcycle,j)

            } else votefori = F
          }        
        }
      }
      if(votefori){
        dropindex = c(dropindex,i)
      }
    }
    
    if(length(dropindex) != 0) currentcycleindex = currentcycleindex[-dropindex]
    if(length(addcycle) != 0) currentcycleindex = c(currentcycleindex, addcycle)             
  }
  
  return(inferredvalues)
}

# Helper script for WithCycles
# Returns relations between cycles

CycleLevelRelations = function(structure, start, end){
  startcycle = NULL
  endcycle = NULL
  for(i in 1:length(structure)){
    nodes = structure[[i]]
    fromrelationindex = which(start %in% nodes)
    torelationindex = which(end %in% nodes)
    
    ends = unique(end[fromrelationindex])
    starts = unique(start[torelationindex])
    
    # Drop used relations
    droprelation = unique(c(fromrelationindex,torelationindex))
    end = end[-droprelation]
    start = start[-droprelation]
    
    startcycles = unlist(lapply(structure, function(x){
      belongs = F
      if(any(x %in% starts)) belongs = T
      return(belongs)
    }))
    endcycles = unlist(lapply(structure, function(x){
      belongs = F
      if(any(x %in% ends)) belongs = T
      return(belongs)
    }))
    
    # Drop out relations from cycle to itself
    endcycles[i] = FALSE
    startcycles[i] = FALSE
    
    startcycle = c(startcycle,which(startcycles))
    endcycle = c(endcycle, rep(i,length(which(startcycles))))
    
    endcycle = c(endcycle,which(endcycles))
    startcycle = c(startcycle, rep(i,length(which(endcycles))))
  }
  return(list(startrelation=startcycle, endrelation=endcycle))
}