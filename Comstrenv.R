############## Test envrionmental factors that infulence community structure #######


Comstr_env=function(comm,env,method,nrepet){
  
  
  vegan::vegdist

  rownames(comm)=comm[,1]
  comm=comm[,-1]
  Comdist=vegdist(comm, method, binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
  
  rownames(env) = env[,1]
  env=as.matrix(env[,-1])
  evdis=list()
  for (i in 1:ncol(env)){
    evdis[[i]]=list()
    evdis[[i]]=dist(env[,i],method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
    
  }
  # evdis=dist(env,method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  ######## here we should induce permutation that test their signiciance
  
  #########in completed ,check again using mantel test
  ade4::mantel.randtest
  comevcor=setNames(as.list(c(1:ncol(env))),c(colnames(env)))
  
  for (i in seq_along(evdis)) {
    comevcor[[i]]=mantel.randtest(Comdist,evdis[[i]], nrepet)
  }
  return(list(corComdis_env = comevcor))
}


