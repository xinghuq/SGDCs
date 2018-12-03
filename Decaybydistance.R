

DBD=function(comm,geo,method,nrepet){
  
  vegan::vegdist
  
  rownames(comm)=comm[,1]
  comm=comm[,-1]
  Comdist=vegdist(comm, method, binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
  
  
  
  geo=as.dist(geo, diag = FALSE, upper = FALSE)
  # evdis=dist(env,method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  ######## here we should induce permutation that test their signiciance
  
  #########in completed ,check again using mantel test
  ade4::mantel.randtest
  
  comDgeo=mantel.randtest(Comdist,geo, nrepet)
  
  return(list(corCom_geo = comDgeo))
}
