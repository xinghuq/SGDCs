####  Auto correlaton or auto corregram for envrionmental factors,Mantel correlograms.

Envcorgram=function(ev,geo,nperm){
  
  rownames(env) = env[,1]
  env=as.matrix(env[,-1])
  evdis=setNames(as.list(c(1:ncol(env))),c(colnames(env)))
  for (i in 1:ncol(env)){
    # evdis[[i]]=list()
    evdis[[i]]=dist(env[,i],method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  }
  
  i=1:ncol(env)
  geo=as.dist(geo)

mgram(ev, geo, breaks, nclass, stepsize, nperm,mrank = FALSE, alternative = "two.sided", trace = FALSE)
  
}


spcorgram=function(comm,geo,method,nperm){
  
  vegan::vegdist
  
  rownames(comm)=comm[,1]
  comm=comm[,-1]
  Comdist=vegdist(comm, method, binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
 
  geo=as.dist(geo)
  
  spcogram=mgram(Comdist, geo,nperm,mrank = FALSE, alternative = "two.sided", trace = FALSE)
  return(list(spcogram=spcorgram,plot=plot(spcorgram)))
}
