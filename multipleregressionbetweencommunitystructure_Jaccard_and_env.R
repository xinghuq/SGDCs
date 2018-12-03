############## Test envrionmental factors that infulence community structure, here Jaccard dissimilarity #######
### comm ,colum is species, row is site/community

multRegcomdis_env=function(comm,env,method,nperm,sig){
  
  
  vegan::vegdist
  
  rownames(comm)=comm[,1]
  comm=comm[,-1]
  Comdist=vegdist(comm, method, binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
  
  rownames(env) = env[,1]
  env=as.matrix(env[,-1])
  
  evdis=setNames(as.list(c(1:ncol(env))),c(colnames(env)))
  for (i in 1:ncol(env)){
    # evdis[[i]]=list()
    evdis[[i]]=dist(env[,i],method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  }
  i=1:ncol(env)
  
  mrForm<- as.formula(paste("Comdist~", paste("evdis$",names(evdis)[i],collapse = " + ", sep = ""), sep=" "))
  #GlmForm<- as.formula(paste("PFst ~", paste("evdis[[",i,"]]",collapse = "+",sep = ""), sep=""))
 
  
  ecodist::MRM
  
  mrdiff_ev=MRM(mrForm,method="linear",nperm=nperm, mrank = FALSE)
  model.full=mrdiff_ev
  all_vifs=rownames(mrdiff_ev$coef)[-1] # names of all X variables
  # Get the non-significant vars
  pvals <- mrdiff_ev$coef[,2]# get all p values
  pavls=sort(pvals,decreasing=T)
  not_significant <- character()  # init variables that aren't statsitically significant
  not_significant <- names(which(pvals > sig))
  not_significant <- not_significant[!not_significant %in% "Int"]  # remove 'intercept'. Optional!
  not_significant=sort(not_significant,decreasing=T)
  # If there are any non-significant variables, 
  # if(length(not_significant)> 0){
  while(length(not_significant) > 0){
    all_vifs <- all_vifs[!all_vifs %in% not_significant[1]]
    if (length(all_vifs)< 1) {mrdiff_ev=c("No significant value")} 
    if (length(all_vifs)< 1)break 
    
    mrForm<- as.formula(paste("Comdist ~ ", paste(all_vifs, collapse=" + "), sep=""))  # new formula
    mrdiff_ev<- MRM(mrForm,method="linear", nperm=nperm,mrank = FALSE)  # re-build model with new formula
    
    # Get the non-significant vars.
    # if(pval[1]<0.05) break
    pvals <- mrdiff_ev$coef[,2]
    # not_significant <- character()
    not_significant <- names(which(pvals > sig))
    not_significant <- not_significant[!not_significant %in% "Int"]
    not_significant=sort(not_significant,decreasing=T)
  }
  # }
  
  final=mrdiff_ev
  return(list(full.model=model.full,final_model=final))
  }



