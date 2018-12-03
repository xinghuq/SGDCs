############## Test envrionmental factors that infulence community structure, here Jaccard dissimilarity #######
### comm ,colum is species, row is site/community

glmcomdis_env=function(comm,env,method){
  
  
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
  GlmForm<- as.formula(paste("Comdist~", paste("evdis$",names(evdis)[i],collapse = " + ", sep = ""), sep=" "))
  #GlmForm<- as.formula(paste("PFst ~", paste("evdis[[",i,"]]",collapse = "+",sep = ""), sep=""))
  glmdiff_ev=glm(GlmForm)
  # names(glmdiff_ev[[1]])=c("Intercept",colnames(env))
  summary_glmdiff_ev=summary(glmdiff_ev)
  #model.null = glm(PFst ~ 1)
  model.full=glmdiff_ev
  stepselection=step(model.full)
  #stepselection=step(model.full,scope = list(lower=model.null,upper=model.full),direction="both")
  glmstepselection=summary(stepselection)
  
  #if (all_vifs <- car::vif(stepselection)>2){
  
 # signif_all <- names(all_vifs)
  
  # Remove vars with VIF> 4 and re-build model until none of VIFs don't exceed 4.
 # while(any(all_vifs > 4)){
  #  var_with_max_vif <- names(which(all_vifs == max(all_vifs)))  # get the var with max vif
  #  signif_all <- signif_all[!(signif_all) %in% var_with_max_vif]  # remove
 #   myForm <- as.formula(paste("Comdist ~ ", paste (signif_all, collapse=" + "), sep=""))  # new formula
 #   stepselection <- glm(myForm)  # re-build model with new formula
 #   all_vifs <- car::vif( stepselection)
#  }
 # summary(stepselection)
  
 # }
  
  
  all_vifs <- names(stepselection[[1]])[-1]  # names of all X variables
  
  # Get the non-significant vars
  summ <- summary( stepselection)  # model summary
  pvals <- summ$coefficients[,4]  # get all p values
  not_significant <- character()  # init variables that aren't statsitically significant
  not_significant <- names(which(pvals > 0.1))
  not_significant <- not_significant[!not_significant %in% "(Intercept)"]  # remove 'intercept'. Optional!
  
  # If there are any non-significant variables, 
  while(length(not_significant) > 0){
    all_vifs <- all_vifs[!all_vifs %in% not_significant[1]]
    myForm<- as.formula(paste("Comdist ~ ", paste(all_vifs, collapse=" + "), sep=""))  # new formula
    stepselection <- glm(myForm)  # re-build model with new formula
    
    # Get the non-significant vars.
    summ <- summary( stepselection)
    pvals <- summ$coefficients[,4]
    # not_significant <- character()
    not_significant <- names(which(pvals > 0.1))
    not_significant <- not_significant[!not_significant %in% "(Intercept)"]
  }
  
  final=summary(stepselection)
  
  return(list(glmcom_ev.full=summary_glmdiff_ev,glmstepselection=glmstepselection, final_model=final))
}


