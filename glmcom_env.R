################## Glm to detect the factors drive comm structure

############## Test envrionmental factors that infulence community structure #######

#### comm data frame is diff from comstrenv function

glmcomD_env=function(comm,env){
  
  DeltaD = function(abun, struc) {
    ## Chao et al, 2017
    n = sum(abun)
    N = ncol(abun)
    ga = rowSums(abun)
    gp = ga[ga > 0]/n
    G = sum(-gp * log(gp))
    H = nrow(struc)
    A = numeric(H - 1)
    W = numeric(H - 1)
    Diff = numeric(H - 1)
    wi = colSums(abun)/n
    W[H - 1] = -sum(wi[wi > 0] * log(wi[wi > 0]))
    pi = sapply(1:N, function(k) abun[, k]/sum(abun[, k]))
    Ai = sapply(1:N, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0])))
    A[H - 1] = sum(wi * Ai)
    if (H > 2) {
      for (i in 2:(H - 1)) {
        I = unique(struc[i, ])
        NN = length(I)
        ai = matrix(0, ncol = NN, nrow = nrow(abun))
        c
        for (j in 1:NN) {
          II = which(struc[i, ] == I[j])
          if (length(II) == 1) {
            ai[, j] = abun[, II]
          } else {
            ai[, j] = rowSums(abun[, II])
          }
        }
        pi = sapply(1:NN, function(k) ai[, k]/sum(ai[, k]))
        wi = colSums(ai)/sum(ai)
        W[i - 1] = -sum(wi * log(wi))
        Ai = sapply(1:NN, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0])))
        A[i - 1] = sum(wi * Ai)
      }
    }
    Diff[1] = (G - A[1])/W[1]
    if (H > 2) {
      for (i in 2:(H - 1)) {
        Diff[i] = (A[i - 1] - A[i])/(W[i] - W[i - 1])
      }
    }
    Diff = Diff
    out = matrix(c(Diff), ncol = 1)
    return(out)
  }
  
  v1 = c("ecosystem", "region1", "pop1")
  v2 = c("ecosystem", "region1", "pop2")
  str = data.frame(v1, v2)
  str = as.matrix(str)
  
  
  ncom = nrow(comm)
  nsp = ncol(comm)
  Dsp = matrix(data = 0, nrow = nsp, ncol = nsp)
  
  for (i in 1:nsp) {
    for (j in 1:nsp) {
      Dsp[i,j] = DeltaD(comm[,c(i, j)], str)[2]  ### select two pops from allelefrequency
    }
  }
  
  colnames(Dsp) = colnames(comm)
  rownames(Dsp) = colnames(comm)
  
  Dcom=as.dist(Dsp, diag = FALSE, upper = FALSE)
  
  rownames(env) = env[,1]
  env=as.matrix(env[,-1])
  evdis=setNames(as.list(c(1:ncol(env))),c(colnames(env)))
  for (i in 1:ncol(env)){
    # evdis[[i]]=list()
    evdis[[i]]=dist(env[,i],method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  }
  i=1:ncol(env)
  GlmForm<- as.formula(paste("Dcom ~", paste("evdis$",names(evdis)[i],collapse = " + ", sep = ""), sep=" "))
  #GlmForm<- as.formula(paste("PFst ~", paste("evdis[[",i,"]]",collapse = "+",sep = ""), sep=""))
  glmdiff_ev=glm(GlmForm)
  # names(glmdiff_ev[[1]])=c("Intercept",colnames(env))
  summary_glmdiff_ev=summary(glmdiff_ev)
  #model.null = glm(PFst ~ 1)
  model.full=glmdiff_ev
  stepselection=step(model.full)
  #stepselection=step(model.full,scope = list(lower=model.null,upper=model.full),direction="both")
  glmstepselection=summary(stepselection)
  
  all_vifs <- car::vif(stepselection)
  
  signif_all <- names(all_vifs)
  
  # Remove vars with VIF> 4 and re-build model until none of VIFs don't exceed 4.
  while(any(all_vifs > 4)){
    var_with_max_vif <- names(which(all_vifs == max(all_vifs)))  # get the var with max vif
    signif_all <- signif_all[!(signif_all) %in% var_with_max_vif]  # remove
    myForm <- as.formula(paste("Dcom ~ ", paste (signif_all, collapse=" + "), sep=""))  # new formula
    stepselection <- glm(myForm)  # re-build model with new formula
    all_vifs <- car::vif( stepselection)
  }
  summary( stepselection)
  
  
  
  
  all_vifs <- names( stepselection[[1]])[-1]  # names of all X variables
  
  # Get the non-significant vars
  summ <- summary( stepselection)  # model summary
  pvals <- summ$coefficients[,4]  # get all p values
  not_significant <- character()  # init variables that aren't statsitically significant
  not_significant <- names(which(pvals > 0.1))
  not_significant <- not_significant[!not_significant %in% "(Intercept)"]  # remove 'intercept'. Optional!
  
  # If there are any non-significant variables, 
  while(length(not_significant) > 0){
    all_vifs <- all_vifs[!all_vifs %in% not_significant[1]]
    myForm<- as.formula(paste("Dcom ~ ", paste(all_vifs, collapse=" + "), sep=""))  # new formula
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



