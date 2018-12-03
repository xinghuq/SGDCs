################## Do multiple regression to detect the factors drive comm structure

############## Test envrionmental factors that infulence community structure , here Delta D#######

#### comm data frame is diff from comstrenv function

multRegcomD_env=function(comm,env,nperm,sig){
  
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
  
  mrForm<- as.formula(paste("Dcom ~", paste("evdis$",names(evdis)[i],collapse = " + ", sep = ""), sep=" "))
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
      
    mrForm<- as.formula(paste("Dcom ~ ", paste(all_vifs, collapse=" + "), sep=""))  # new formula
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



