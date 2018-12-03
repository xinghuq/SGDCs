####Multiple Mantel regression model to select the significnat factors

############## This function test envrionmental factors that infulence Belta SGDCs#######

# envrionmental factors that infulence population structure (Fst)



multRegFst_env=function(g,env,ncode,nperm,sig){
  read.genepop1 <- function(file, ncode, quiet = TRUE) {
    adegenet::.readExt
    adegenet::.genlab
    adegenet::df2genind
    adegenet::is.genind
    adegenet::pop
    adegenet::repool
    adegenet::Hs
    adegenet::seppop
    adegenet::popNames
    ade4::mantel.randtest
    if (toupper(.readExt(file)) != "GEN")
      stop("File extension .gen expected")
    if (!quiet)
      cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")
    prevcall <- match.call()
    txt <- scan(file, sep = "\n", what = "character", quiet = TRUE)
    if (!quiet)
      cat("\nFile description: ", txt[1], "\n")
    txt <- txt[-1]
    txt <- gsub("\t", " ", txt)
    NA.char <- paste(rep("0", ncode), collapse = "")
    locinfo.idx <- 1:(min(grep("POP", toupper(txt))) - 1)
    locinfo <- txt[locinfo.idx]
    locinfo <- paste(locinfo, collapse = ",")
    loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
    loc.names <- trimws(loc.names)
    nloc <- length(loc.names)
    txt <- txt[-locinfo.idx]
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
    npop <- length(pop.idx)
    nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
    splited <- nocomma[which(!nocomma %in% pop.idx)]
    if (length(splited) > 0) {
      for (i in sort(splited, decreasing = TRUE)) {
        txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
      }
      txt <- txt[-splited]
    }
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
    txt[length(txt) + 1] <- "POP"
    nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))) - 1
    pop <- factor(rep(1:npop, nind.bypop))
    txt <- txt[-c(pop.idx, length(txt))]
    temp <- sapply(1:length(txt), function(i) strsplit(txt[i], ","))
    ind.names <- vapply(temp, function(e) e[1], character(1))
    ind.names <- trimws(ind.names)
    vec.genot <- vapply(temp, function(e) e[2], character(1))
    vec.genot <- trimws(vec.genot)
    X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), ncol = nloc, byrow = TRUE)
    if (any(duplicated(ind.names))) {
      rownames(X) <- .genlab("", nrow(X))
    } else {
      rownames(X) <- ind.names
    }
    colnames(X) <- loc.names
    pop.names.idx <- cumsum(table(pop))
    pop.names <- ind.names[pop.names.idx]
    levels(pop) <- pop.names
    if (!all(unique(nchar(X)) == (ncode * 2)))
      stop(paste("some alleles are not encoded with", ncode, "characters\nCheck 'ncode' argument"))
    res <- df2genind(X = X, pop = as.character(pop), ploidy = 2, ncode = ncode, NA.char = NA.char)
    res@call <- prevcall
    if (!quiet)
      cat("\n...done.\n\n")
    return(res)
  }
  g0 = read.genepop1(g, ncode, quiet = TRUE)
  # genind file
  hierfstat::genind2hierfstat
  g1=genind2hierfstat(g0)
  hierfstat::pairwise.WCfst
  PFst = pairwise.WCfst(g1)
  PFst=as.dist(PFst, diag = FALSE, upper = FALSE)
  rownames(env) = env[,1]
  env=as.matrix(env[,-1])
  evdis=setNames(as.list(c(1:ncol(env))),c(colnames(env)))
  for (i in 1:ncol(env)){
    # evdis[[i]]=list()
    evdis[[i]]=dist(env[,i],method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  }
  i=1:ncol(env)
  mrForm<- as.formula(paste("PFst ~", paste("evdis$",names(evdis)[i],collapse = " + ", sep = ""), sep=" "))
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
    
    mrForm<- as.formula(paste("PFst ~ ", paste(all_vifs, collapse=" + "), sep=""))  # new formula
    mrdiff_ev<- MRM(mrForm,method="linear", nperm=nperm,mrank = FALSE)  # re-build model with new formula
    
    # Get the non-significant vars.
    # if(pval[1]<0.05) break
    pvals <- mrdiff_ev$coef[,2]
    # not_significant <- character()
    not_significant <- names(which(pvals > sig))
    not_significant <- not_significant[!not_significant %in% "Int"]
    not_significant=sort(not_significant,decreasing=T)
  }
  final=mrdiff_ev
  
  return(list(full.model=model.full,final_model=final))
}
