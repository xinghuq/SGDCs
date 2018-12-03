####glm model to select the significnat factors

############## This function test envrionmental factors that infulence beta SGDCs#######

# envrionmental factors that infulence population structure (Fst)



glmFst_env=function(g,env,ncode){
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
  GlmForm<- as.formula(paste("PFst ~", paste("evdis$",names(evdis)[i],collapse = " + ", sep = ""), sep=" "))
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
    myForm <- as.formula(paste("PFst ~ ", paste (signif_all, collapse=" + "), sep=""))  # new formula
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
    myForm<- as.formula(paste("PFst ~ ", paste(all_vifs, collapse=" + "), sep=""))  # new formula
     stepselection <- glm(myForm)  # re-build model with new formula
    
    # Get the non-significant vars.
    summ <- summary( stepselection)
    pvals <- summ$coefficients[,4]
   # not_significant <- character()
    not_significant <- names(which(pvals > 0.1))
    not_significant <- not_significant[!not_significant %in% "(Intercept)"]
  }
  
  final=summary( stepselection)
  
  return(list(glmFst_ev.full=summary_glmdiff_ev,glmstepselection=glmstepselection, final_model=final))
}
