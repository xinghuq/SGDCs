############## This function calculate beta SGDCs, which correlates community dissimilarity measured by different methods and genetic differentiation measured by Fst #######


ComDis_Fst=function(comm,g,ncode,method,nrepet){
  
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
  vegan::vegdist
  
 rownames(comm)=comm[,1]
 comm=comm[,-1]
  Comdist=vegdist(comm, method, binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
  Commdist_Fst=mantel.randtest(Comdist, PFst,nrepet)
  return(list(PairwiseFst = PFst, pairwiseDcom = Comdist, Commdist_Fst = Commdist_Fst))
}


