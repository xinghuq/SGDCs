############## This function calculate beta SGDCs, which correlates community species differentiation measured by Delta D and genetic differentiation measured by Fst #######


betaSGDCFst=function(comm,x,ncode,nrepet){
  
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
  g = read.genepop1(x, ncode, quiet = TRUE)
  # genind file
  hierfstat::genind2hierfstat
  g1=genind2hierfstat(g)
  hierfstat::pairwise.WCfst
  PFst = pairwise.WCfst(g1)
  PFst=as.dist(PFst, diag = FALSE, upper = FALSE)
  
  #npops = length(levels(g1$pop))
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
  
  v1 = c("ecosystem", "region1", "com1")
  v2 = c("ecosystem", "region1", "com2")
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
  
  betaSGDC_Fst=mantel.randtest(Dcom, PFst,nrepet)
  
  return(list(PairwiseFst = PFst, pairwiseDcom = Dcom, betaSGDC = betaSGDC_Fst))
}


