##### function for calculating community species differentiation, pairwise community species differentiation########


pwspdiff=function(comm){

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

return(Dcom)

}