smooth.FEM.basis.Elastic <- function(FEMloc, FEMdata, FEMbasis, fg, bound, lambda=1e-12) {
  
  FEMdata   <- as.matrix(FEMdata)
  nobs      <- 2*length(FEMdata)
  order     <- FEMbasis$params$order
  nodes     <- FEMbasis$params$nodes
  nodeindex <- FEMbasis$params$t
  Jvec      <- FEMbasis$params$J
  metric    <- FEMbasis$params$metric
  numnodes  <- 2*nrow(nodes)
  indnodes  <- 1:numnodes
  
  # Mass matrix
  k0 <- mass.FEM(FEMbasis)
  k0 <- k0[rep(seq_len(ncol(k0)), each = 2),rep(seq_len(ncol(k0)), each = 2)]
  
  # Stiffness matrix
  k1 <- stiff.FEM(FEMbasis)
  
  eval   <- eval.FEM.basis(FEMloc,FEMbasis,nderivs=c(0,0))
  Phimat <- eval[,rep(seq_len(ncol(eval)), each = 2)]
  
  PhitPhimat <- t(Phimat)%*%Phimat
  
  indnodes2  <- 1:nrow(nodes)
  indnodes2  <- indnodes2[-bound]
  indnodes2  <- sort(c(indnodes2*2,indnodes2*2-1))
  
  TK         <- length(indnodes2)
  PhitYmat<-matrix(0,2*TK,1)
  PhitYmat[1:TK,] <- t(Phimat[,indnodes2])%*%FEMdata
  Phit1<-matrix(0,2*TK,nobs)
  Phit1[1:TK,] <- t(Phimat[,indnodes2])
  
  IPhi <- int_phi.FEM(FEMbasis)
  PhitYmat[(TK+1):(2*TK),]<-(IPhi[indnodes2]*fg[indnodes2])
  
  # Create A
  Amat  <- rbind(cbind(PhitPhimat[indnodes2,indnodes2], -lambda*k1[indnodes2,indnodes2]),
                 cbind(k1[indnodes2,indnodes2],            k0[indnodes2,indnodes2]))
  
  # solve the linear equations using lsfit
  PhitYmatfit <-solve(Amat,PhitYmat)

  coefmat <- as.matrix(PhitYmatfit)
  coef1 <- coefmat[1:TK,]
  coef2 <- coefmat[(length(coef1)+1):length(coefmat),]
  
  datahat <- Phimat[,indnodes2] %*% coef1
  
  datahatu<-Phimat[,indnodes2[seq(1,TK,by=2)]] %*% coef1[seq(1,TK,by=2)]
  datahatv<-Phimat[,indnodes2[seq(2,TK,by=2)]] %*% coef1[seq(2,TK,by=2)]
  
  eval2<-eval.FEM.basis(FEMloc,FEMbasis,nderivs=c(1,0))
  Phimat2 <- eval2[,rep(seq_len(ncol(eval2)), each = 2)]
  
  datahatu2<-Phimat2[,indnodes2[seq(1,TK,by=2)]] %*% coef1[seq(1,TK,by=2)]
  datahatv2<-Phimat2[,indnodes2[seq(2,TK,by=2)]]%*% coef1[seq(2,TK,by=2)]
  
  SSE <- sum((FEMdata - datahat)^2)
 
  S   = t(Phit1)%*%solve(Amat,Phit1)
  edf<-2*Trace(S[1:nrow(nodes),1:nrow(nodes)])-Trace(S[1:nrow(nodes),1:nrow(nodes)]%*%t(S[1:nrow(nodes),1:nrow(nodes)]))
  
  stderr2 = SSE / ( length(datahat) -edf)
  GCV = ( length(datahat) / ( length(datahat) - edf )) * stderr2
  
  varf<-stderr2*diag((S[1:nrow(nodes),1:nrow(nodes)]%*%S[1:nrow(nodes),1:nrow(nodes)]^T))
  smoothList <- list(f=datahat,SSE=SSE,coef=coef1,Amat=Amat,edf=edf,GCV=GCV,k0=k0,k1=k1,u=datahatu,v=datahatv,du=datahatu2,dv=datahatv2,var=varf)

  return(smoothList)
  
}



