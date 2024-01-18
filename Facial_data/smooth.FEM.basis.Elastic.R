smooth.FEM.basis.Elastic<- function(FEMloc, FEMdata, FEMbasis, fg, bound,lambda=1e-12) {
  
  FEMdata   <- as.matrix(FEMdata)
  FEMloc<-as.matrix(FEMloc)
  nobs      <- 2*length(FEMdata)
  order     <- FEMbasis$params$order
  nodes     <- FEMbasis$params$nodes
  nodeindex <- FEMbasis$params$t
  Jvec      <- FEMbasis$params$J
  metric    <- FEMbasis$params$metric
  numnodes  <- 2*nrow(nodes)
  indnodes  <- 1:numnodes
  
  eval<-eval.FEM.basis(FEMloc,FEMbasis,nderivs=c(0,0))
  Phimat <- eval[,rep(seq_len(ncol(eval)), each = 2)]
  
  eval2<-eval.FEM.basis(FEMloc,FEMbasis,nderivs=c(1,0))
  Phimat2 <- eval2[,rep(seq_len(ncol(eval2)), each = 2)]
  
  PhitPhimat <- t(Phimat)%*%Phimat
  
  indnodes2  <- 1:(numnodes-bound)
  TK         <- length(indnodes2)
  PhitYmat<-matrix(0,2*length(indnodes2),1)
  PhitYmat[1:TK,] <- t(Phimat[,1:TK])%*%FEMdata
  Phit1<-matrix(0,2*length(indnodes2),nobs)
  Phit1[1:TK,] <- t(Phimat[,1:TK])
  
  IPhi <- int_phi.FEM(FEMbasis)
  PhitYmat[(TK+1):(2*TK),]<-IPhi[1:TK]*fg[1:TK]
  
 
  Amat  <- rbind(cbind(PhitPhimat[1:TK,1:TK], -lambda*k1[1:TK,1:TK]),
                 cbind(k1[1:TK,1:TK],             k0[1:TK,1:TK]))
  
  # solve the linear equations using lsfit
  PhitYmatfit <-solve(Amat,PhitYmat)
  S=t(Phit1)%*%solve(Amat,Phit1)
  edf<-2*Trace(S[1:nrow(nodes),1:nrow(nodes)])-Trace(S[1:nrow(nodes),1:nrow(nodes)]%*%t(S[1:nrow(nodes),1:nrow(nodes)]))
  coefmat <- as.matrix(PhitYmatfit)
  coef1 <- coefmat[indnodes2,]
  coef2 <- as.matrix(coefmat[TK+indnodes2,])
  
  datahat <- Phimat[,indnodes2] %*% coef1
  datahatu<-Phimat[,indnodes2[seq(1,TK,by=2)]] %*% coef1[seq(1,TK,by=2)]
  datahatv<-Phimat[,indnodes2[seq(2,TK,by=2)]]%*% coef1[seq(2,TK,by=2)]
  
  datahatu2<-Phimat2[,indnodes2[seq(1,TK,by=2)]] %*% coef1[seq(1,TK,by=2)]
  datahatv2<-Phimat2[,indnodes2[seq(2,TK,by=2)]]%*% coef1[seq(2,TK,by=2)]
  
  SSE <- sum((FEMdata - datahat)^2)
  
  stderr2 = SSE / ( length(datahat) -edf)
  GCV = ( length(datahat) / ( length(datahat) - edf )) * stderr2
  
  #varf<-stderr2*diag((S%*%S^T))
  varf<-stderr2*diag((S[1:nrow(nodes),1:nrow(nodes)]%*%S[1:nrow(nodes),1:nrow(nodes)]^T))
  #varff=varf[seq(1,nobs,by=2)]+varf[seq(2,nobs,by=2)]
  smoothList <- list(f=datahat,SSE=SSE,coef=coef1,Amat=Amat,edf=edf,GCV=GCV,k0=k0,k1=k1,u=datahatu,v=datahatv,du=datahatu2,dv=datahatv2,var=varf)
  
  return(smoothList)
  
}



