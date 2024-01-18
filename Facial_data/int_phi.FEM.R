int_phi.FEM <- function(FEMbasis) {
  
  order     <- FEMbasis$params$order
  nodes     <- FEMbasis$params$nodes
  nodeindex <- FEMbasis$params$t
  Jvec      <- FEMbasis$params$J

  nele  <- dim(nodeindex)[[1]]
  nnod  <- dim(nodes)[[1]]
  
  K0M <- matrix( c( 1,  3,  3), ncol=3, nrow=1, byrow=T) / 6
  
  #  assemble the mass matrix
  
  k0 <- matrix(0,nrow=nnod,ncol=1)
  for (el in 1:nele) {
    ind <- nodeindex[el,]
    k0[ind,1] <- k0[ind,1] + K0M * Jvec[el]
  }
  
  k0<-c(rep(k0, each = 2))
  k0
}

