phi_phihat<-function(eval){
  
  nnod  <- dim(eval)[[1]]
  
  
  phi <- matrix(c(rep(c(1,0), each = nnod), rep(c(0,1), each = nnod)), nrow = 2*nnod, ncol = 2*nnod)
  # split the columns of m into a list based on the column index modulo 3
  l <- split(eval, rep(1:nnod, each=nnod))
  
  # replicate each element of the list twice
  l <- lapply(l, function(x) rep(list(x), 2))
  
  # unlist the list to obtain a flattened list and convert it to a matrix
  result <- matrix(unlist(l), nrow=2*nnod, byrow=TRUE)
  
  # print the result
  phi[phi == 1] <- matrix(t(result) ,ncol = 1)
  
  phi_final<-phi%*%t(phi)
}

phi<-function(eval){
  
  nnod  <- dim(eval)[[1]]
  
  
  phi <- matrix(c(rep(c(1,0), each = nnod), rep(c(0,1), each = nnod)), nrow = 2*nnod, ncol = 2*nnod)
  # split the columns of m into a list based on the column index modulo 3
  l <- split(eval, rep(1:nnod, each=nnod))
  
  # replicate each element of the list twice
  l <- lapply(l, function(x) rep(list(x), 2))
  
  # unlist the list to obtain a flattened list and convert it to a matrix
  result <- matrix(unlist(l), nrow=2*nnod, byrow=TRUE)
  
  # print the result
  phi[phi == 1] <- matrix(t(result) ,ncol = 1)
  
  phi
}



