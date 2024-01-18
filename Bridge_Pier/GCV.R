gcv<-function(x,y,FEMbasis,fg,bound,gcv_values){
  # Initialize variables
  gcv_results <- numeric(length(gcv_values))
# Iterate over candidate GCV values
for (i in seq_along(gcv_values)) {
  
  res<-smooth.FEM.basis.Elastic(x,y,FEMbasis,fg,bound,lambda=gcv_values[i])
  
  # Store GCV value in results vector
  gcv_results[i] <- round(res$GCV,10)
}
  gcv_res<-list(gcv_results=gcv_results,gcv_values=gcv_values)
  return(gcv_res)
}
