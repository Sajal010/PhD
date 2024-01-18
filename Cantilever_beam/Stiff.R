stiff.FEM <- function(FEMbasis) {
  # STIFF.FEM produces the nnod*nnod stiffness matrix K1
  # defined (K1)jk = int(dpsik/da*dpsij/da + dpsik/db*dpsij/db).
  #
  # Input: FEMbasis is a List object produced by function makenodes.
  #    It contains:
  #        ORDER     ... The order of the element (1 or 2)
  #        NODES     ... Coordinates of node points
  #        NODEINDEX ... indices of node points for each element
  #        JVEC      ... Jacobian of the affine transformation of each
  #                      element to the master element
  #        METRIC    ... The crossproduct of the inverse of the linear
  #                      part of the transformation
  #
  # Output: K1 is an nnod*nnod matrix out which is
  #        the sum of the nele element stiffness matrices
  #        and the penalty stiffness matrix.
  #        These i'th element matrix has (ij)'th element defined
  #        as follows:
  #        Let psita and psitb be the partial derivatives of the
  #        t'th shape function with respect to a and b (1<=t<=6).
  #        Then the integral of the sum of products
  #        (psija*psika+psijb+psikb) over the i'th element is
  #        computed.  Then that value is assigned to the
  #        (nodeindex(i,j),nodeindex(i,k))'th entry of the i'th elemental
  #        stiffness matrix and the other elements are given the value zero.
  #
  #
  #  Last modified 19 November 2021 by Jim Ramsay.
  
  ###For Coarse data
  E = 200000.;     # Elastic modulus in MPa
   vu = 0.3;        # Poisson's ratio 
  thick = 5.;      # Beam thickness in mm
  
  ###For Cantilver beam
  #E = 200000.;     # Elastic modulus in MPa
 #vu = 0.3;        # Poisson's ratio 
  #thick = 5.;      # Beam thickness in mm
  
  # For facial data
  #E = 0.015;     # Elastic modulus in MPa
  #vu = 0.49;        # Poisson's ratio 
 # thick = 1.51;      # Beam thickness in mm
  
  #Form the elastic matrix for plane stress 
  
  # Function
  C.FEM=function(E,nu){
    
    c=E/(1-nu^2)
    
    d=c* matrix( c(  1,  nu, 0,
                     nu,  1,  0,
                     0,  0,  0.5*(1-nu)), ncol=3, nrow=3, byrow=T)
    d
  }
  
  dee = C.FEM(E,vu)
  
  #  retrieve arrays from FEMbasis
  nne = 3
  nel<-dim(FEMbasis$params$t)[[1]]
  geom<-    FEMbasis$params$p
  connec<-FEMbasis$params$t
  order     <- FEMbasis$params$order
  nodes     <- FEMbasis$params$nodes
  nodeindex <- FEMbasis$params$nodeindex
  Jvec      <- FEMbasis$params$J
  metric    <- FEMbasis$params$metric
  
  nele  <- dim(nodeindex)[[1]]
  nnod  <- dim(nodes)[[1]]
  nodof <-  dim(nodes)[[2]]
  eldof = nne*nodof;
  
  nf = matrix(1,nnod, nodof) 
  
  # Counting of the free degrees of freedom
  n=0;
  for (i in 1:nnod){
    for (j in 1:nodof) {
      if(nf[i,j]!=0){
        n=n+1
        nf[i,j]=n
      }
    }
  }
  
  elem<- function(i){
    
    x1 = geom[connec[i,1],1];   y1 = geom[connec[i,1],2];
    x2 = geom[connec[i,2],1];   y2 = geom[connec[i,2],2];
    x3 = geom[connec[i,3],1];   y3 = geom[connec[i,3],2];
    
    x=matrix(c(1,x1,y1,
               1,x2,y2,
               1,x3,y3),3,3);
    A = (0.5)*det(x);
    
    
    m11 = (x2*y3 - x3*y2)/(2*A);
    m21 = (x3*y1 - x1*y3)/(2*A);
    m31 = (x1*y2 - y1*x2)/(2*A);
    m12 = (y2 - y3)/(2*A);
    m22 = (y3 - y1)/(2*A);
    m32 = (y1 - y2)/(2*A);
    m13 = (x3 - x2)/(2*A);
    m23 = (x1 - x3)/(2*A);
    m33 = (x2 -x1)/(2*A);
    
    bee = matrix(c(m12,0 ,m22,0,m32,0,
                   0 ,m13  , 0 ,m23  , 0 , m33,
                   m13 ,m12, m23 ,m22  ,m33  ,m32), 3, 6,byrow=TRUE);
    
    l=1;
    g = matrix(NA,1,nne*nodof); 
    for (k in 1:nne) {
      for (j in 1:nodof) {
        g[l]=nf[connec[i,k],j];
        l=l+1;
      }
    }
    return(list(bee=bee,g=g,A=A))
  }
  
  ## form_kk function
  kk = matrix(0,n,n);
  form_kk=function(kk,kg,g){
    for (i in 1:eldof) {
      if(g[i]!=0)
        for (j in 1:eldof) {
          if(g[j]!=0)
            kk[g[i],g[j]]= kk[g[i],g[j]]+kg[i,j];
        }
    }
    return(kk)
  }
  
  ## KK matrix
  for (i in 1:nel){
    k=elem(i);
    ke=thick*k$A*t(k$bee)%*%dee%*%k$bee;
    kk=form_kk(kk,ke,k$g);
  }
  
  kk
  
  
}



