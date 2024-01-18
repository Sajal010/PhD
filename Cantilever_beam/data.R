


Length = 60; # Length of the model
Width =20;    # Width
NXE = 24;      # Number of rows in the x direction
NYE = 10;      # Number of rows in the y direction
dhx = Length/NXE; # Element size in the x direction
dhy = Width/NYE;  # Element size in the x direction
X_origin = 0 ;  # X origin of the global coordinate system
Y_origin = Width/2 ;   # Y origin of the global coordinate system
#
nne = 3;
nodof = 2;
eldof = nne*nodof;
connec = matrix(0,480,3);
geom<-matrix(0,275,2)
nnd = 0
k = 0

for( i in 1:NXE){
  for (j in 1:NYE)
  {
    k = k + 1;
    n1 = j + (i-1)*(NYE + 1);
    geom[n1,] = cbind((i-1)*dhx - X_origin, (j-1)*dhy - Y_origin );
    n2 = j + i*(NYE+1);
    geom[n2,] = cbind(i*dhx - X_origin, (j-1)*dhy - Y_origin );
    n3 = n1 + 1;
    geom[n3,] = cbind((i-1)*dhx - X_origin, j*dhy - Y_origin );
    n4 = n2 + 1;
    geom[n4,] = cbind(i*dhx- X_origin, j*dhy - Y_origin );
    nel = 2*k;
    m = nel -1;
    connec[m,] = cbind(n1 ,n2 ,n3);
    connec[nel,] = cbind(n2 ,n4 ,n3);
    nnd = n4;
  }
}
nf = matrix(1,nnd, nodof)  # Initialise the matrix nf to 1

for (i in 1:nnd) {
  if (geom[i,1] == Length)
    nf[i,] = c(0, 0)
}



# Counting of the free degrees of freedom
n=0;
for (i in 1:nnd){
  for (j in 1:nodof) {
    if(nf[i,j]!=0){
      n=n+1
      nf[i,j]=n
    }
  }
}

# Initialising load vector
Nodal_loads= matrix(0,nnd, 2)

Force = 1000;  # N

for (i in 1:nnd) {
  if (geom[i,1] == 0. && geom[i,2] == 0)
    Nodal_loads[i,] = cbind(0,  -Force); 
}

fg=matrix(0,n,1)
for (i in 1:nnd) {
  if(nf[i,1]!=0){
    fg[nf[i,1]]=Nodal_loads[i,1]
  }
  if(nf[i,2]!=0){
    fg[nf[i,2]]=Nodal_loads[i,2]
  }
}





