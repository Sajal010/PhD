# Setting working directory
setwd("~/Desktop/Phd Simulations/Cantilever_beam")

# Loading all the relevant libraries
library(geometry)
library(plotly)
library(fdaPDE)
library(sp)
library(concaveman)
library("fda")
library("pracma")

# Sourcing functions
source("make_nodes.R")
source("create.FEM.basis.R")
source("eval.FEM.basis.R")
source("insideIndex.R")
source("int_phi.FEM.R")
source("Mass.R")
source("phi_phihat.R")
source("Stiff.R")
source("tricoefCal.R")
source("GCV.R")
source("smooth.FEM.basis.Elastic.R")
source("data.R")


# Loading dataset
Rdata<-read.csv("true_z.csv",header = FALSE)
obs<-Rdata[,1]+Rdata[,2]

# Plotting original data
d<-data.frame(x=geom[,1],y=geom[,2],z=obs)
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig


# Scaling the true z values for the simulation
FEMdata = matrix(0,2*nrow(Rdata),1)
FEMdata[seq(1,2*nrow(Rdata),by=2),1]=Rdata[,1]
FEMdata[seq(2,2*nrow(Rdata),by=2),1]=Rdata[,2]


xx<-FEMdata[seq(1,2*nrow(Rdata),by=2),1]+FEMdata[seq(2,2*nrow(Rdata),by=2),1]+
  rnorm(nrow(Rdata),0,0.02*diff(range(FEMdata[seq(1,2*nrow(Rdata),by=2),1]+FEMdata[seq(2,2*nrow(Rdata),by=2),1])))

# Plotting scaled z values
d<-data.frame(x=geom[,1],y=geom[,2],z=xx)
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# Create a SpatialPolygons object representing the boundary polygon
boundary <- concaveman(geom)
boundary_polygon <- SpatialPolygons(list(Polygons(list(Polygon(boundary)), ID = "boundary")))
geom<-data.frame(geom)

# Check if points are inside or on the boundary of the shape
is_inside <- point.in.polygon(geom[,1], geom[,2], boundary_polygon@polygons[[1]]@Polygons[[1]]@coords[,1],
                              boundary_polygon@polygons[[1]]@Polygons[[1]]@coords[,2])
inside_points <- geom[is_inside == 1, ]
boundary_points <- geom[is_inside >=2, ]


# Creating mesh of the domain
mesh=create.mesh.2D(nodes = geom, segments = convhulln(geom))
plot(mesh)

# Creating finite element basis functions
FEMbasis<-create.FEM.basis(mesh$nodes,mesh$segments,mesh$triangles,1,0)

# Selection of optimal lambda using GCV
gcv_values <- 10^(seq(-9, 10, by = 1))

gcv_optimise<-gcv(as.matrix(geom),xx,FEMbasis,fg,22,gcv_values)

# Print the optimal GCV value
cat("Optimal GCV value:", gcv_optimise$gcv_values[which.min(gcv_optimise$gcv_results)])
optimal_gcv<-gcv_optimise$gcv_values[which.min(gcv_optimise$gcv_results)]

# Plot the GCV values
plot(gcv_optimise$gcv_values[1:5], gcv_optimise$gcv_results[1:5], type = "b", pch = 16, 
     xlab = "GCV Value", ylab = "GCV", main = "GCV Values")
points(optimal_gcv,gcv_optimise$gcv_results[which.min(gcv_optimise$gcv_results)],pch = 16,col="red",cex=1.5)

# Running the smooth.FEM function
fd<-smooth.FEM.basis.Elastic(as.matrix(geom),xx,FEMbasis,fg,22,lambda=optimal_gcv)

# Plotting fit of the cantilever beam in 3D
d<-data.frame(x=geom[,1],y=geom[,2],z= fd$f)
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# Plotting true data of the cantilever beam in 3D
d<-data.frame(x=geom[,1],y=geom[,2],z= Rdata[,1]+Rdata[,2])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# Plotting the fit of displacement u in x-direction of the cantilever beam in 3D
d<-data.frame(x=geom[,1],y=geom[,2],z= fd$u)
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# Plotting the truth of displacement u in x-direction of the cantilever beam in 3D
d<-data.frame(x=geom[,1],y=geom[,2],z= Rdata[,1])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# Plotting the fit of displacement v in y-direction of the cantilever beam in 3D
d<-data.frame(x=geom[,1],y=geom[,2],z= fd$v)
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# Plotting the truth of displacement v in y-direction of the cantilever beam in 3D
d<-data.frame(x=geom[,1],y=geom[,2],z= Rdata[,2])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# Plotting the fit of sum of the two displacements u and v of the cantilever beam in 3D
d<-data.frame(x=geom[,1],y=geom[,2],z= fd$u+fd$v)
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

###~~~~ Calculating stresses and Strains
###For Cantilver beam
E = 200000.;     # Elastic modulus in MPa
vu = 0.3;        # Poisson's ratio 
thick = 5.;      # Beam thickness in mm

# Function
C.FEM=function(E,nu){
  
  c=E/(1-nu^2)
  
  d=c* matrix( c(  1,  nu, 0,
                   nu,  1,  0,
                   0,  0,  0.5*(1-nu)), ncol=3, nrow=3, byrow=T)
  d
}

dee = C.FEM(E,vu)

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


delta<-matrix(0,2*length(fd$u),1)
delta[seq(1,2*length(fd$u),by=2)]=fd$u[c(1:length(fd$u))]
delta[seq(2,2*length(fd$v),by=2)]=fd$v[c(1:length(fd$v))]

nne=3
nodes     <- FEMbasis$params$nodes
nodeindex <- FEMbasis$params$nodeindex

nnd  <- dim(nodes)[[1]]
nodof <-  dim(nodes)[[2]]
nf = matrix(1,nnd, nodof) 
nele  <- dim(nodeindex)[[1]]
connec<-FEMbasis$params$t
eldof = nne*nodof

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


node_disp<-matrix(0,nnd,2)
for (i in 1:nnd) {
  {if(nf[i,1]==0){
    x_disp =0
  }
    else {
      x_disp = delta[nf[i,1]]
    }
  }
  {if(nf[i,2]==0){
    y_disp =0
  }
    else {
      y_disp = delta[nf[i,2]]
    }
  }
  node_disp[i,] =c(x_disp, y_disp);
}

EPS<-matrix(0,nele,3)
SIGMA<-matrix(0,nele,3)
for (i in 1:nele){
  k=elem(i);
  eld<-matrix(0,eldof,1)
  for (m in 1:eldof) {
    if(k$g[m]==0){
      eld[m]=0
    }
    else
    {
      eld[m]=delta[k$g[m]] 
    }
  }
  eps=k$bee%*%eld # Compute strains
  EPS[i,]=eps # Store strains for all elements
  sigma=dee%*%eps; # Compute stresses
  SIGMA[i,]=sigma ; # Store stresses for all elements
}


#~~~~ Stresses
# x-direction
d<-data.frame(x=nodes[nodeindex[,1],1],y=nodes[nodeindex[,1],2],z=SIGMA[,1])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# y-direction
d<-data.frame(x=nodes[nodeindex[,1],1],y=nodes[nodeindex[,1],2],z=SIGMA[,2])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# xy-direction
d<-data.frame(x=nodes[nodeindex[,1],1],y=nodes[nodeindex[,1],2],z=SIGMA[,3])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

#~~~~ Strains
# x-direction
d<-data.frame(x=nodes[nodeindex[,1],1],y=nodes[nodeindex[,1],2],z=EPS[,1])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# y-direction
d<-data.frame(x=nodes[nodeindex[,1],1],y=nodes[nodeindex[,1],2],z=EPS[,2])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig

# xy-direction
d<-data.frame(x=nodes[nodeindex[,1],1],y=nodes[nodeindex[,1],2],z=EPS[,3])
fig <- plot_ly(d, x = ~x, y = ~y, z = ~z, color = ~z)
fig



