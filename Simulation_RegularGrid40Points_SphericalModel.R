#---------------------------------------------------------
# R code to carry out the Mantel test for functional
# data based on simulated data
# Grid: Regular
# Points:40
# Covariance Model: Spherical
# ---------------------------------------------------------

rm(list=ls())
library(ade4)
library(geoR)
library(MASS)
source("integral.R")      # Trapecio's rule for evaluating the integral

a=10             # median parameter a (this is fixed in the simulation)
mb=0.7           # median b parameter  (this would be simulated from a independent Gaussian field)
sdb=0.15         # standard deviation b
sites=40         # number of sites
mt=180           # max time in minutes

grid<-expand.grid(seq(0,80, length=5),seq(0,140,length=8))
plot(grid, xlab="", ylab="", pch=19, main="")
title(main=list("Grid 1 (Regular 40 Points)", cex=2, font =1) , xlab=list("x", cex=2, font =1), ylab=list("y", cex=2, font =1))


#  Parameters of the spatial correlation

sigma2=sdb^2
mean=rep(mb,sites)
distance_coords=dist(grid, diag=TRUE, upper=TRUE)
distance_coords=as.matrix(distance_coords)

phi=seq (60, 80, by=20)
for (j in 1: length(phi))
     {
     covariance=cov.spatial(distance_coords, cov.pars = c(sigma2, phi[j]), cov.model = "sph") 
     #############################################
     # Obtaining curves with spatial correlation
     ############################################
        mantel_test=NULL
        rep=1000
        for (k in 1:rep)
             {
             I=matrix(0,mt, sites)     # Infiltration data matrix
             logI=matrix(0,mt, sites)  # Logaritms of infiltration data matrix
             b=mvrnorm(1, mean, covariance) # Simulating b parameters with spatial correlation

             for (i in 1:sites)
                  {
                   t=seq(0,mt, length=mt)
                   I[,i]=a*t^b[i]
                   logI[,i]=log(a)+b[i]*log(t)
                  }

     ##############################################################
     # Integrating the diferences between infiltartion curves
     # and calculating the matrix of distances between curves
     ##############################################################

             distcurv=matrix(0, nrow=sites, ncol=sites)
             for (i in 1:sites)
                  {
                  for  (j in 1:sites)
                        {
                         dif=(I[,i]-I[,j])^2
                         y=intgnum(dif,t,equispaced=TRUE)
                         distcurv[i,j]=sqrt(y)
                        }
                  }

     #######################################
     # Calculating Mantel's Test
     #######################################

             distance_coords=dist(grid, diag=TRUE, upper=TRUE)
             distance_curves=as.dist(distcurv, diag=TRUE, upper=TRUE)
             mantel_test[k]=mantel.rtest(distance_coords, distance_curves, nrepet=999)$pvalue
            }
      res=table(mantel_test<0.05)
      prob[j]=res[[2]]/rep
     }

cbind(phi,prob)