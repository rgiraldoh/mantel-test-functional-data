#---------------------------------------------------------
# R code to carry out the Mantel test for functional
# data based on infiltration data recorded at an experimental
# in Puerto Lopez, Colombia
# ---------------------------------------------------------

#---------------------------
# R libraries requiered
#---------------------------

library(ade4)
library(geoR)
library(MASS)
source("integral.R")            # Allows to calculate an integral by using Trapecio´s rule

#-------------------------------------------------------------------------
# Data: Coeficcients of Kostiakov models (y = at^b) fitted to infiltration 
# data recorded at 40 sites of an experimental of 1.12 ha located 
# in Puerto López, Meta, Colombia
#--------------------------------------------------------------------------

a=c(17.621, 18.380, 4.312, 11.561, 30.009,
      4.829, 18.583, 13.192, 54.117, 4.429,
      0.570, 4.833, 0.549, 9.224, 12.535, 
      0.201, 0.763, 1.183, 8.583, 8.764,
      0.502, 0.461, 0.494, 0.924, 0.933, 
      0.171, 0.251, 0.323, 0.510, 0.287, 
      1.395, 0.987, 0.432, 0.345, 1.755, 
      3.882, 19.653, 0.208, 1.618, 0.876)

b=c(0.684, 0.799, 0.743, 0.717, 0.670,
     0.730, 0.780, 0.698, 0.466, 0.438,
     0.128, 0.776, 0.265, 0.919, 0.729,
     0.055, 0.064, 0.380, 0.857, 0.495,
     0.022, 0.898, 0.093, 0.643, 0.124,
     0.209, 0.178, 0.360, 0.267, 0.410,
     0.404, 0.278, 0.267, 0.419, 0.326,
     0.491, 0.515, 0.424, 0.528, 0.455)


#----------------------------------------------------------
# Grid of sites
#----------------------------------------------------------
grid<-expand.grid(seq(5,85, length=5),seq(35,175,length=8))
plot(grid, xlab="x", ylab="y", type="p", pch=19)
t=seq(0,180, 1)


#----------------------------------------------------------
# Kostiakov models
#----------------------------------------------------------

Kos=matrix(0, nrow=length(t), ncol=length(a1))
for (j in 1: length(a1))
     {
     Kos[,j]=a[j]*t^(b[j])
     }
    
matplot(Kos, type="l", lty=1, xlab="Time (min)", ylab="Infiltration (m)", col=1)


#----------------------------------------------------------
# Distance between curves obtaine with the Kostiakov model
#----------------------------------------------------------

sites=length(a)
distcurv=matrix(0, nrow=sites, ncol=sites)
for (i in 1:sites)
     {
      for  (j in 1:sites)
            {
             dif=(Kos[,i]-Kos[,j])^2
             y=intgnum(dif,t,equispaced=TRUE)
             distcurv[i,j]=sqrt(y)
             }
      }

#---------------------------------------------------------
# Distance between sites
#----------------------------------------------------------

distance_coords=dist(grid, diag=TRUE, upper=TRUE)
distance_curves=as.dist(distcurv)



#-----------------------------------------------------------------------------------
# Mantel test by using the library ade4
#------------------------------------------------------------------------------------

mantel_test=mantel.rtest(distance_coords, distance_curves, nrepet=9999)
plot(mantel_test, xlab="M statistic",  ylab ="Frequency", main="")
hist(mantel_test$sim, xlab="M statistic", xlim=c(0, 0.3), ylab ="Frequency", main="")
abline(v=mantel_test$obs)



