#################################################################################################################
# Numerical integration 
#
# Input:
#            f  (nx1) valores de la función que se quiere integrar en los puntos de I.
#        xeval  (nx1) Puntos en los que se ha evaluado la función. I(1)=a, I(n)=b.
#   equispaced  TRUE if evalpts are equispaced (used to do faster numerical integration). FALSE is the default. 
# Output:
#   Integ Valor de la integral.
#
# Use:
#     Integ <- intgnum(f,xeval)
#
# (C) Pedro Delicado, Octubre 2001
###################################################################################################################

intgnum <- function(f,xeval,equispaced=FALSE){
   n <- length(f)
   if (equispaced){
      dx <- xeval[2]-xeval[1]
      Integ <- sum(f, na.rm=TRUE)*dx
   } else{
      Integ <- sum( .5*(f[2:n]+f[1:(n-1)])*diff(xeval), na.rm=TRUE )
   }
   Integ
}

