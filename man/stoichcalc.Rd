\name{stoichcalc-package}
\alias{stoichcalc-package}
\alias{stoichcalc}
\docType{package}
\title{
  R-Functions for Solving Stoichiometric Equations
}
\description{
Given a list of substance compositions, a list of substances
involved in a process, and a list of constraints in addition to mass 
conservation of elementary constituents, the package contains functions 
to build the substance composition matrix, to analyze the uniqueness 
of process stoichiometry, and to calculate stoichiometric coefficients 
if process stoichiometry is unique 
(see reference given below for more details).
}
\details{
\tabular{ll}{
Package: \tab stoichcalc\cr
Type: \tab Package\cr
Version: \tab 1.1-5\cr
Date: \tab 2023-08-28\cr
License: \tab GPL >= 2\cr
LazyLoad: \tab yes\cr
}
The package contains the following three functions:\cr
\code{\link{calc.comp.matrix}} 
constructs the substance composition matrix 
froma a list of substance composition vectors,\cr
\code{\link{calc.stoich.basis}} calculates the basis of the stoichiometry space 
that is compatible with mass balances of elementary constituents 
and additional constraints,
\code{\link{calc.stoich.coef}} calculates the stoichiometric coefficients
of a process from involved substances, their composition and constraints.
}
\author{
Peter Reichert <peter.reichert@emeriti.eawag.ch>
}
\references{
Reichert, P. and Schuwirth, N., 
A generic framework for deriving process stoichiometry in environmental models,
Environmental Modelling and Software 25, 1241-1251, 2010.
}
\keyword{process stoichiometry}
\seealso{
  \code{\link{calc.comp.matrix}},
  \code{\link{calc.stoich.basis}},
  \code{\link{calc.stoich.coef}}
}
\examples{
subst.comp <- 
  list(NH4  = c(H      = 4*1/14,  # gH/gNH4-N
                N      = 1,       # gN/gNH4-N
                charge = 1/14),   # chu/gNH4-N
       NO3  = c(O      = 3*16/14, # gO/gNO3-N
                N      = 1,       # gN/gNO3-N
                charge = -1/14),  # chu/gNO3-N
       HPO4 = c(O      = 4*16/31, # gO/gHPO4-P
                H      = 1*1/31,  # gH/gHPO4-P
                P      = 1,       # gP/gHPO4-P
                charge = -2/31),  # chu/gHPO4-P
       HCO3 = c(C      = 1,       # gC/gHCO3-C
                O      = 3*16/12, # gO/gHCO3-C
                H      = 1*1/12,  # gH/gHCO3-C
                charge = -1/12),  # chu/gHCO3-C
       O2   = c(O      = 1),      # gO/gO2-O
       H    = c(H      = 1,       # gH/molH
                charge = 1),      # chu/molH
       H2O  = c(O      = 1*12,    # gO/molH2O
                H      = 2*1),    # gH/molH2O
       ALG  = c(N      = 0.06,    # gN/gALG
                P      = 0.005,   # gP/gALG
                O      = 0.50,    # gO/gALG
                H      = 0.07,    # gH/gALG
                C      = 0.365),  # gC/gALG
       ZOO  = c(N      = 0.06,    # gN/gZOO
                P      = 0.01,    # gP/gZOO
                O      = 0.50,    # gO/gZOO
                H      = 0.07,    # gH/gZOO
                C      = 0.36),   # gC/gZOO
       POM  = c(N      = 0.04,    # gN/gPOM
                P      = 0.007,   # gP/gPOM
                O      = 0.40,    # gO/gPOM
                H      = 0.07,    # gH/gPOM
                C      = 0.483),  # gC/gPOM
       DOM  = c(N      = 0.04,    # gN/gDOM
                P      = 0.007,   # gP/gDOM
                O      = 0.40,    # gO/gDOM
                H      = 0.07,    # gH/gDOM
                C      = 0.483))  # gC/gDOM
                
Y.ZOO <- 0.2; f.POM <- 0.2; f.DOM <- 0.1

alpha <- calc.comp.matrix(subst.comp)

subst.gro.ALG.NO3 <- c("NO3","HPO4","HCO3",
                       "O2","H","H2O","ALG")

basis.gro.ALG.NO3 <- 
  calc.stoich.basis(alpha,subst.gro.ALG.NO3)

nu.gro.ALG.NO3    <- 
  calc.stoich.coef(alpha       = alpha,
                   name        = "gro.ALG.NO3",
                   subst       = subst.gro.ALG.NO3,
                   subst.norm  = "ALG",
                   nu.norm     = 1)

subst.gro.ZOO <- c("NH4","HPO4","HCO3","O2","H",
                   "H2O","ALG","ZOO","POM","DOM")

basis.gro.ZOO <- 
  calc.stoich.basis(alpha,subst.gro.ZOO)

const.gro.ZOO <- list(c("ZOO" = 1,"ALG" = Y.ZOO),
                      c("POM" = 1,"ALG" = f.POM),
                      c("DOM" = 1,"ALG" = f.DOM)) 

nu.gro.ZOO    <- 
  calc.stoich.coef(alpha       = alpha,
                   name        = "gro.ZOO",
                   subst       = subst.gro.ZOO,
                   subst.norm  = "ZOO",
                   nu.norm     = 1,
                   constraints = const.gro.ZOO)

nu <- rbind(nu.gro.ALG.NO3,
            nu.gro.ZOO)

print(nu,digits=2)
}
