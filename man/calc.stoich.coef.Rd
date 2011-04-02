\name{calc.stoich.coef}
\alias{calc.stoich.coef}
\title{Calculate Stoichiometric Coefficients}
\description{
Calculate stoichiometric coefficients of a process from involved substances, 
their composition and constraints
}
\usage{
calc.stoich.coef(alpha, name, subst, subst.norm, nu.norm = 1, constraints = list(), eps = 1e-5)
}
\arguments{
  \item{alpha}{
Substance composition matrix of all substances (labelled columns) 
and mass fractions of elementary constituents (labelled rows). 
Typically calculated by the function \code{\link{calc.comp.matrix}}.}
  \item{name}{
Name of the process
}
  \item{subst}{
Character vector of names of substances affected by the process
(this must be a subset of the column names of alpha)
}
  \item{subst.norm}{
Name of the substance that should have a normalized (given)
stoichiometric coefficient
}
  \item{nu.norm}{
Stoichiometric coefficient of the substance the name of which
is specified in the argument subst.norm
}
  \item{constraints}{
List of stoichiometric constraints in addition to mass conservation of 
elementary constituents. 
Each stoichiometric constraint must be stored as a vector containing the 
coefficients of the linear equation in elementary constituents that defines
the constraint.
The elements of this vector must be labelled by the names of the corresponding
elementary constituents.}
  \item{eps}{
relative tolerance for checking ratios of stoichiometric coefficients 
(only used for informing user about substance pairs with fixed stoichiometric ratio)
}
}
\details{
This is the key function of the package for the calculation of stoichiometric
coefficients of individual processes. 
The results for different processes can easily be bound to the comprehensive 
stoichiometric matrix of all processes by using rbind.
}
\value{
Matrix consisting of one row of stoichiometric coefficients of the process
or an error message if the process stoichiometry is not uniquely defined.
The row name of the matrix is equal to the process name specified as an
argument (to allow binding the stoichiometries of several processes to a
comprehensive stoichiometric matrix), the column names are equal to the
substance names provided by the substance composition matrix alpha.
}
\references{
Reichert, P. and Schuwirth, N.,
A generic framework for deriving process stoichiometry in environmental models,
Environmental Modelling and Software 25, 1241-1251, 2010.
}
\author{
Peter Reichert <peter.reichert@gmx.ch>
}

\seealso{
\code{\link{calc.comp.matrix}}, \code{\link{calc.stoich.basis}}
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
\keyword{process stoichiometry}
