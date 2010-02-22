\name{calc.comp.matrix}
\alias{calc.comp.matrix}
\title{Construct Composition Matrix}
\description{
Construct substance composition matrix 
from list of substance composition vectors
}
\usage{
calc.comp.matrix(subst.comp)
}
\arguments{
  \item{subst.comp}{
Named list of named composition vectors.  
The list must contain entries labelled by the substance names
containing vectors of the mass fractions of elementary 
constituents (typically chemical elements, charge or COD resp.
ThOD) that characterize the composition of the substance. 
Each element of these vectors must be labelled by the name of  
the corresponding elementary constituent.
}
}
\details{
This function compiles the substance composition matrix used in the other
functions of the STOICHCALC library. 
It can alternatively be composed manually or by a user-defined function. 
The main advantage of the use of this function is that substance compositions 
can be maintained in lists.
This makes it much easier to remove and add substances and elementary
constituents.
}
\value{
Composition matrix of all substances (labelled columns) and mass fractions
of elementary constituents (labelled rows).
}
\references{
Reichert, P. and Schuwirth, N.,
A generic framework for deriving process stoichiometry in environmental models,
Environmental Modelling and Software, in review, 2010.
}
\author{
Peter Reichert <peter.reichert@gmx.ch>
}
\seealso{
\code{\link{calc.stoich.basis}}, \code{\link{calc.stoich.coef}}
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
                
alpha <- calc.comp.matrix(subst.comp)

print(alpha)
}
\keyword{process stoichiometry}
