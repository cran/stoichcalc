# ==============================================================================
# STOICHCALC - R-Routines for Solving Stoichiometric Equations
# ==============================================================================
#
# Version 1.1-0                Peter Reichert, Feb. 08, 2010 (reichert@eawag.ch)
# =============
#
# Literature: Peter Reichert and Nele Schuwirth
#             A generic framework for deriving process stoichiometry 
#             in environmental models
#             submitted to Environmental Modelling and Software
#
# ==============================================================================


# calc.comp.matrix
# ================

# Description
# -----------

# Construct composition matrix from list of substance composition vectors

# Usage
# -----

# calc.comp.matrix(subst.comp)

# Arguments
# ---------

# subst.comp   named list of named composition vectors.  
#              The list must contain entries labelled by the substance names
#              containing vectors of the "mass" fractions of "elementary 
#              constituents" (typically chemical elements, charge or COD resp.
#              ThOD) that characterize the composition of the substance. 
#              Each element of these vectors must be labelled by the name of  
#              the corresponding "elementary constituent".

# Details
# -------

# This function compiles the substance composition matrix used in the other
# functions of the STOICHCALC library. It can alternatively be composed
# manually or by a user-defined function. The main advantage of the use of
# this function is that substance compositions can be maintained in lists.
# This makes it much easier to remove and add substances and "elementary
# constituents".

# Value
# -----

# Composition matrix of all substances (labelled columns) and "mass" fractions
# of "elementary constituents" (labelled rows)

calc.comp.matrix <- function(subst.comp)
{
   # check input:
   
   if ( length(subst.comp) == 0 )
   {
      print("Unable to construct composition matrix:")
      print("no substance composition vectors provided")
      return(NA);
   }
   elements <- character(0)
   for ( i in 1:length(subst.comp) )
   {
      elements <- c(elements,names(subst.comp[[i]]))
   }
   elements <- unique(elements)
   if ( length(elements) == 0 )
   {
      print("Unable to construct composition matrix:")
      print("no elementary constituents found")
      return(NA);
   }
   
   # set up composition matrix:
   
   alpha <- matrix(data=0,nrow=length(elements),ncol=length(subst.comp))
   colnames(alpha) <- names(subst.comp)
   rownames(alpha) <- elements
   
   # fill in elements of composition matrix:
   
   for ( i in 1:length(subst.comp) )
   {
      for ( k in 1:length(subst.comp[[i]]) )
      {
         alpha[names(subst.comp[[i]])[k],i] <- subst.comp[[i]][k]
      }
   }
   
   # print summary and return composition matrix:
   
   print(paste("Composition matrix (",nrow(alpha),"x",ncol(alpha),
               ") successfully constructed:",sep=""))
   print(paste("el. const.:",paste(rownames(alpha),collapse=",")))
   print(paste("substances:",paste(colnames(alpha),collapse=",")))
   return(alpha)
}


# ==============================================================================


# calc.stoich.basis
# =================

# Description
# -----------

# Calculate the basis of the stoichiometry space that is compatible with
# mass balances of "elementary constituents" and additional constraints

# Usage
# -----

# calc.stoich.basis(alpha,subst=NA,constraints=list())

# Arguments
# ---------

# alpha        substance composition matrix of all substances (labelled columns) 
#              and "mass" fractions of "elementary constituents" (labelled 
#              rows). Typically calculated by the function "calc.comp.matrix". 
# subst        character vector of names of substances to be used for analysis
#              (this must be a subset of the column names of alpha)
# constraints  list of stoichiometric constraints in addition to "mass"
#              conservation of "elementary constituents". Each stoichiometric
#              constraint must be stored as a vector containing the coefficients
#              of the linear equation in "elementary constituents" that defines
#              the constraint. The elements of this vector must be labelled
#              by the names of the corresponding "elementary constituents".

# Details
# -------

# This function is primarily used in the function "calc.stoich.coef".
# However, it can also be used to check the number of required stoichiometric
# constraints in addition to "mass" conservation of "elementary constituents"
# for a given process. In this case the composition matrix should only contain
# the substances relevant for this process. The number of required constraints 
# is then equal to the row dimension of the output matrix minus 1.

# Value
# -----

# Matrix of basis vectors (in rows) that span the compatible stoichiometric 
# space.


calc.stoich.basis <- function(alpha,subst=NA,constraints=list())
{
   # calculate reduced composition matrix:
   
   alpha.reduced <- alpha
   if ( !is.na(subst[1]) )
   {
      substances <- unique(subst)
      if ( sum(ifelse(is.na(match(substances,colnames(alpha))),1,0)) > 0 )
      {
         print(paste("Substance name vector contains name",
                     "that does not correspond to column name of alpha"))
         return(NA)
      }
      alpha.reduced <- alpha[,substances]
      if ( is.vector(alpha.reduced) ) 
      {
         names(alpha.reduced) <- substances
         alpha.reduced <- t(as.matrix(alpha.reduced))
      }
      else
      {
         colnames(alpha.reduced) <- substances
      }
   }

   # transpose composition matrix:
   
   a <- t(alpha.reduced)
   
   # add constraints as additional columns to t(alpha):
   
   if ( length(constraints) > 0 )
   {
      for ( i in 1:length(constraints) )
      {
         a <- cbind(a,rep(0,nrow(a)))
         for ( j in 1:length(constraints[[i]]) )
         {
            a[names(constraints[[i]])[j],ncol(a)] <- constraints[[i]][j]
         }
      }
   }
   
   # extend t(alpha) by additional columns of zeros if necessary:
   
   if ( nrow(a) > ncol(a) )
   {
      a <- cbind(a,matrix(0,nrow=nrow(a),ncol=nrow(a)-ncol(a)))
   }
   
   # perform singular value decomposition:
   
   svd.res <- svd(a)
   
   # extract basis of null space:
   
   ut <- t(svd.res$u)
   d  <- svd.res$d
   d.max <- max(d)
   fact  <- 1e-10
   res <- matrix(nrow=0,ncol=ncol(ut))
   colnames(res) <- colnames(alpha.reduced)
   for ( i in 1:length(d) )
   {
      if ( d[i] < fact*d.max )
      {
         res <- rbind(res,ut[i,])
      }
   }
   
   # extend to original dimension and order of composition matrix:
   
   basis <- matrix(0,nrow=nrow(res),ncol=ncol(alpha))
   colnames(basis) <- colnames(alpha)
   basis[,colnames(res)] <- res 
   
   # print summary and return extracted basis:

   print(paste("Number of substances:             ",ncol(alpha.reduced)))
   print(paste("Number of elementary constituents:",nrow(alpha.reduced)))
   print(paste("Number of constraints:            ",length(constraints)))   
   print(paste("Number of independent processes:  ",nrow(basis)))
   return(basis)
}


# ==============================================================================


# calc.stoich.coef
# ================

# Description
# -----------

# Calculate stoichiometry of a process from involved substances and constraints

# Usage
# -----

# calc.stoich.coef(alpha,name,subst,subst.norm,nu.norm=1,constraints=list())

# Arguments
# ---------

# alpha        substance composition matrix of all substances (labelled columns) 
#              and "mass" fractions of "elementary constituents" (labelled 
#              rows). Typically calculated by the function "calc.comp.matrix".
# name         name of the process
# subst        character vector of names of substances affected by the process
#              (this must be a subset of the column names of alpha)
# subst.norm   name of the substance that should have a normalized (given)
#              stoichiometric coefficient
# nu.norm      stoichiometric coefficient of the substance the name of which
#              is specified in the argument "subst.norm" 
# constraints  list of stoichiometric constraints in addition to "mass"
#              conservation of "elementary constituents". Each stoichiometric
#              constraint must be stored as a vector containing the coefficients
#              of the linear equation in "elementary constituents" that defines
#              the constraint. The elements of this vector must be labelled
#              by the names of the corresponding "elementary constituents".

# Details
# -------

# This is the key function of the library for the calculation of stoichiometric
# coefficients of individual processes. The results for different processes can
# easily be bound to the comprehensive stoichiometric matrix of all processes
# by using "rbind".

# Value
# -----

# Matrix consisting of one row of stoichiometric coefficients of the process
# or an error message if the process stoichiometry is not uniquely defined.
# The row name of the matrix is equal to the process name specified as an
# argument (to allow binding the stoichiometries of several processes to a
# comprehensive stoichiometric matrix), the column names are equal to the
# substance names provided by the substance composition matrix alpha.


calc.stoich.coef <- function(alpha,name,subst,subst.norm,nu.norm=1,
                             constraints=list())
{
   # calculate basis of stoichiometric coefficients:
   
   nu.basis <- calc.stoich.basis(alpha,c(subst.norm,subst),constraints)
   if ( is.na(nu.basis[1]) )
   {
      print("No process possible")
      return(NA)
   }
   
   # check if one dimensional space:
   
   if ( nrow(nu.basis) == 0 )
   {
      print("No process possible.")
      return(NA)
   }
   if ( nrow(nu.basis) >  1 ) 
   {
      if ( length(constraints) == 0 )
      {
         print(paste("Process not unique,",nrow(nu.basis)-1,
                     "constraints needed."))
      }
      else
      {
         print(paste("Process not unique,",nrow(nu.basis)-1,
                     "additional constraints needed."))
      }
      return(NA)
   }
   
   # calculate stoichiometry vector as one row matrix:
   
   nu <- matrix(data=0,nrow=1,ncol=ncol(alpha))
   colnames(nu) <- colnames(alpha)
   nu[1,colnames(nu.basis)] <- nu.basis[1,]/nu.basis[1,subst.norm]*nu.norm
   rownames(nu) <- name
   
   # print and return stoichiometry:
   
   print("Stoichiometric coefficients successfully calculated.")
   return(nu)
}


# ==============================================================================
