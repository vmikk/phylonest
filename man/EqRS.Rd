% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EqRS.R
\name{EqRS}
\alias{EqRS}
\title{First decomposition of diversity introduced in Pavoine et al. 2016.}
\usage{
EqRS(
  df,
  dis = NULL,
  structures = NULL,
  option = c("eq", "normed1", "normed2"),
  formula = c("QE", "EDI"),
  tol = 1e-08
)
}
\arguments{
\item{df}{Dataframe or matrix with sites as rows and species as columns. Entries are abundances of species within sites.}

\item{dis}{Dissimilarity among species (NULL or class 'dist').}

\item{structures}{Data frame that contains the name of the group (row) of an level (column) to which the site belongs. Sites in structures should be in the same order as in df. Default is NULL.}

\item{option}{Rescaling type ("eq", "normed1" or "normed2").}

\item{formula}{Quadratic entropy formula ("QE", "EDI"). "QE" is default.}

\item{tol}{A tolerance threshold (a value less than tol is considered equal to zero).}
}
\value{
A data frame with each component of the selected diversity decomposition.

Rescaled estimates of the beta components (option = "normed1" or "normed2") reach the maximum value of 1 when sites within a region are maximally dissimilar whatever the level of diversity within sites (e.g., within any region, sites do not share species and any species from any site always is maximally dissimilar from all species in all other sites); beta components should be equal 0 when sites are identical within regions.
}
\description{
Use when the interest is in multiple-site and multiple-region beta diversity, the sampling design is even (same number of sites within regions), sites within a region can be given equal weights and regions can be given equal weights.
}
\details{
Rescaling types (argument 'option'):
\itemize{
\item \strong{"eq"} - the diversity components are given in terms of equivalent number of species, sites, regions etc.;
\item \strong{"normed1"} - the normed components of diversity will be returned with formula (1 – 1 / E) / (1 - 1 / Emax) [de Bello et al., 2010, Eqn 14];
\item \strong{"normed2"} - the normed components of diversity will be returned with formula (E – 1) / (Emax - 1) [Villéger et al., 2012].
}
For Eβ, Emax = M (the number of sites). For Eα and Eγ, Emax=S (the number of species in the data set).

Formula type (argument 'formula'):
\itemize{
\item \strong{"QE"} - the definition of the quadratic entropy is following Rao (1982);
\item \strong{"EDI"} - the Euclidean Diversity Index of Champely and Chessel (2002).
}

For the associated permutation test see \code{\link{randtestEqRS}}.
}
\examples{
data(macroloire)
# Taxonomic dissimilarities among species:
dTaxo <- dist.taxo(macroloire$taxo)^2/2
dTaxo <- dTaxo/max(dTaxo)
# Size-based dissimilarities among species
dSize <- dist.prop(macroloire$traits[ ,1:4], method = 2)
# Dissimilarities among species in terms of feeding categories
dFeed <- dist.prop(macroloire$traits[ ,5:11], method = 2)
# Dissimilarities among species in terms of both size and feeding categories
dSF <- (dSize+dFeed)/2

# Table with sites as rows (stations), species as columns and abundances as entries
ab <- as.data.frame(t(macroloire$fau))
# Table with sites as rows and one column only. Entries indicate the geological region associated with each site
stru <- macroloire$envir["Morphoregion"]

EqRS(ab, , stru, option="eq")
EqRS(ab, dTaxo, stru, formula = "QE", option="eq")
EqRS(ab, dSize, stru, formula = "QE", option="eq")
EqRS(ab, dFeed, stru, formula = "QE", option="eq")
EqRS(ab, dSF, stru, formula = "QE", option="eq")

EqRS(ab, , stru, option="normed2")
EqRS(ab, dTaxo, stru, formula = "QE", option="normed2")
EqRS(ab, dSize, stru, formula = "QE", option="normed2")
EqRS(ab, dFeed, stru, formula = "QE", option="normed2")
EqRS(ab, dSF, stru, formula = "QE", option="normed2")

}
\references{
Pavoine S., Marcon E., Ricotta C. (2016) ‘Equivalent numbers’ for species, phylogenetic or functional diversity in a nested hierarchy of multiple scales. Methods Ecol Evol, 7: 1152-1163. DOI:10.1111/2041-210X.12591

de Bello F., Lavergne S., Meynard C.N., Lepš J., Thuiller W. (2010) The partitioning of diversity: showing Theseus a way out of the labyrinth. Journal of Vegetation Science, 21: 992-1000. DOI: 10.1111/j.1654-1103.2010.01195.x

Villéger S., Miranda J.R., Hernandez D.F., Mouillot, D. (2012) Low functional beta-diversity despite high taxonomic beta-diversity among tropical estuarine fish communities. PLoS One, 7, e40679. DOI: 10.1371/journal.pone.0040679

Champely S., Chessel D. (2002) Measuring biological diversity using Euclidean metrics. Environmental and Ecological Statistics, 9, 167-177. DOI: 10.1023/A:1015170104476

Rao C.R. (1982) Diversity and dissimilarity coefficients: a unified approach. Theoretical Population Biology, 21, 24-43. DOI: 10.1016/0040-5809(82)90004-1
}
\seealso{
\code{\link{randtestEqRS}}, \code{\link{EqRSintra}}, \code{\link{EqRao}}.
}
\author{
Sandrine Pavoine, Eric Marcon, Carlo Ricotta.
}
