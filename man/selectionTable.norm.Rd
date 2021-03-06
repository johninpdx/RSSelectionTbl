% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectionTables.r
\name{selectionTable.norm}
\alias{selectionTable.norm}
\title{Calculates the location of the social norm and its standard error
in the attraction function}
\usage{
selectionTable.norm(x, xd, name, vname)
}
\arguments{
\item{x}{sienaFit: Results from a single group analysis}

\item{xd}{siena: the RSiena  Data Object or Data Set that
was used to generate x}

\item{name}{character: the name of the endogenous variable
of interest (a network or behavior variable)}

\item{vname}{character: the actor variable name (i.e. predictor
variable)}
}
\value{
list: [1] The calculated norm, [2] the SE of the norm
}
\description{
Calculates the location of the social norm and its standard error
in the attraction function
}
\details{
Assumes an attraction function with terms
'altX','altSqX','egoX','egoSqX','diffSqX'. Note that the attraction function
 should not include the 'egoXaltX' effect. See Snijders & Lomi(c2018),
 "Beyond Homophily", for more info.
}
