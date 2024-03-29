% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tugHall_Model_functions.R
\name{trial_complex}
\alias{trial_complex}
\title{Function trial for complex case of models}
\usage{
trial_complex(clone1, onco1)
}
\arguments{
\item{clone1}{Object of class 'Clone'}

\item{onco1}{Object of class 'OncoGene'}
}
\value{
Number of new clones originated by clone1
}
\description{
Function trial for complex case of models
}
\examples{
clone1 = tugHall_dataset$clones[[ 1 ]]
onco1 = tugHall_dataset$onco
trial_complex( clone1, onco1 )
unlist( lapply( X = 1:20, FUN = function( x ) trial_complex( clone1, onco1 ) ) )
}
