% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tugHall_Model_functions.R
\name{init_onco_clones}
\alias{init_onco_clones}
\title{Function to make list of objects of class 'OncoGene' and generate initial onco settings for all clones (onco_clones)}
\usage{
init_onco_clones(onco1, clones)
}
\arguments{
\item{onco1}{Object of class 'OncoGene'}

\item{clones}{List of objects of class 'Clone'}
}
\value{
List of objects of class 'OncoGene'
}
\description{
Function to make list of objects of class 'OncoGene' and generate initial onco settings for all clones (onco_clones)
}
\examples{
copy_files_to_Input()
load_tugHall.Environment( tugHall_dataset )
clone1 = tugHall_dataset$clones[[ 1 ]]
withr::with_environment( env = pck.env, code = { clones = init_clones(clonefile, clone1)  } )
withr::with_environment( env = pck.env, code = { onco_clones = init_onco_clones( onco1 = onco, clones ) } )
}
