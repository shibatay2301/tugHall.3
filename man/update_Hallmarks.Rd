% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tugHall_Classes.R
\name{update_Hallmarks}
\alias{update_Hallmarks}
\title{Function to update Hallmark and variable after division or under initialization}
\usage{
update_Hallmarks(clone1)
}
\arguments{
\item{clone1}{Object of class 'Clone'}
}
\value{
The same object of class 'Clone' with updated fields
}
\description{
Function to update Hallmark and variable after division or under initialization
}
\examples{
clone = tugHall_dataset$clones[[ 1 ]]
load_tugHall.Environment( tugHall_dataset )
withr::with_environment( env = pck.env, code = update_Hallmarks( clone ) )
}
