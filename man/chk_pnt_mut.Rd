% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tugHall_Model_functions.R
\name{chk_pnt_mut}
\alias{chk_pnt_mut}
\title{Function to check point mutations match or don't match into duplication or deletion}
\usage{
chk_pnt_mut(pnt1, Ref_start, Ref_end, Chr, prntl)
}
\arguments{
\item{pnt1}{Object of class 'Point_Mutations'}

\item{Ref_start}{Initial position of CNA}

\item{Ref_end}{Final position of CNA}

\item{Chr}{Chromosome name}

\item{prntl}{Parental chromosome 1 or 2}
}
\value{
Logical: TRUE if point mutation matches CNA, FALSE if it doesn't match
}
\description{
Function to check point mutations match or don't match into duplication or deletion
}
\examples{
pnt1 = tugHall_dataset$pnt_clones[[ 5 ]]
pstn = pnt1$Phys_pos[1]
message( pstn )
prntl = pnt1$Parental_1or2
Chr = pnt1$Chr
chk_pnt_mut( pnt1 , Ref_start = pstn - 200, Ref_end = pstn + 200, Chr, prntl )
chk_pnt_mut( pnt1 , Ref_start = pstn - 200, Ref_end = pstn - 100, Chr, prntl )
}
