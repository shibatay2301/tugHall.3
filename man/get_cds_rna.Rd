% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tugHall_Mutation_functions.R
\name{get_cds_rna}
\alias{get_cds_rna}
\title{Function to get length of CDS and of genes from data.frame gene_map and related probabilities}
\usage{
get_cds_rna(gm)
}
\arguments{
\item{gm}{data.frame gene_map with info about genes' location}
}
\value{
list( names, CDS, RNA, PROB, SUM, P0 )
}
\description{
Function to get length of CDS and of genes from data.frame gene_map and related probabilities
}
\examples{
gene_map  =  tugHall_dataset$gene_map
load_tugHall.Environment( tugHall_dataset )
withr::with_environment( env = pck.env, code = get_cds_rna( gm = gene_map ) )
}
