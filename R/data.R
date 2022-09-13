
#' tugHall dataset named 'tugHall_dataset'
#'
#' Dataset contains all the necessary data.frames and objects to check functions of tugHall.
#' Description of each data.frame and object could be found in documentation to tugHall package.
#'
#' @format A data frame with 15 data.frames/lists and 33 objects:
#' \describe{
#'   \item{Input parameters}{Input parameters: censor_cells_number, censor_time_step, Compaction_factor, \cr
#'   d0, E0, F0, k0, m0, lambda_del, lambda_dup, m_del, m_dup, model_name, monitor, n_repeat, \cr
#'   s0, real_time_stop, uo_del, uo_dup, us_del, us_dup, uo, us }
#'   \item{CF}{Data frame of compaction factor }
#'   \item{Names of files and folder}{Names of files to input and output data: clonefile, cloneoutfile, \cr
#'   file_monitor, genefile, geneoutfile, logoutfile, mainDir }
#'   \item{data_flow}{simulation data for all time steps, data from file cloneout.txt}
#'   \item{data_last}{simulation data for the last time step, data from file cloneout.txt}
#'   \item{data_avg}{simulation data averaged for the each time step, data from file cloneout.txt}
#'   \item{pnt_clones}{list of all the point mutations}
#'   \item{cna_clones}{list of all the CNA mutations}
#'   \item{clones}{list of all the clones}
#'   \item{env}{list of average data for the last timestep (environment of clones)}
#'   \item{gene_map}{data.frame with genes' locations information}
#'   \item{hall}{Object of class 'HallMark'}
#'   \item{onco}{Object of class 'OncoGene'}
#'   \item{time_max}{Value of maximal time step in an example simulation}
#'   \item{mut_order}{Value of integer indicator of current mutation order in the simulation}
#'   \item{vf}{data.frame of preliminary data for VAF calculations}
#'   \item{VAF}{data.frame with VAF values for different rho}
#'   \item{rdr_dysf}{data.frame of order of genes dysfunction for each clone}
#' }
#'
#'
'tugHall_dataset'

