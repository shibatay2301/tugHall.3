% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tugHall_Model_functions.R
\name{write_log}
\alias{write_log}
\title{Function to write log file}
\usage{
write_log(
  genefile,
  clonefile,
  geneoutfile,
  cloneoutfile,
  logoutfile,
  E0,
  F0,
  m0,
  uo,
  us,
  s0,
  k0,
  ctmax,
  m_dup,
  m_del,
  lambda_dup,
  lambda_del,
  uo_dup,
  us_dup,
  uo_del,
  us_del,
  censor_cells_number,
  censor_time_step,
  d0,
  Compaction_factor,
  model_name,
  real_time_stop,
  n_repeat,
  monitor,
  tumbler_for_metastasis_trial,
  tumbler_for_apoptosis_trial,
  tumbler_for_immortalization_trial,
  tumbler_for_angiogenesis_trial,
  tumbler_for_drug_intervention_trial
)
}
\arguments{
\item{genefile}{File name of initial OncoGene information}

\item{clonefile}{File name of info about initial clones}

\item{geneoutfile}{File name for output info about OncoGene information}

\item{cloneoutfile}{File name for output info with clone evolution data}

\item{logoutfile}{Name of log file with all the parameters}

\item{E0}{Parameter in the division probability, numeric type only}

\item{F0}{Parameter in the division probability, numeric type only}

\item{m0}{Mutation probability for point mutation, numeric type only}

\item{uo}{Oncogene mutation probability, numeric type only}

\item{us}{Suppressor mutation probability, numeric type only}

\item{s0}{Parameter in the sigmoid function, numeric type only}

\item{k0}{Environmental death probability, numeric type only}

\item{ctmax}{Hayflick limitation for cell division, integer type}

\item{m_dup}{Mutation probability for duplication, numeric type only}

\item{m_del}{Mutation probability for deletion, numeric type only}

\item{lambda_dup}{CNA duplication average length (of the geometrical distribution for the length), integer type only}

\item{lambda_del}{CNA deletion average length (of the geometrical distribution for the length), integer type only}

\item{uo_dup}{Gene malfunction probability by CNA duplication for oncogene, numeric type only}

\item{us_dup}{Gene malfunction probability by CNA duplication for suppressor, numeric type only}

\item{uo_del}{Gene malfunction probability by CNA deletion    for oncogene, numeric type only}

\item{us_del}{Gene malfunction probability by CNA deletion    for suppressor, numeric type only}

\item{censor_cells_number}{Max cell number where the program forcibly stops, integer type only}

\item{censor_time_step}{Max time where the program forcibly stops, integer type only}

\item{d0}{Initial probability to divide cells, numeric type only}

\item{Compaction_factor}{Compaction factor, logical type only. True means 'to use', False means 'do not use' Compaction factor for hallmarks variables}

\item{model_name}{Name of the model to use. Can be  'proportional_metastatic' or 'threshold_metastatic' or 'simplified'}

\item{real_time_stop}{Max time in seconds of running after that the program forcibly stops, integer type only}

\item{n_repeat}{Max number of repetition of the program until the NON-ZERO output will be getting, integer type only}

\item{monitor}{The indicator to make monitor file during a simulation or do not make, logical type only}

\item{tumbler_for_metastasis_trial}{Logical parameter to turn on/off invasion/metastasis transformation trial}

\item{tumbler_for_apoptosis_trial}{Logical parameter to turn on/off the apoptosis trial}

\item{tumbler_for_immortalization_trial}{Logical parameter to turn on/off the immortalization trial}

\item{tumbler_for_angiogenesis_trial}{Logical parameter to turn on/off angiogenesis trial}

\item{tumbler_for_drug_intervention_trial}{Logical parameter to turn on/off drug intervention trial}
}
\value{
NULL, write log file to Output folder
}
\description{
Function to write log file
}
\examples{
copy_files_to_Input()
define_files_names()
load_tugHall.Environment( tugHall_dataset )
if ( !dir.exists('./Output') ) dir.create('./Output')
\dontrun{
write_log(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile,
E0, F0, m0, uo, us, s0, k0, ctmax, m_dup, m_del, lambda_dup, lambda_del,
uo_dup, us_dup, uo_del, us_del, censor_cells_number, censor_time_step, d0,
Compaction_factor, model_name, real_time_stop, n_repeat, monitor )
}
}
