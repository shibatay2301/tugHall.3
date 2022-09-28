## code to prepare `tugHall_dataset` dataset goes here

smlt     =  simulation( seed    =  123456, work_dir =  '../Simulation/',
                                        verbose =  FALSE,  to_plot  =  FALSE)

tugHall_dataset  =  pck.env
# rename files to local working directory:
tugHall_dataset$clonefile     =  './Input/cloneinit.txt'
tugHall_dataset$cloneoutfile  =  './Output/cloneout.txt'

tugHall_dataset$genefile      =  './Input/gene_hallmarks.txt'
tugHall_dataset$geneoutfile   =  './Output/geneout.txt'

tugHall_dataset$logoutfile   =  './Output/log.txt'

tugHall_dataset$mainDir      =  './'



usethis::use_data( tugHall_dataset, overwrite = TRUE )

devtools::install( build_vignettes = TRUE )



