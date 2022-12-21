
library('tugHall.3')

set.seed( seed = 1 )
verbose  =  TRUE
# if ( verbose ) print('This code will be executed: ')
# if ( verbose ) print( simulation )

# Attach packages from import list
check_packages()

############## INITIAL Simulation


#copy_files_to_Input( files = c( 'CCDS.current.txt', 'CF.txt',
#                                'cloneinit.txt','gene_hallmarks.txt',
#                                'gene_map.txt','parameters.txt' ) ,
#                                 dir = 'Input' )

check_previous_data( )
define_files_names()
define_gene_location()
define_parameters( read_fl = TRUE , file_name = './Input/parameters.txt' )
define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' )
if ( verbose ) print_parameters(  )

smlt = model( )

# clones        =  pck.env$clones       # smlt[[ 1 ]]
# onco_clones   =  pck.env$onco_clones  # smlt[[ 2 ]]

if ( verbose ) print('Save point mutations and CNA to the files in Output folder. ')
write_pnt_clones( pck.env$pnt_clones, file_out = 'Output/point_mutations.txt' )
write_pnt_clones( pck.env$cna_clones, file_out = 'Output/CNA_mutations.txt' )

if ( verbose ) print('Get data of the last simulation from cloneout file. ')
dtst = get_flow_data( pck.env$cloneoutfile, pck.env$genefile )
pck.env$data_avg   =  dtst$data_avg
pck.env$data_flow  =  dtst$data_flow
pck.env$time_max   =  dtst$time_max
pck.env$data_last  =  dtst$data_last
cna_mut = dtst$cna_mut
pnt_mut   =  dtst$pnt_mut
pnt_mut_B =  dtst$pnt_mut_B
# onco = dtst$onco
# hall = dtst$hall



if ( verbose ) print('Also get VAF data and save them into file Output/VAF_data.txt ' )
pck.env$vf = get_VAF( pnt_mut, pck.env$data_last )
pck.env$VAF  =  get_rho_VAF( vf = pck.env$vf, rho = c( 0.0, 0.1, 0.2, 0.5, 0.7, 0.9 ) ,
                             file_name = './Output/VAF.txt' )

pck.env$rdr_dysf  =  get_order_of_genes_dysfunction( pnt_mut = pnt_mut,
                                                     pck.env$data_last, cna_mut,
                                                     file_name = './Output/order_genes_dysfunction.txt' )


plot_order_dysfunction( pck.env$rdr_dysf , pos = c(5,800), logscale = 'y', cex = 1. )

plot_average_simulation_data( pck.env$data_avg, pck.env$time_max )


plot_clone_evolution( pck.env$data_flow, threshold = c(0.01, 1 ), lwd = 2.0,
                          hue = c(" ", "random", "red", "orange", "yellow",
                                  "green", "blue", "purple", "pink", "monochrome")[1],
                          luminosity = c(" ", "random", "light", "bright", "dark")[4],
                          yr = NA , add_initial = TRUE, log_scale = TRUE )


plot_clone_evolution( pck.env$data_flow, threshold = c(0.0, 0.01), lwd = 2.0,
                          hue = c(" ", "random", "red", "orange", "yellow",
                                  "green", "blue", "purple", "pink", "monochrome")[1],
                          luminosity = c(" ", "random", "light", "bright", "dark")[4],
                          yr = NA , add_initial = FALSE, log_scale = TRUE )


res  =  get_tugHall.Environment()
vf_init  =  pck.env$vf

plot_VAF( VAF = pck.env$VAF, rho = 0 , violin = FALSE )

# Split original vf for drivers and passengers:
vf_drivers  =  vf_init[ which(  vf_init$MalfunctionedByPointMut ), ]
vf_pass     =  vf_init[ which( !vf_init$MalfunctionedByPointMut ), ]
VAF_drivers  =  get_rho_VAF( vf = vf_drivers, save_to_file = FALSE )
VAF_pass  =  get_rho_VAF( vf = vf_pass, save_to_file = FALSE )

if (length(VAF_drivers) > 0 ) plot_VAF( VAF = VAF_drivers, rho = 0 , violin = FALSE )
if (length(VAF_pass) > 0 )    plot_VAF( VAF = VAF_pass,    rho = 0 , violin = FALSE )


# Cancer care steps:
CC_steps    =  7   # Define number of steps for cancer care:
kill_prob   =  0
block_prob  =  0.5
gene        =  'PIK3CA'

for ( st in pck.env$time_max : (pck.env$time_max + CC_steps ) ){
    drug_intervention( kill_prob = kill_prob, block_prob = block_prob, gene = gene )
    restart_simulation( loadRDS = FALSE, fileRDS = '', loadInput = FALSE,
                        change_parameters = list(censor_cells_number = 1E06,
                                                 censor_time_step = st + 1 ),
                        seed = NA, work_dir = './', digits = 6, to_plot = FALSE, verbose = FALSE )

}



if ( verbose ) print('Get data of the last simulation from cloneout file. ')
dtst = get_flow_data( pck.env$cloneoutfile, pck.env$genefile )
pck.env$data_avg   =  dtst$data_avg
pck.env$data_flow  =  dtst$data_flow
pck.env$time_max   =  dtst$time_max
pck.env$data_last  =  dtst$data_last
cna_mut = dtst$cna_mut
pnt_mut   =  dtst$pnt_mut
pnt_mut_B =  dtst$pnt_mut_B
# onco = dtst$onco
# hall = dtst$hall



if ( verbose ) print('Also get VAF data and save them into file Output/VAF_data.txt ' )
pck.env$vf = get_VAF( pnt_mut, pck.env$data_last )
pck.env$VAF  =  get_rho_VAF( vf = pck.env$vf, rho = c( 0.0, 0.1, 0.2, 0.5, 0.7, 0.9 ) ,
                             file_name = './Output/VAF.txt' )

pck.env$rdr_dysf  =  get_order_of_genes_dysfunction( pnt_mut = pnt_mut,
                                                     pck.env$data_last, cna_mut,
                                                     file_name = './Output/order_genes_dysfunction.txt' )


plot_order_dysfunction( pck.env$rdr_dysf , pos = c(5,800), logscale = 'y', cex = 1. )

plot_average_simulation_data( pck.env$data_avg, pck.env$time_max )


plot_clone_evolution( pck.env$data_flow, threshold = c(0.01, 1 ), lwd = 2.0,
                      hue = c(" ", "random", "red", "orange", "yellow",
                              "green", "blue", "purple", "pink", "monochrome")[1],
                      luminosity = c(" ", "random", "light", "bright", "dark")[4],
                      yr = NA , add_initial = TRUE, log_scale = TRUE )


plot_clone_evolution( pck.env$data_flow, threshold = c(0.0, 0.01), lwd = 2.0,
                      hue = c(" ", "random", "red", "orange", "yellow",
                              "green", "blue", "purple", "pink", "monochrome")[1],
                      luminosity = c(" ", "random", "light", "bright", "dark")[4],
                      yr = NA , add_initial = FALSE, log_scale = TRUE )






vf_drug  =  get_VAF_clones( env = pck.env$env, clones = pck.env$clones,
                            pnt_clones = pck.env$pnt_clones )
VAF_drug  =  get_rho_VAF( vf = vf_drug, save_to_file = FALSE )
plot_VAF( VAF = VAF_drug, rho = 0, violin = FALSE )


# check the blocking drivers (moved to passengers):
vf_drivers  =  vf_drug[ which(  vf_drug$MalfunctionedByPointMut ), ]
vf_pass     =  vf_drug[ which( !vf_drug$MalfunctionedByPointMut ), ]
VAF_drivers  =  get_rho_VAF( vf = vf_drivers, save_to_file = FALSE )
VAF_pass  =  get_rho_VAF( vf = vf_pass, save_to_file = FALSE )

if (length(VAF_drivers) > 0 ) plot_VAF( VAF = VAF_drivers, rho = 0 , violin = FALSE )
if (length(VAF_pass) > 0 ) plot_VAF( VAF = VAF_pass, rho = 0 , violin = FALSE )








