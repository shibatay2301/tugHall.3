### Pipeline to prepare and execute parallel simulations

# GET Input Dataset -------------------------------------------------------

library( 'tugHall.3' )

# Please, use here your own simulation example
#         to get appropriate values of constant parameters

# Please, just delete the comment symbol to run a simulation
# simulation( verbose = FALSE, to_plot = FALSE, seed = 123456 , work_dir = '../Simulation/')

### OR JUST load dataset as package environment without simulation like:
check_packages()
clear_tugHall.Environment()  # clear package environment before loading
load_tugHall.Environment( results = tugHall_dataset )


# Parameters that should be fixed, so, all the other parameters will be variate:

par_exclude = c(  'censor_cells_number', 'censor_time_step', 'clonefile',
                  'cloneoutfile', 'ctmax', 'genefile', 'geneoutfile',
                  'lambda_del', 'lambda_dup', 'logoutfile', 'model_name',
                  'monitor', 'n_repeat', 'real_time_stop',
                  'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
                  'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
                  'tumbler_for_drug_intervention_trial' )

# That is the list all the parameters:
par_all  =  c( 'Compaction_factor', 'E0', 'F0', 'censor_cells_number',
               'censor_time_step', 'clonefile', 'cloneoutfile', 'd0', 'ctmax',
               'genefile', 'geneoutfile', 'k0',
               'lambda_del', 'lambda_dup', 'logoutfile', 'm0',
               'm_del', 'm_dup', 'model_name', 'monitor',
               'n_repeat', 's0', 'real_time_stop',
               'uo', 'uo_del', 'uo_dup', 'us', 'us_del', 'us_dup',
               'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
               'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
               'tumbler_for_drug_intervention_trial' )

# Parameters to variate in parallel simulations:
par_var  =  par_all[ !( par_all %in% par_exclude ) ]

# Make the format of dataset of input parameters:
frmt  =  make_input_format( par_exclude = par_exclude )

# Make the ranges for all the input parameters:
rng  =  make_input_range( frmt = frmt )
rng$m0  =  c( 1E-12, 1E-8 )
rng$m_del  =  c( 1E-12, 1E-8 )
rng$m_dup  =  c( 1E-12, 1E-8 )

# Make the dataset of input parameters:
DF   =  make_input_dataset( frmt = frmt, rng = rng, n_simulations = 1000,
                            discrete = TRUE, n_graduations = 11 )

### Save dataset to a file:
write.table( x = DF, file = './Input_DataSet.txt', append = FALSE,
             sep = '\t', row.names = FALSE, col.names = TRUE )

# Make data frame with constant parameters and their values:
DF_val  =  lapply( X = par_exclude, FUN = function( x ) pck.env[[ x ]] )

DF_constant  =  data.frame( DF_val, stringsAsFactors = FALSE )
colnames( DF_constant )  =  par_exclude

write.table( x = as.data.frame( t ( DF_constant ) ), file = './Input_const_parameters.txt', append = FALSE,
             sep = '\t', row.names = TRUE, col.names = FALSE )


# Functions for parallel simulation ---------------------------------------

### These functions are not included in the tugHall.3 package to allow
###                the user to modify them according to their requests

#' Function to save input parameters of parallel simulations to the input file as well as
#' gene_hallmarks.txt file with gene-hallmarks weights
#'
#' @description \code{save_to_input()} function copies /Input/ folder to the each folder
#' of parallel trial like /Parallel_simulations/i/Input where i is any integer number or ID of trial.
#' Getting parameters from DF dataset, this function changes related input files like: \cr
#' - parameters.txt \cr
#' - gene_hallmarks.txt \cr
#' - CF.txt
#'
#' @param DF Data frame with variable input parameters
#' @param DF_constant data frame with constant parameters
#' @param i  ID of a simulation to pick up parameters from DF
#' @param main_dir  Name of main folder. Folder to save file is \code{main_dir/i/Input/file_save}
#' @param par_var Vector of names of parameters which will be varied
#' @param file_save Name of file to save. By default \code{file_save = 'parameters.txt'}
#'
#' @return
#'
#' @export
#'
#' @examples
#' NULL
save_to_input  <-  function( DF_constant, DF, i = 1, main_dir, par_var,
                             file_save  =  'parameters.txt' ){

    DF_save  =  cbind( DF_constant, DF[ i, par_var ] )

    if ( !dir.exists( main_dir ) ) dir.create( main_dir )
    if ( !dir.exists( file.path( main_dir, i ) ) ) dir.create( file.path( main_dir, i ) )

    # Copy the template Input from working folder:
    file.copy( from = './Input', to =  file.path( main_dir, i ),
               recursive = TRUE )

    # Save parameters.txt file:
    file_i_save  =  file.path( main_dir, i, 'Input', file_save )

    write.table( x = as.data.frame( t ( DF_save ) ),
                 file = file_i_save, append = FALSE,
                 sep = '\t', row.names = TRUE, col.names = FALSE )

    # Save CF.txt file:
    file_i_save  =  file.path( main_dir, i, 'Input', 'CF.txt' )
    DF_save  =  DF[ i, c( "CompFactor_Ha", "CompFactor_Hd", "CompFactor_Hi",
                          "CompFactor_Hb", "CompFactor_Him"   ) ]
    names( DF_save )  =  c( 'apoptosis', 'growth', 'immortalization', 'angiogenesis', 'invasion' )
    write.table( x = as.data.frame( t ( DF_save ) ),
                 file = file_i_save, append = FALSE,
                 sep = '\t', row.names = TRUE, col.names = FALSE )

    # Save gene_hallmarks.txt file:
    file_i_save  =  file.path( main_dir, i, 'Input', 'gene_hallmarks.txt' )

    hlmrs       =  c( 'Ha', 'Hb', 'Hd', 'Hi', 'Him' )
    hall_names  =  c( 'apoptosis', 'angiogenesis', 'growth',
                     'immortalization', 'invasion' )
    DF_save  =  NULL
    for ( hlmr in hlmrs ){

        hall_name  =  hall_names[ which( hlmr == hlmrs ) ]

        nms  =  pck.env$onco$name[ pck.env$hall[[ hlmr ]] ]

        for( nm in nms ){
            nm_weight  =  paste0( hlmr, '_', nm )
            DF_1  =  data.frame( gene      =  nm,
                                 hallmark  =  hall_name,
                                 onsp      =  pck.env$onco$onsp[ nm == pck.env$onco$name ],
                                 weight    =  DF[ i, nm_weight ]
                                     )
            DF_save  =  rbind( DF_save, DF_1 )

        }
    }

    write.table( x = DF_save,
                 file = file_i_save, append = FALSE,
                 sep = '\t', row.names = FALSE, col.names = FALSE )


}

#' Function to prepare simulation for parallel calculations
#'
#' @description \code{simulation_parallel()} makes a simulation with parameters from Input folder
#' and results save in pck.env as well in './Results_of_simulation.RDS' file in \code{work_dir} folder.
#' Please, PAY ATTENTION user has to change input for all functions like: \cr
#' - define_files_names( );
#' - define_parameters( read_fl = TRUE , file_name = './Input/parameters.txt' );
#' - define_gene_location( genes_list = pck.env$ls_genes );
#' - define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' ).
#'
#' @param verbose Logical type to show or do not show messages during execution
#' @param to_plot Logical type to plot or do not plot graphical results of a simulation
#' @param seed Numeric type to set seed for a simulation, if seed = NA then it will be skipped
#' @param work_dir Working directory for a simulation, by default \code{ work_dir = getwd() }
#' @param copy_input Logical parameter to copy or do not copy default Input folder to the simulation folder
#'
#' @return List of results of simulation with default values for all the parameters
#'
#' @examples
#'
#' # it takes a time for a simulation and then it will demonstrates results, \cr
#' # so, please, wait for a while
#' simulation( verbose = FALSE , to_plot = FALSE )
simulation_parallel  <-  function( verbose = TRUE , to_plot = TRUE,
                          seed = 123456, work_dir = getwd(),
                          copy_input  =  TRUE ){

    local_dir( new = work_dir )

    if ( !is.na( seed ) ) set.seed( seed = seed )

    # if ( verbose ) print('This code will be executed: ')
    # if ( verbose ) print( simulation )

    # Attach packages from import list
    check_packages()

    if ( copy_input ){
        copy_files_to_Input( files = c( 'CCDS.current.txt', 'CF.txt',
                                        'cloneinit.txt','gene_hallmarks.txt',
                                        'gene_map.txt','parameters.txt' ) ,
                             dir = 'Input' )
    }

    check_previous_data( )

    ### BLOCK TO CHANGE BY USER ###############################################

    define_files_names( )
    define_parameters( read_fl = TRUE , file_name = './Input/parameters.txt' )
    define_gene_location( genes_list = pck.env$ls_genes )
    define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' )

    ### FINISH OF BLOCK TO CHANGE BY USER #####################################

    if ( verbose ) print_parameters()

    n_c  =  0
    repeat{
        n_c  =  n_c + 1
        if ( verbose ) print('Start simulation, please, wait for a while ... or see Sim_monitoring.txt file in working folder ')
        smlt = model( )

        if ( file.exists( pck.env$cloneoutfile ) ) break
        if ( n_c  >=  n_repeat )           break
    }
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

    if ( to_plot ){
        plot_order_dysfunction( pck.env$rdr_dysf , pos = c(5,800), logscale = 'y', cex = 1. )

        plot_average_simulation_data( pck.env$data_avg, pck.env$time_max )

        # Main clones
        plot_clone_evolution( pck.env$data_flow, threshold = c(0.01, 1 ), lwd = 2.0,
                              hue = c(" ", "random", "red", "orange", "yellow",
                                      "green", "blue", "purple", "pink", "monochrome")[1],
                              luminosity = c(" ", "random", "light", "bright", "dark")[4],
                              yr = NA , add_initial = TRUE, log_scale = TRUE )

        # Minor clones but large amount of them
        readline('Next? ')

        plot_clone_evolution( pck.env$data_flow, threshold = c(0.0, 0.01), lwd = 2.0,
                              hue = c(" ", "random", "red", "orange", "yellow",
                                      "green", "blue", "purple", "pink", "monochrome")[1],
                              luminosity = c(" ", "random", "light", "bright", "dark")[4],
                              yr = NA , add_initial = FALSE, log_scale = TRUE )

        # Plot VAF for rho = 0
        readline('Next? ')
        plot_VAF( VAF = pck.env$VAF, rho = 0 , violin = FALSE )
    }

    res  =  get_tugHall.Environment()
    saveRDS( object = res, file = './Results_of_simulation.RDS' )

    return( res )
}


# The local function to implement parallel simulations from prepared input folders

### PAY ATTENTION that copy_input  =  TRUE makes all the simulation from default INPUT data
###
### copy_input  =  FALSE  makes simulation from generated data frame,
###                       but it should be REALISTIC values of parameters
###
SIM_PARALLEL <- function( i, copy_input  =  FALSE ){

    fldr  =  file.path( getwd(), 'Parallel_simulations', i )
    if ( !dir.exists(fldr) ) dir.create( fldr )
    res  =  simulation_parallel( verbose = FALSE , to_plot = FALSE, seed = NA,
                        work_dir = fldr,
                        copy_input  =  copy_input )

    # Return VAF from a simulation
    return( res$VAF )
}

# Make parallel simulations -----------------------------------------------

library( 'parallel' )

main_dir  =  file.path( getwd() , 'Parallel_simulations' )

print( 'Please, prepare in the working directory the folder Input/ with all the input files. ' )
print( 'All the files from Input folder will be copied to the input folders for parallel calculations,' )
print( 'and some of them will be modified in accordance with dataset of the input parameters.' )

N_simulations  =  8

for( j in 1:N_simulations ){
    save_to_input( DF_constant = DF_constant, DF = DF, i = j,
                    main_dir = main_dir, par_var = par_var,
                                 file_save  =  'parameters.txt' )
}



# Get number of processors:
numCores  =  detectCores()
print( paste0(' Number of processors in the simulation equals  ', numCores ) )

# Define the convenient number of cores
numCores  = round( numCores / 2 )
print( paste0(' Number of using processors in the simulation equals  ', numCores ) )

id_simulations  =  1:N_simulations

# Check required libraries:
check_packages()

vafs  =  mclapply( id_simulations, SIM_PARALLEL, mc.cores = numCores )

VAF  = lapply( X = 1:length(vafs) , FUN = function( x ){
                                                    VF_1  =  vafs[[ x ]]
                                                    VF_1$ID  =  x
                                                    return( VF_1 )
                                                        } )

### Combine results from the list to a single data frame:
###   TWO methods:

# The first one:
MergeListOfDf = function( data  ){
    if ( length( data ) == 2 ) {
        return( merge( data[[ 1 ]] , data[[ 2 ]], all = TRUE ) )
    }
    return( merge( MergeListOfDf( data[ -1 ] ) , data[[ 1 ]], all = TRUE ) )
}

RES = MergeListOfDf( VAF )


# OR
# The second method:
RES_2  =  Reduce( function(x, y) merge(x, y, all=TRUE), VAF )

### Save data frame to a file:
write.table( x = RES, file = './VAF_parallel.txt', append = FALSE,
             sep = '\t', row.names = FALSE, col.names = TRUE )


# to see all the data of VAF:
plot_VAF( VAF = RES, violin = FALSE, y_lim = c(0,0.55) )

# to see one simulation:
plot_VAF( VAF = VAF[[ 8 ]], violin = FALSE, y_lim = c(0,0.55) )
