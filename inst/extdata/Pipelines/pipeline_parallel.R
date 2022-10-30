### Pipeline to prepare and execute parallel simulations

# GET Input Dataset -------------------------------------------------------

library( 'tugHall.3' )

# Please, use here your own simulation example
#         to get appropriate values of constant parameters

# Please, just delete the comment symbol to run a simulation
# simulation( verbose = FALSE, to_plot = FALSE, seed = 123456 , work_dir = '../Simulation/')

### OR JUST load dataset as package environment without simulation like:

clear_tugHall.Environment()  # clear package environment before loading
load_tugHall.Environment( results = tugHall_dataset )

par_exclude = c(  'censor_cells_number', 'censor_time_step', 'clonefile',
                  'cloneoutfile', 'ctmax', 'genefile', 'geneoutfile',
                  'lambda_del', 'lambda_dup', 'logoutfile', 'model_name',
                  'monitor', 'n_repeat', 'real_time_stop',
                  'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
                  'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
                  'tumbler_for_drug_intervention_trial' )

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

frmt  =  make_input_format( par_exclude = par_exclude )

rng  =  make_input_range( frmt = frmt )

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

# The local function to implement parallel simulations from prepared input folders

SIM_PARALLEL <- function( i ){

    res  =  simulation( verbose = FALSE , to_plot = FALSE, seed = NA,
                        work_dir = file.path( getwd(), 'Parallel_simulations', i ),
                        copy_input  =  TRUE )

    # Return VAF from a simulation
    return( res$VAF )
}

# Make parallel simulations -----------------------------------------------

library( 'parallel' )

main_dir  =  file.path( getwd() , 'Parallel_simulations' )

print( 'Please, prepare in the working directory the folder Input/ with all the input files. ' )
print( 'All the files from Input folder will be copied to the input folders for parallel calculations,' )
print( 'and some of them will be modified in accordance with dataset of the input parameters.' )

N_simulations  =  16

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












# THE OLD FUNCTIONS -------------------------------------------------------

# The local function to implement parallel simulations

SIM_PARALLEL <- function( i ){

                Initialization( i, name_model = name_model, name_weights = name_weights, name_init = name_init )

                Get_parameters( file = file_param, i , name_weights = name_weights, name_init = name_init, name_model = name_model )

                ### Get_initial_clones( clonefile = clonefile )   ### This function to generate the file for initial clones

                VAF  <-  SIM( i )

                if ( !is.null( VAF ) ) {
                    df   <-  Safe_SIM_MAX( VAF, i , onco , num_of_max = 5, file = paste0( ResultDir, 'data_sim_all_', floor( i / 1000 ), '.txt' ) )
                } else {
                    df <- NULL
                }

                if ( !is.null( VAF ) ) {
                    df   <-  Safe_SIM_MAX_Primary( VAF, i , onco , num_of_max = 5, file = paste0( ResultDir, 'data_sim_primary_', floor( i / 1000 ), '.txt' ) )
                } else {
                    df <- NULL
                }
                ### REMOVE the files and folder after calculations
                unlink( paste0( 'Output/', name_weights, "/", name_model, "/", name_init, "/", i ) , recursive = TRUE, force = TRUE )

    return( df )
}



# GET Results of simulations ----------------------------------------------





