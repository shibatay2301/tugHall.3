### Pipeline to prepare and execute parallel simulations

# GET Input Dataset -------------------------------------------------------

library( 'tugHall.3' )

# Please, use here your own simulation example
#         to get appropriate values of constant parameters
simulation( verbose = FALSE, to_plot = FALSE,
            seed = 123456 , work_dir = '../Simulation/')

par_exclude = c(  'censor_cells_number', 'censor_time_step', 'clonefile',
                  'cloneoutfile', 'ctmax', 'genefile', 'geneoutfile',
                  'lambda_del', 'lambda_dup', 'logoutfile', 'model_name',
                  'monitor', 'n_repeat', 'real_time_stop',
                  'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
                  'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
                  'tumbler_for_drug_intervention_trial' )

frmt  =  make_input_format( par_exclude = par_exclude )

rng  =  make_input_range( frmt = frmt )

DF   =  make_input_dataset( frmt = frmt, rng = rng, n_simulations = 1000,
                            discrete = TRUE, n_graduations = 11 )

### Save dataset to a file:
write.table( x = DF, file = './Input_DataSet.txt', append = FALSE,
             sep = '\t', row.names = FALSE, col.names = TRUE )

# Make data frame with constant parameters and their vakues:

DF_val  =  lapply( X = par_exclude, FUN = function( x ) pck.env[[ x ]] )

DF_constant  =  data.frame( DF_val, stringsAsFactors = FALSE )
colnames( DF_constant )  =  par_exclude

write.table( x = as.data.frame( t ( DF_constant ) ), file = './Input_const_parameters.txt', append = FALSE,
             sep = '\t', row.names = TRUE, col.names = FALSE )




# Functions for parallel simulation ---------------------------------------





# Make parallel simulations -----------------------------------------------

library( 'parallel' )

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



# Get number of processors:
numCores  =  detectCores()
print( paste0(' Number of processors in the simulation equals  ', numCores ) )

numCores  = round( numCores / 2 )

print( paste0(' Number of using processors in the simulation equals  ', numCores ) )

id_simulations  =  1:nrow( DF )

# Check requared libraries:
check_packages()

mclapply( id_simulations, SIM_PARALLEL, mc.cores = numCores )

# GET Results of simulations ----------------------------------------------





