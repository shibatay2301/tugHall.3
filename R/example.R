#' Simulation for lazy start with parameters from Input folder
#'
#' @description \code{simulation()} makes a simulation with parameters from Input folder
#' and results save in pck.env as well in './Results_of_simulation.RDS' file in \code{work_dir} folder
#'
#' @param verbose Logical type to show or do not show messages during execution
#' @param to_plot Logical type to plot or do not plot graphical results of a simulation
#' @param seed Numeric type to set seed for a simulation, if seed = NA (by default) then it will be skipped
#' @param work_dir Working directory for a simulation, by default \code{ work_dir = getwd() }
#' @param digits Number of digits in the numeric format of pck.env environment of tugHall
#'
#' @return List of results of simulation with default values for all the parameters
#' @export
#'
#' @examples
#'
#' # it takes a time for a simulation and then it will demonstrates results, \cr
#' # so, please, wait for a while
#' simulation( verbose = FALSE , to_plot = FALSE )
simulation  <-  function( verbose = TRUE , to_plot = TRUE,
                                  seed = NA, work_dir = getwd(), digits = 6 ){

    local_dir( new = work_dir )
    pck.env$digits  =  digits

    if ( !is.na( seed ) ) set.seed( seed = seed )

    # if ( verbose ) print('This code will be executed: ')
    # if ( verbose ) print( simulation )

    # Attach packages from import list
    check_packages()

    copy_files_to_Input( files = c( 'CCDS.current.txt', 'CF.txt',
                                   'cloneinit.txt','gene_hallmarks.txt',
                                   'gene_map.txt','parameters.txt' ) ,
                         dir = 'Input' )

    define_files_names()
    define_gene_location()
    define_parameters( read_fl = TRUE , file_name = './Input/parameters.txt' )
    define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' )
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
    }

    res  =  get_tugHall.Environment()
    saveRDS( object = res, file = './Results_of_simulation.RDS' )

    return( res )
}


#' @describeIn simulation \code{simulation_keep_run()} is needed to start simulation from previous results with new parameter set
#'
#' @description \code{simulation_keep_run()} is needed to start simulation from previous results with new parameter set.
#' Parameter set can be defined as usually from Input folder or keep all the parameters excluding
#' input list of parameters.
#'
#' @param loadRDS logical to load data of previous simulation from file fileRDS.
#' If \code{loadRDS = FALSE} then it loads data from pck.env that should contain the data of a simulation.
#' @param fileRDS file name to load data of previous simulation, only if \code{loadRDS = TRUE}
#' @param change_parameters List of parameters to change from the previous simulation,
#' each parameter should be corresponding to variable name.
#' For example, \code{change_parameters = list(censore_n = 1E06, censore_t = 60 ) }
#'
#' @return List of (clones, onco_clones), where clones - list of objects of class 'Clone', and onco_clones - list of objects of class 'OncoGene'. During a simulation it saves data to geneoutfile.
#' @export
#'
#' @examples
#' NULL
#' \dontrun{
#'
#'
#' }
simulation_keep_run <- function( loadRDS = TRUE,
                                 fileRDS  =  './Results_of_simulation.RDS',
                                 change_parameters = list(),
                                 seed = NA, work_dir = getwd(), digits = 6 ){

    local_environment( env = pck.env )
    local_dir( new = work_dir )

    if ( loadRDS ) {
        res = readRDS( file = fileRDS )
        load_tugHall.Environment( results = res )
    }
    pck.env$digits  =  digits
    if ( !is.na( pck.env$digits ) ){
        defer(options( digits = pck.env$digits ), envir = pck.env )
    }
    if ( !is.na( seed ) ) set.seed( seed = seed )

    # Attach packages from import list
    check_packages()

    write_log(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile,
              E0, F0, m0, uo, us, s0, k0,
              m_dup, m_del, lambda_dup, lambda_del, # CNA parameters
              uo_dup, us_dup, uo_del, us_del,       # CNA parameters
              censore_n, censore_t, d0, Compaction_factor, model_name, time_stop,
              n_repeat, monitor )   # write input parameters

    # Define trial() function: trial_complex or trial_simple
    if ( model_name != 'simplified' ){
        trial  =  trial_complex
    } else {
        trial  =  trial_simple
    }

    write_geneout(geneoutfile, pck.env$hall, Compaction_factor, pck.env$CF)                  # write the geneout.txt file with initial hallmarks
    write_weights("Output/Weights.txt", pck.env$hall)                 # write the weights of genes for hallmarks
    if ( !file.exists( cloneoutfile )){
        write_header( cloneoutfile, pck.env$env, pck.env$onco )                   #
    }
    if ( monitor ) write_monitor( outfile = pck.env$file_monitor,
                                  start = FALSE, clones = clones, env = env )
    cells_number <- sum_N_P_M( pck.env$env, clones )                 # to calculate cells numbers - N,M
    # init_pnt_clones( clones, pck.env$onco_clones )              # initialization of pnt_clones for point mutations

    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    pck.env$hall$updateEnviron( pck.env$env, clones )                     # make averaging for cells
    isFirst = TRUE
    if ( model_name == 'simplified' ) lapply( clones, FUN = function( clone1 ) clone1$invasion = TRUE )
    write_cloneout( cloneoutfile, pck.env$env, clones, isFirst, pck.env$onco_clones )     #  write initial clones

    print( paste0("The probability of an absence of the mutations is p0 = ", as.character( pck.env$onco$p0_1 ) ))
    time_start  =  Sys.time()
    time_current  =  Sys.time()

    clones  =  pck.env$clones

    while(length(clones) > 0 && pck.env$censore_n > cells_number &&
          pck.env$env$T < pck.env$censore_t  &&
          ( as.numeric( difftime( time_current, time_start, units = 'secs') ) < pck.env$time_stop ) ){

        k_old = length(clones)          # the number of clones from last step

        clones_new <- NULL
        onco_clones_new <- NULL

        N_clones_new = unlist( mapply( trial, clones, pck.env$onco_clones ) )

        survived_clones = NULL

        for (i in 1:k_old) {
            if (N_clones_new[i] > 0) {
                for (j in 1:N_clones_new[i])  {
                    clones_new = c(clones_new,clone_copy(clones[[i]]) )
                    onco_clones_new = c(onco_clones_new, onco_copy( pck.env$onco_clones[[i]] ))
                    onco_clones_new[[length(onco_clones_new)]]$id = clones_new[[length(clones_new)]]$id
                }
            }

            # To delete the clones with N_cells == 0 because they are died
            if (clones[[i]]$N_cells == 0 )  survived_clones = c(survived_clones, FALSE)  else survived_clones = c(survived_clones, TRUE )
        }


        # The number of mutations for each NEW clone
        N_new <- length(clones_new)

        if ( N_new > 0) {
            num_mut <- numeric(0)
            sm <- unlist( lapply(onco_clones_new, function(x) (x$sum_prob_1+x$sum_prob_2)/2 ) ) # sum of prob for new clones
            num_mut <- rztpois(N_new, sm ) # Numbers of mutations for each new clone
            # To apply the mutagenesis only to new clones with a number of mutations:
            for ( nn in 1:N_new )  {
                trial_mutagenesis( clones_new[[nn]], num_mut[nn], onco_clones_new[[nn]]  )
            }
        }

        # the new generation = the survived clones + new_clones
        clones = c(clones[survived_clones],clones_new)
        pck.env$onco_clones = c( pck.env$onco_clones[survived_clones], onco_clones_new )

        cells_number <- sum_N_P_M( pck.env$env, clones )                 # to calculate cells numbers - N,M for next step
        lapply(clones,update_Hallmarks)
        pck.env$hall$updateEnviron( pck.env$env, clones )                      # to average probabilities and hallmarks

        pck.env$env$T = pck.env$env$T + 1                                    # to next step

        write_cloneout( cloneoutfile, pck.env$env, clones, isFirst, pck.env$onco_clones )
        #print(c(env$T,env$N,env$M,env$last_id, length(clones), "N_clones_new = ", N_clones_new))
        if ( monitor ) write_monitor( outfile = pck.env$file_monitor, start = FALSE, env = pck.env$env, clones = clones )
        time_current  =  Sys.time()

    }

    # write_pnt_clones( pnt_clones, file = 'Output/point_mutations.txt' )
    # write_pnt_clones( cna_clones, file = 'Output/CNA_mutations.txt' )

    pck.env$clones  =  clones

    return( list( clones , pck.env$onco_clones ) )
}
