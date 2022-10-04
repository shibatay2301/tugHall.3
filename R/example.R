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


#' @describeIn simulation \code{restart_simulation()} is needed to start simulation from previous results with new parameter set
#'
#' @description \code{restart_simulation()} is needed to start simulation from previous results with new parameter set.
#' Parameter set can be defined as usually from Input folder or keep all the parameters excluding
#' input list of parameters.
#'
#' @param loadRDS logical to load data of previous simulation from file fileRDS.
#' If \code{loadRDS = FALSE} then it loads data from pck.env that should contain the data of a simulation.
#' @param fileRDS file name to load data of previous simulation, only if \code{loadRDS = TRUE}
#' @param change_parameters List of parameters to change from the previous simulation,
#' each parameter should be corresponding to variable name.
#' For example, \code{change_parameters = list(censor_cells_number = 1E06, censor_time_step = 60 ) }
#' @param loadInput Logical to load parameters from Input folder or not.
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
restart_simulation <- function( loadRDS = TRUE,
                                 fileRDS = './Results_of_simulation.RDS',
                                 loadInput = FALSE,
                                 change_parameters = list(censor_cells_number = 1E06, censor_time_step = 60 ),
                                 seed = NA, work_dir = getwd(), digits = 6,
                                 to_plot = TRUE, verbose = FALSE ){

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

    # Load changed parameters:
    if ( length(change_parameters) > 0 ){
        for( i in 1:length(change_parameters) ){
            nm  =  names( change_parameters)[ i ]
            pck.env[[ nm ]]  =  change_parameters[[ i ]]

        }
    }

    # Attach packages from import list
    check_packages()

    # Simulate from new parameters
    model_keep_run()

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
    saveRDS( object = res, file = './Results_of_restart_simulation.RDS' )

    return( res )
}
