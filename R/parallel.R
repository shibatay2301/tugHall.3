### Code to prepare parallel simulations


# Get input dataset -------------------------------------------------------



#' @describeIn make_input_dataset Function to prepare a format of dataset of input parameters for parallel calculations
#'
#' @param par_exclude List of parameters to exclude from data frame of
#' input parameters because they will be constant for all the simulations
#'
#' @description \code{make_input_format()} function allows to prepare a format of dataset of input parameters from
#' results of a trial simulation.
#'
#' @return \code{make_input_format()} returns data frame with a single row corresponding to
#' a set of current input parameters
#'
#' @export
#'
#' @examples
#' NULL
make_input_format  <-  function( par_exclude = c(  'censor_cells_number', 'censor_time_step', 'clonefile', 'cloneoutfile', 'ctmax',
                                                   'genefile', 'geneoutfile', 'lambda_del', 'lambda_dup', 'logoutfile', 'model_name', 'monitor',
                                                   'n_repeat', 'real_time_stop', 'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
                                                   'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial', 'tumbler_for_drug_intervention_trial' )  ){

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

    if ( !all( par_exclude  %in% par_all ) ) stop( 'Not all the parameters in the list. Check the list.' )

    if ( length( par_exclude ) > 0 ){
        wch  =  sapply( par_exclude, FUN = function( x ) which( x == par_all ) )
        prmtrs  =  par_all[ - wch ]
    } else {
        prmtrs  =  par_all
    }

    if ( length( pck.env ) < 30 ){
        stop( 'Please, use this function after trial simulation to get a format of data.' )
    }

    val_prmtrs  =  lapply( prmtrs, FUN = function( x )   pck.env[[ x ]] )

    DF_par  =  data.frame( val_prmtrs, stringsAsFactors = FALSE )
    colnames( DF_par ) = prmtrs

    # Block related to hallmarks-genes weights
    genes  =  pck.env$onco$name
    hallm  =  c( 'Ha', 'Hi', 'Him', 'Hb', 'Hd' )

    Ha     =  pck.env$hall$Ha_w
    names( Ha )  =  str_c(  'Ha_', genes[ pck.env$hall$Ha ] )

    Hi     =  pck.env$hall$Hi_w
    names( Hi )  =  str_c(  'Hi_', genes[ pck.env$hall$Hi ] )

    Him     =  pck.env$hall$Him_w
    names( Him )  =  str_c( 'Him_', genes[ pck.env$hall$Him ] )

    Hb     =  pck.env$hall$Hb_w
    names( Hb )  =  str_c(  'Hb_', genes[ pck.env$hall$Hb ] )

    Hd     =  pck.env$hall$Hd_w
    names( Hd )  =  str_c(  'Hd_', genes[ pck.env$hall$Hd ] )

    DF_hall  =  t( data.frame( c ( Ha, Hb, Hd, Hi, Him ),
                            row.names = names( c ( Ha, Hb, Hd, Hi, Him ) ) ) )
    rownames( DF_hall)  =  1

    DF_CF  =  data.frame( pck.env$CF )
    colnames( DF_CF )  =  str_c(  'CompFactor_', colnames( DF_CF ) )

    DF  =  cbind( DF_hall, DF_CF, DF_par )

    return( list( format = DF, Const_parameters = par_exclude ) )
}



#' @describeIn make_input_dataset Function to make the range for each input parameter in the data frame
#'
#' @param frmt List of results of function \code{make_input_format()} as input format for the range of each parameter
#'
#' @return \code{make_input_range()} returns a data frame with two rows, the first row is minimal values,
#' and the second row is maximal values of parameters.
#'
#' @export
#'
#' @examples
#' NULL
make_input_range  <-  function( frmt ){

    rng  =  frmt$format
    rng[ 2, ] = rng[ 1, ]
    tps  =  sapply( frmt$format, class)

    for( i in 1:ncol( rng ) ){
        if ( tps[ i ] == 'numeric' ){
            rng[ 1, i ] = 0
            rng[ 2, i ] = 1
        } else {
            if ( tps[ i ] == 'logical' ){
                rng[ 1, i ] = FALSE
                rng[ 2, i ] = TRUE
            } else {
                if (tps[ i ] == 'character'  ){
                    rng[ 1, i ] = NA
                    rng[ 2, i ] = NA
                    print( 'NA in the range data frame for string parameter. Please, change it in the future.' )
                } else stop( paste0( 'Type of variable ',  names( rng)[ i ], ' is not define in the function.' ) )
            }
        }
    }

    return( rng )
}

#' Function to prepare dataset of input parameters for parallel calculations
#'
#' @param rng Data frame was gotten as a result of function \code{make_input_range()}
#' @param n_simulations Number of rows for output data frame corresponding to a number of simulations.
#' @param discrete Logical parameter, if TRUE then random values will be generated from discrete set of values,
#' if FALSE then random values will be generated from continuous range.
#' @param n_graduations Number of discrete values for parameter generation.
#' Applicable only if discrete is TRUE.
#'
#' @return \code{make_input_dataset()} returns data frame with different sets of input parameters
#'
#' @export
#'
#' @examples
#' NULL
make_input_dataset  <-  function( frmt, rng, n_simulations = 10,
                                  discrete = TRUE, n_graduations = 11 ){

    if ( n_graduations < 2 ) stop( 'n_graduations sholud be more than 1. ' )
    DF  = frmt$format
    DF[ 1:n_simulations, ]  =  DF[ 1, ]

    for( j in 1:ncol( DF ) ){

        if ( class( rng[ , j ]) == 'numeric' ){

            mn  =  rng[ 1, j ]
            mx  =  rng[ 2, j ]

            if ( discrete ){
                st  =  ( mx - mn ) / ( n_graduations - 1 )
                DF[ 1:n_simulations, j ]  =  sample( ( 1:n_graduations - 1 ), size = n_simulations,
                                                     replace = TRUE )  *  st
            } else {
                DF[ 1:n_simulations, j ]  =  runif( n_simulations, min = mn, max = mx )
            }
        } else {
            if ( class( rng[ , j ]) == 'logical' ){
                DF[ 1:n_simulations, j ]  =  sample( c( FALSE, TRUE ), size = n_simulations,
                                                     replace = TRUE )
            } else {
                print( paste0( 'The values for the column ', names( DF )[ j ], ' is not generated, so, please, do it by hands.' ) )
            }
        }
    }

    return( DF )
}







# Run parallel simulations ------------------------------------------------


