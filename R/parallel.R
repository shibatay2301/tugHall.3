### Code to prepare parallel simulations

#' Function to prepare a format of dataset of input parameters for parallel calculations
#'
#' @param par_exclude List of parameters to exclude from data frame of
#' input parameters because they will be constant for all the simulations
#'
#' @description \code{make_input_format()} function allows to prepare a format of dataset of input parameters from
#' results of a trial simulation. By default it equals all the names of files with some additional arguments: \cr
#' \code{par_exclude = c(  'censor_cells_number', 'censor_time_step', 'clonefile', 'cloneoutfile', 'ctmax',
#' 'genefile', 'geneoutfile', 'lambda_del', 'lambda_dup', 'logoutfile', 'model_name', 'monitor',
#' 'n_repeat', 'real_time_stop', 'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
#' 'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial', 'tumbler_for_drug_intervention_trial' ) }
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



#' Function to prepare dataset of input parameters for parallel calculations
#'
#' @return \code{make_input_dataset()} returns data frame with different sets of input parameters
#'
#' @export
#'
#' @examples
#' NULL
make_input_dataset  <-  function( ){

    if ( length( pck.env ) < 30 ){
        stop( 'Please, use this function after trial simulation to get a format of data.' )
    }

}