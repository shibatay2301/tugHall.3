#' Environment of the package 'tugHall.3' to store all the objects of a simulation
#'
#' @description \code{pck.env } is environment of the package 'tugHall.3'
#' where all the objects of a simulation are stored and used
#'
#' @export
pck.env = new.env( parent = emptyenv() )
attr( pck.env, "name" ) = "tugHall.Environment"


#' @describeIn pck.env Get results of simulation stored in pck.env or tugHall.Environment environment
#'
#' @description \code{get_tugHall.Environment} function returns all the objects
#' in the pck.env environment of the package tugHall.3
#'
#' @return \code{get_tugHall.Environment} returns all the objects in the pck.env
#' or tugHall.Environment environment
#' @export
#'
#' @examples
#' NULL
get_tugHall.Environment  <-  function(){

    # local_environment( env = pck.env )

    # print( ls( pck.env ) )
    # print('---------------')

    results  =  list()
    l = ls( pck.env )
    if ( length( l ) > 0 ){
        for( i in 1:length( l ) ){
            results[[ l[i] ]]  =  pck.env[[ l[i] ]]
        }
    }

    return( results )

}


#' @describeIn  pck.env Load previous results of simulation to the environment pck.env or tugHall.Environment
#'
#' @description \code{load_tugHall.Environment} loads list
#' 'results' that is results of simulation to the environment pck.env or tugHall.Environment
#'
#' @param results List of results of a simulation to load to the environment pck.env or tugHall.Environment
#'
#' @return \code{load_tugHall.Environment} returns NULL and loads
#' results of simulation to the environment pck.env or tugHall.Environment
#' @export
#'
#' @examples
#' NULL
load_tugHall.Environment  <-  function( results ){

    l = ls( results )
    if ( length( l ) > 0 ){
        for( i in 1:length( l ) ){
            pck.env[[ l[i] ]]  =  results[[ l[i] ]]
        }
    }

    return( NULL )

}


#' @describeIn  pck.env Remove all the objects from the environment pck.env or tugHall.Environment
#'
#' @description \code{clear_tugHall.Environment} clears
#' the environment pck.env or tugHall.Environment
#'
#' @return \code{clear_tugHall.Environment} returns NULL and clears
#' the environment pck.env or tugHall.Environment
#' @export
#'
#' @examples
#' NULL
clear_tugHall.Environment  <-  function( ){

    # l = ls( envir = pck.env )
    # print( l )
    remove( list =  ls( envir = pck.env ), envir = pck.env, inherits = FALSE )

    return( NULL )
}

# Define global variables in tugHall.3:

utils::globalVariables( c( 'Compaction_factor', 'E0', 'F0', 'censor_cells_number',
                           'censor_time_step', 'clonefile', 'cloneoutfile', 'd0', 'ctmax',
                           'gene_map', 'genefile', 'geneoutfile', 'k0',
                           'lambda_del', 'lambda_dup', 'logoutfile', 'm0',
                           'm_del', 'm_dup', 'model_name', 'monitor',
                           'n_repeat', 's0', 'real_time_stop',
                           'uo', 'uo_del', 'uo_dup', 'us', 'us_del', 'us_dup',
                           'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
                           'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
                           'tumbler_for_drug_intervention_trial' ) )







#' Function to read file
#'
#' @param file_name Name of file to read
#' @param stringsAsFactors Parameter for read.table function, by default stringsAsFactors = FALSE
#' @param header Logical type to read or do not read head of a file
#'
#' @return data.frame of data from a file
#' @export
#'
#' @examples
#' fl = system.file('extdata/Input', 'gene_map.txt',package = 'tugHall.3', mustWork = TRUE )
#' read_file(file_name = fl, stringsAsFactors = FALSE )
#' fl = system.file('extdata/Input', 'CF.txt',package = 'tugHall.3', mustWork = TRUE )
#' read_file(file_name = fl, stringsAsFactors = FALSE, header = FALSE )
read_file  <-  function( file_name = '', stringsAsFactors = FALSE, header = TRUE ){
    if ( file.size( file_name )  < 10 ) return( NULL )
    return( read.table( file = file_name, stringsAsFactors  =  stringsAsFactors ,
                        sep="\t", header = header ))
}


#' Function to copy the files by default from extdata folder in the library to Input folder in the working directory
#'
#' @param files Files to copy, vector of names of files by default: \cr
#' files = c( 'CCDS.current.txt', 'CF.txt', 'cloneinit.txt', 'gene_hallmarks.txt','gene_map.txt','parameters.txt' )
#' @param dir Folder to where files should be save, by default dir = 'Input'
#'
#' @return List of logic numbers for each copied file, TRUE - success, FALSE - not success
#' @export
#'
#' @examples
#' files = c('CF.txt', 'cloneinit.txt','gene_hallmarks.txt','gene_map.txt','parameters.txt' )
#' copy_files_to_Input( files, dir = 'Input' )
copy_files_to_Input  <-  function( files = c('CCDS.current.txt', 'CF.txt', 'cloneinit.txt',
                                             'gene_hallmarks.txt','gene_map.txt','parameters.txt' ) ,
                                   dir = 'Input' ){

    fls  =  lapply( X = files,
            FUN = function( x ) system.file('extdata/Input', x, package = 'tugHall.3', mustWork = TRUE ) )
    # fls  =  unlist( fls )

    if ( !file.exists( dir ) ) dir.create( dir )
    lapply( X = 1:length( fls ) , FUN = function( x ){
        file.copy( fls[[ x ]],  dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE )
        } )

}

#' Function to copy the pipelines from extdata folder in the library to /Pipelines/ folder in the working directory
#'
#' @param dir Folder to where files should be save, by default dir = './'
#'
#' @return List of logic numbers for each copied file, TRUE - success, FALSE - not success
#' @export
#'
#' @examples
#' copy_pipelines( dir = 'Input' )
copy_pipelines  <-  function( dir = './' ){

    dir_pck =  system.file('extdata', 'Pipelines', package = 'tugHall.3', mustWork = TRUE )
    files  =  list.files( dir_pck )
    fls  =  lapply( X = files,
                    FUN = function( x ) system.file('extdata/Pipelines', x, package = 'tugHall.3', mustWork = TRUE ) )
    # fls  =  unlist( fls )

    if ( !file.exists( dir ) ) dir.create( dir )
    lapply( X = 1:length( fls ) , FUN = function( x ){
        file.copy( fls[[ x ]],  dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE )
    } )

}

#' Function to copy the files of an example of simulation or from '/extdata/Output/' folder in the library to '/Output/' folder in the working directory
#'
#' @param files Files to copy, vector of names of files by default: \cr
#' files = c( 'cloneout.txt', 'CNA_mutations.txt', 'point_mutations.txt', 'gene_MAP.txt','geneout.txt','log.txt', 'order_genes_dysfunction.txt', 'VAF_data.txt', 'VAF.txt', 'weights.txt' )
#' @param dir Folder to where files should be save, by default dir = 'Output'
#'
#' @return List of logic numbers for each copied file, TRUE - success, FALSE - not success
#' @export
#'
#' @examples
#' files = c( 'cloneout.txt', 'CNA_mutations.txt', 'point_mutations.txt', 'gene_MAP.txt')
#' copy_files_to_Output( files )
#' files = c('geneout.txt','log.txt', 'VAF_data.txt', 'VAF.txt', 'weights.txt' )
#' copy_files_to_Output( files )
copy_files_to_Output  <-  function( files = c( 'cloneout.txt', 'CNA_mutations.txt', 'point_mutations.txt',
                                               'gene_MAP.txt','geneout.txt','log.txt', 'order_genes_dysfunction.txt',
                                               'VAF_data.txt', 'VAF.txt', 'weights.txt' ),
                                   dir = 'Output' ){

    fls  =  lapply( X = files,
                    FUN = function( x ) system.file('extdata/Output', x, package = 'tugHall.3', mustWork = TRUE ) )
    # fls  =  unlist( fls )

    if ( !file.exists( dir ) ) dir.create( dir )
    lapply( X = 1:length( fls ) , FUN = function( x ){
        file.copy( fls[[ x ]],  dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE )
    } )

}


#' Function to make a large number of colors
#'
#' @param nm Number of colors
#'
#' @return Vector of colors with length more than nm
#'
#' @export
#' @examples
#' clrs = gen_colors( nm = 120 )
gen_colors  <-  function(nm = 12){
    # nm is a number of colors
    w <- (nm^(1/3)) %/% 1 +1

    st <-  w^3 %/% nm

    sq <- seq(0,1-1/w,1/w)

    cr <- 1:nm

    l <- 0
    R <- 1:(w^3)
    G <- R
    B <- R

    for (i in 1:w) {
        for (j in 1:w) {
            for (k in 1:w) {
                l <- l+1
                R[l] <- sq[i]
                G[l] <- sq[j]
                B[l] <- sq[k]
            }
        }
    }

    # seq(1,w^3,st) # is consequence of each color to make a high diversity of colors
    jColor <- data.frame( number = 1:length( seq( 1,w^3, st ) ),
                          color  = rgb( R[seq( 1, w^3, st ) ], G[seq( 1, w^3, st)],
                                        B[seq( 1, w^3, st ) ] ), stringsAsFactors = FALSE )

    return(jColor)

}


#' Check the installation of a package for some functions
#'
#' @param pkg Package name
#'
#' @return if the package is installed then it returns NULL else it returns error message
#'
#' @export
#' @examples
#' check_pkg( pkg = 'grDevices' )
check_pkg  <-  function( pkg ){
    msg  =  paste0( 'Package ', pkg, ' must be installed to use this function. \n ' )
    if ( !requireNamespace( pkg , quietly = TRUE ) )    stop( msg, call. = FALSE )
}


#' Check the installation of packages and attach them with corresponding functions
#'
#' @param pkgs List of package names with related function names,
#' by default (or when pkgs = NULL) the list of packages are described in Namespace file of the package or
#' 'R/MaxWiK-package.R' file
#'
#' @return if the packages are installed then it returns NULL else it returns error message
#'
#' @export
#'
#' @examples
#' check_packages(  )
check_packages  <-  function( pkgs = NULL ){

    if ( is.null( pkgs ) ) {
        pkgs  =  list(  actuar = 'rztpois',
                        dplyr  = '%>%',
                        ggplot2  =  c( 'aes', 'geom_boxplot', 'geom_point', 'geom_violin',
                                       'ggplot', 'ggtitle', 'theme', 'xlab', 'ylab', 'ylim' ),
                        graphics = c('axis', 'legend', 'lines', 'par', 'plot', 'text', 'title' ),
                        grDevices = c('dev.off', 'pdf', 'rgb'),
                        methods = 'new',
                        randomcoloR = 'randomColor',
                        stats = c('aggregate', 'rbinom', 'rexp', 'rnorm', 'runif' ),
                        stringr = c('str_length', 'str_split', 'str_sub', 'str_trim', 'str_remove', 'str_replace_all', 'str_c', 'str_match', 'str_trunc' ),
                        utils = c('read.delim', 'read.table', 'write.table', 'globalVariables' ),
                        withr  =  c('local_environment', 'local_par', 'local_dir', 'local_options', 'defer' )
                        )
    }

    ### Attach the packages
    for( pck in names( pkgs ) ){
        check_pkg( pkg = pck )
        require( package = pck, character.only = TRUE, include.only = pkgs[[ pck ]])
    }
}


#' Foolproof function allows to checking the consistency of all the input parameters
#'
#' @description Foolproof function allows to checking the consistency of all the input parameters.
#' It should be used just before a loop of simulation. So, all the parameters should be defined,
#' objects onco and hall should be initialized as well. The function checks the list of parameters' names,
#' absence of NA and NULL in the input data, self-consistency of genes names,
#' correctness of hallmarks values, and finally, that all the necessary information is defined.
#'
#' @return NULL
#'
#' @export
#'
#' @examples
#' NULL
foolproof  <-  function(){

    ### Have to check:
    # 1. List of genes in different objects: onco, hall, gene_map
    #2. list of variables - whole or particular
    #3. NA and NULL check in onco and hall

    # 1. List of genes in different objects: onco, hall, gene_map
    genes_onco  =  sort( unique( pck.env$onco$name) )
    genes_hall  =  sort( unique( c( pck.env$hall$Ha, pck.env$hall$Hi, pck.env$hall$Hd,
                                        pck.env$hall$Hb, pck.env$hall$Him )  ) )
    genes_map   =  sort( unique( pck.env$gene_map$Gene ) )

    if ( length( genes_onco ) == length( genes_map ) & length( genes_onco ) == length( genes_hall ) ){
        if ( !all.equal( genes_onco, genes_map ) ){
            stop( 'Names of genes defined in onco and gene_map are different.' )
        }

        if ( ! length( genes_onco) == length( genes_hall ) ){
            stop( 'Not all the genes use in hallmrks definition.' )
        }
    } else {
        stop( 'Names of genes are different in onco, hall and gene_map.' )
    }

    #3. NA and NULL check in onco and hall:

    ### Check onco object:
    flds  =  c( "name", "cds_1", "cds_2", "rna_1", "rna_2", "p0_1", "p0_2", "prob_1", "prob_2",
                "sum_prob_1", "sum_prob_2", "onsp", "len" )

    for( v in flds ){
        if ( length(  pck.env$onco[[ v ]] ) == 0 |
             any( is.na(   pck.env$onco[[ v ]] ) )   |
             is.null( pck.env$onco[[ v ]] ) ){
            stop( paste0( 'The parameter ', v, ' in the onco object is not defined.'))
        }
    }

    ### Check hall object:
    flds  =  c( "Ha", "Hi", "Hd", "Hb", "Him", "Ha_w", "Hi_w", "Hd_w", "Hb_w", "Him_w" )

    for( v in flds ){
        if ( length(  pck.env$hall[[ v ]] ) == 0 |
             any( is.na(   pck.env$hall[[ v ]] ) )    |
             is.null( pck.env$hall[[ v ]] ) ){
            stop( paste0( 'The parameter ', v, ' in the hall object is not defined.'))
        }
    }

    #2. list of variables - whole or particular
    env_names  =  c( 'Compaction_factor', 'E0', 'F0', 'censor_cells_number',
                     'censor_time_step', 'clonefile', 'cloneoutfile', 'd0', 'ctmax',
                     'gene_map', 'genefile', 'geneoutfile', 'k0',
                     'lambda_del', 'lambda_dup', 'logoutfile', 'm0',
                     'm_del', 'm_dup', 'model_name', 'monitor',
                     'n_repeat', 's0', 'real_time_stop',
                     'uo', 'uo_del', 'uo_dup', 'us', 'us_del', 'us_dup',
                     'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
                     'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
                     'tumbler_for_drug_intervention_trial' )

    # Check that all the parameters are defined and have values
    for( v in env_names ){
        if ( length( pck.env[[ v ]] ) == 0 | any( is.na( pck.env[[ v ]] ) ) | is.null( pck.env[[ v ]] ) ){
            stop( paste0( 'The parameter ', v, ' is not defined. '))
        }
    }

    return( NULL )
}
