###   MODEL: main function 'model' --------------------------------------------

#' Main function 'model' to simulate clones' evolution
#'
#' @param genefile Name of file with input data of hallmarks
#' @param clonefile Name of file with initial clones data
#' @param geneoutfile Name of file with hallmark data
#' @param cloneoutfile Name of file with cloneout data
#' @param logoutfile File name of log file
#' @param E0 Parameter in the division probability, numeric type only
#' @param F0 Parameter in the division probability, numeric type only
#' @param m0 Mutation probability for point mutation, numeric type only
#' @param uo Oncogene mutation probability, numeric type only
#' @param us Suppressor mutation probability, numeric type only
#' @param s0 Parameter in the sigmoid function, numeric type only
#' @param k0 Environmental death probability, numeric type only
#' @param d0 Initial probability to divide cells, numeric type only
#' @param censore_n Max cell number where the program forcibly stops, integer type only
#' @param censore_t Max time where the program forcibly stops, integer type only
#'
#' @return List of (clones, onco_clones), where clones - list of objects of class 'Clone', and onco_clones - list of objects of class 'OncoGene'. During a simulation it saves data to geneoutfile.
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' define_files_names()
#' define_gene_location()
#' define_parameters( read_fl = TRUE , file_name = './Input/parameters.txt' )
#' define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' )
#' time_stop = 3  #  Duration of simulation time is 3 sec
#' \dontrun{
#' res = model( genefile, clonefile, geneoutfile, cloneoutfile, logoutfile,
#' E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0 )
#' }
model <- function(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile,
                  E0, F0, m0, uo, us, s0, k0, censore_n, censore_t, d0) {
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

    onco = oncogene$new()        # make the vector onco about the hallmarks
    onco$read(genefile)          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall = hallmark$new()        # make a vector hall with hallmarks parameters
    hall$read( genefile, onco$name, normalization = TRUE )     # read from the genefile - 'gene_hallmarks.txt'
    env = environ$new(F0)               # new vector for average values of cells
    pnt = Point_Mutations$new()
    pnt_clones = NULL
    cna = CNA_Mutations$new()
    cna_clones = NULL
    mut_order  <<-  0            #  mutation order to reproduce gene map
    assign("mut_order", mut_order, envir = .GlobalEnv)
    assign("onco", onco, envir = .GlobalEnv)
    assign("env", env, envir = .GlobalEnv)
    assign("pnt", pnt, envir = .GlobalEnv)
    assign("pnt_clones", pnt_clones, envir = .GlobalEnv)
    assign("cna", cna, envir = .GlobalEnv)
    assign("cna_clones", cna_clones, envir = .GlobalEnv)
    assign("hall", hall, envir = .GlobalEnv)
    assign("uo", uo, envir = .GlobalEnv)
    assign("us", us, envir = .GlobalEnv)
    clone1 = clone$new(gene_size=length(onco$cds_1),
                       m=m0, s=s0, k=k0, E=E0)          # clone1  -  empty object of clone
    clones = init_clones(clonefile, clone1)           # clones - the clones with hallmarks from cellfile - cellinit.txt - initial cells
    onco_clones = init_onco_clones( onco, clones )    # onco_clones - the onco related to each clone in clones
    write_geneout(geneoutfile, hall, Compaction_factor, CF)                  # write the geneout.txt file with initial hallmarks
    write_weights("Output/Weights.txt", hall)                 # write the weights of genes for hallmarks
    write_header( cloneoutfile, env, onco )                   #
    if ( monitor ) write_monitor( start = TRUE )
    cells_number <- sum_N_P_M(env, clones)                 # to calculate cells numbers - N,M
    init_pnt_clones( clones, onco_clones )              # initialization of pnt_clones for point mutations

    lapply(clones,update_Hallmarks)                     # to calculate the Hallmarks and probabilities for initial cells
    hall$updateEnviron(env, clones)                     # make averaging for cells
    isFirst = TRUE
    if ( model_name == 'simplified' ) lapply( clones, FUN = function( clone1 ) clone1$invasion = TRUE )
    write_cloneout( cloneoutfile, env, clones, isFirst, onco_clones )     #  write initial clones

    print( paste0("The probability of an absence of the mutations is p0 = ", as.character(onco$p0_1) ))
    time_start  =  Sys.time()
    time_current  =  Sys.time()
    while(length(clones) > 0 && censore_n > cells_number &&
          env$T < censore_t  &&
          ( as.numeric( time_current - time_start ) < time_stop ) ){

        k_old = length(clones)          # the number of clones from last step

        clones_new <- NULL
        onco_clones_new <- NULL

        N_clones_new = unlist( mapply( trial, clones, onco_clones ) )

        survived_clones = NULL

        for (i in 1:k_old) {
            if (N_clones_new[i] > 0) {
                for (j in 1:N_clones_new[i])  {
                    clones_new = c(clones_new,clone_copy(clones[[i]]) )
                    onco_clones_new = c(onco_clones_new, onco_copy(onco_clones[[i]]))
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
        onco_clones = c(onco_clones[survived_clones],onco_clones_new)

        cells_number <- sum_N_P_M(env, clones)                 # to calculate cells numbers - N,M for next step
        lapply(clones,update_Hallmarks)
        hall$updateEnviron(env, clones)                      # to average probabilities and hallmarks

        env$T = env$T + 1                                    # to next step

        write_cloneout( cloneoutfile, env, clones, isFirst, onco_clones )
        #print(c(env$T,env$N,env$M,env$last_id, length(clones), "N_clones_new = ", N_clones_new))
        if ( monitor ) write_monitor( start = FALSE, env = env, clones = clones )
        time_current  =  Sys.time()

    }

    # write_pnt_clones( pnt_clones, file = 'Output/point_mutations.txt' )
    # write_pnt_clones( cna_clones, file = 'Output/CNA_mutations.txt' )

    return( list( clones , onco_clones ) )
}


