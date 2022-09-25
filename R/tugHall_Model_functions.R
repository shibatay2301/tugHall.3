

# Functions related to trial  ---------------------------------------------



#' Aggregate data of a clone for environment object
#'
#' @param env Object of class 'Environ'
#' @param clones List of all the objects of class 'Clone'
#'
#' @return NULL, but global variable env is updated
#' @export
#'
#' @examples
#' clones = tugHall_dataset$clones
#' env = tugHall_dataset$env
#' sum_cell(env, clones)
#' message( paste0('Number of primary tumor cells in the pool of clones is ', env$P ) )
#' message( paste0('Number of normal cells in the pool of clones is ', env$N ) )
#' message( paste0('Number of metastatic cells in the pool of clones is ', env$M ) )
sum_cell <- function(env, clones) {
    if (length(clones) > 0) {
        avg = apply(matrix(unlist(lapply(clones, sum_mutation)),ncol=length(clones)),1,sum)  #  /length(clones)
        sNMP  =  pck.env$env$N + pck.env$env$M + pck.env$env$P
        pck.env$env$c = avg[1] / sNMP
        pck.env$env$d = avg[2] / sNMP
        pck.env$env$i = avg[3] / sNMP
        pck.env$env$a = avg[4] / sNMP
        pck.env$env$k = avg[5] / sNMP
        pck.env$env$E = avg[6] / sNMP
        pck.env$env$Nmax = avg[7] / sNMP
        pck.env$env$im = avg[8] / sNMP
        pck.env$env$Ha = avg[9] / sNMP
        pck.env$env$Him = avg[10] / sNMP
        pck.env$env$Hi = avg[11] / sNMP
        pck.env$env$Hb = avg[12] / sNMP
        pck.env$env$Hd = avg[13] / sNMP
        pck.env$env$type = avg[14] / sNMP
        pck.env$env$mutden = avg[15] / sNMP
    } else {
        pck.env$env$M = 0
        pck.env$env$N = 0
        pck.env$env$P = 0
        pck.env$env$c = 0
        pck.env$env$d = 0
        pck.env$env$i = 0
        pck.env$env$a = 0
        pck.env$env$k = 0
        pck.env$env$E = 0
        pck.env$env$Nmax = 0
        pck.env$env$im = 0
        pck.env$env$Ha = 0
        pck.env$env$Him = 0
        pck.env$env$Hi = 0
        pck.env$env$Hb = 0
        pck.env$env$Hd = 0
        pck.env$env$type = 0
    }
    # env$posdriver = rep("", length(onco$name))
    # env$pospasngr = rep("", length(onco$name))
}

#' Serve function for sum_cell() function
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return vector of clone1 variables to aggregate in sum_cell() function
#' @export
#'
#' @examples
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' sum_mutation(clone1)
sum_mutation <- function(clone1) {
    return(c(clone1$c*clone1$N_cells,      clone1$d*clone1$N_cells,    clone1$i*clone1$N_cells,
             clone1$a*clone1$N_cells,      clone1$k*clone1$N_cells,    clone1$E*clone1$N_cells,
             clone1$Nmax*clone1$N_cells,   clone1$im*clone1$N_cells,   clone1$Ha*clone1$N_cells,
             clone1$Him*clone1$N_cells,    clone1$Hi*clone1$N_cells,   clone1$Hb*clone1$N_cells,
             clone1$Hd*clone1$N_cells,     ifelse(clone1$invasion,1,0),
             clone1$mutden*clone1$N_cells)) # ,
    # clone1$N_cells * ifelse(clone1$invasion,0,1),   clone1$N_cells * ifelse(clone1$invasion,1,0) ))

    #           clone1$gene*clone1$N_cells) )
}

#' Function to calculate N and M numbers - normal and metastatic cells
#'
#' @param env Object of class 'Environ'
#' @param clones List of all the objects of class 'Clone'
#'
#' @return Number of all the cells in a simulation (normal + primary tumor + metastatic)
#' @export
#'
#' @examples
#' clones = tugHall_dataset$clones
#' env = tugHall_dataset$env
#' env$M = 0
#' env$P = 0
#' env$N = 0  # View( env )
#' sum_N_P_M(env, clones)  # View( env )
#' message( paste(env$N, env$P, env$M ) )
sum_N_P_M <- function(env, clones) {
    if (length(clones) > 0) {
        avg = apply(matrix(unlist(lapply(clones, number_N_P_M)),ncol=length(clones)),1,sum)  #  /length(clones)
        pck.env$env$N = avg[ 1 ]
        pck.env$env$P = avg[ 2 ]
        pck.env$env$M = avg[ 3 ]
        return(pck.env$env$N + pck.env$env$P + pck.env$env$M)
    }
}

#' Function to get number of cells of a clone with indicator (normal, primary tumor or metastatic)
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return Vector c( N_normal, N_primary, N_metastatic )
#' @export
#'
#' @examples
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' number_N_P_M(clone1)
#' message( paste('Format is as follow: ', 'N_normal', 'N_primary', 'N_metastatic' ) )
number_N_P_M <- function(clone1) {
    indicator = c(  ifelse(clone1$invasion,0,1) * ifelse(clone1$primary,0,1),
                    ifelse(clone1$invasion,0,1) * ifelse(clone1$primary,1,0),
                    ifelse(clone1$invasion,1,0) )
    return( clone1$N_cells * indicator )
}


# trial functions for complex and simplified cases:
# The result is a number of new cells, if N_New < 0 it means that the number is decreased.

#' Function trial for complex case of models
#'
#' @param clone1 Object of class 'Clone'
#' @param onco1 Object of class 'OncoGene'
#'
#' @return Number of new clones originated by clone1
#' @export
#'
#' @examples
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' onco1 = tugHall_dataset$onco
#' trial_complex( clone1, onco1 )
#' unlist( lapply( X = 1:20, FUN = function( x ) trial_complex( clone1, onco1 ) ) )
trial_complex <- function( clone1, onco1 ) {

    # trial for Environmental death of cell
    N_die  =  calc_binom( 1, clone1$N_cells, clone1$k )   # The number of cells to die due to the Environmental death of cells in clone

    # Apoptosis trial
    if ( pck.env$tumbler_for_apoptosis_trial ){
        N_die  =  N_die + calc_binom( 1, clone1$N_cells, clone1$a )
    }

    # invasion / metastasis trial
    if (clone1$im > 0 & pck.env$tumbler_for_metastasis_trial ) {
        if (!clone1$invasion) {
            N_die  =  N_die + calc_binom( 1, clone1$N_cells, ( 1 - clone1$im ) )
            if ( model_name == 'proportional_metastatic') {
                clone1$invasion = ifelse( clone1$im > 0, TRUE, FALSE )  # condition here is related to the model
            } else {
                clone1$invasion = ifelse( clone1$im == 1, TRUE, FALSE )
            }
        }
    }

    # The new number of cells in the clone:
    clone1$N_cells = ifelse( (clone1$N_cells - N_die) > 0, (clone1$N_cells - N_die) , 0)

    N_new = clone1$N_cells   # the initial number to split / before trial / - all cells "want" to split

    # Fragmentation restriction trial
    if ( clone1$c > pck.env$ctmax  &  pck.env$tumbler_for_immortalization_trial ) {
        N_new = calc_binom( 1, N_new, (1 - clone1$i) )
    }

    # Divide trial
    N_new = calc_binom( 1, N_new, clone1$d )   # The ! FINAL ! number of cells to split, this number is the output value of function "trial"

    if ( clone1$N_cells > 0 ) clone1$c = clone1$c + N_new / clone1$N_cells  # How to calculate the counter of division ? - As an average value of the counters for all cells !

    clone1$N_cells = clone1$N_cells + N_new       # The number of cells are increased due to the splitting

    N_new_clones = 0             # The number of new clones arising from clone1

    # p = clone1$m * sum(onco$cds)
    # if (p < 1) N_new_clones = calc_binom(1, 2 * N_new, p ) else N_new_clones = 2 * N_new

    N_new_clones = calc_binom( 1, 2 * N_new, 1 - (onco1$p0_1 + onco1$p0_2 ) / 2  )   # 1 and 2 chr

    clone1$N_cells = ifelse( ( clone1$N_cells - N_new_clones ) > 0, ( clone1$N_cells - N_new_clones ) , 0 )

    return( N_new_clones )
}

#' Function trial for simplified case of model
#'
#' @param clone1 Object of class 'Clone'
#' @param onco1 Object of class 'OncoGene'
#'
#' @return Number of new clones originated by clone1
#' @export
#'
#' @examples
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' onco1 = tugHall_dataset$onco
#' trial_simple( clone1, onco1 )
#' unlist( lapply( X = 1:20, FUN = function( x ) trial_simple( clone1, onco1 ) ) )
trial_simple <- function( clone1, onco1 ) {

    # trial for Environmental death of cell
    N_die  =  calc_binom( 1, clone1$N_cells, clone1$k )   # The number of cells to die due to the Environmental death of cells in clone

    # DELETED: # Apoptosis trial
    # DELETED: N_die  =  N_die + calc_binom( 1, clone1$N_cells, clone1$a )

    # DELETED: # invasion / metastasis trial
    # if (clone1$im > 0) {
    #     if (!clone1$invasion) {
    #         N_die  =  N_die + calc_binom( 1, clone1$N_cells, ( 1 - clone1$im ) )
    #         if ( model_name == 'proportional_metastatic') {
    #             clone1$invasion = ifelse( clone1$im > 0, TRUE, FALSE )  # condition here is related to the model
    #         } else {
    #             clone1$invasion = ifelse( clone1$im == 1, TRUE, FALSE )
    #         }
    #     }
    # }

    # The new number of cells in the clone:
    clone1$N_cells = ifelse( (clone1$N_cells - N_die) > 0, (clone1$N_cells - N_die) , 0)

    # N_new = clone1$N_cells   # the initial number to split / before trial / - all cells "want" to split

    # DELETED: # Fragmentation restriction trial
    # if (clone1$c > 50) {
    #     N_new = calc_binom(1, N_new, (1 - clone1$i))
    # }

    # Divide trial
    N_new = calc_binom( 1, clone1$N_cells, clone1$d )   # The ! FINAL ! number of cells to split, this number is the output value of function "trial"

    if ( clone1$N_cells > 0 ) clone1$c = clone1$c + N_new / clone1$N_cells  # How to calculate the counter of division ? - As an average value of the counters for all cells !

    clone1$N_cells = clone1$N_cells + N_new       # The number of cells are increased due to the splitting

    N_new_clones = 0             # The number of new clones arising from clone1

    # p = clone1$m * sum(onco$cds)
    # if (p < 1) N_new_clones = calc_binom(1, 2 * N_new, p ) else N_new_clones = 2 * N_new

    N_new_clones = calc_binom( 1, 2 * N_new, 1 - (onco1$p0_1 + onco1$p0_2 ) / 2  )   # 1 and 2 chr

    clone1$N_cells = ifelse( ( clone1$N_cells - N_new_clones ) > 0, ( clone1$N_cells - N_new_clones ) , 0 )

    return( N_new_clones )
}

# mutagenesis trial

#' Function for mutagenesis trial
#'
#' @param clone1 Object of class 'Clone'
#' @param num_mut Number of mutations in this NEW clone1
#' @param onco1 Object of class 'OncoGene' corresponding to clone1 (with the same ID)
#'
#' @return Changed object clone1, add related mutations to the lists of point mutations and/or CNA mutations
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' copy_files_to_Output()
#' load_tugHall.Environment( tugHall_dataset )
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' onco1 = tugHall_dataset$onco_clones[[ 1 ]]
#' onco = tugHall_dataset$onco
#' define_gene_location()
#' pnt_clones = tugHall_dataset$pnt_clones
#' cna_clones = tugHall_dataset$cna_clones
#' mut_order = 234  # Just an example number
#' message( c('CNA mutation IDs ', paste(clone1$CNA_ID, collapse = ' ') ) )
#' message( c('Point mutation IDs ', paste(clone1$PointMut_ID, collapse = ' ') ) )
#' \dontrun{
#' trial_mutagenesis( clone1, num_mut = 1, onco1 )  # it adds info to clone1
#' message( c('CNA mutation IDs ', paste(clone1$CNA_ID, collapse = ' ') ) )
#' message( c('Point mutation IDs ', paste(clone1$PointMut_ID, collapse = ' ') ) )
#'
#' trial_mutagenesis( clone1, num_mut = 10, onco1 )  # it adds info to clone1
#' message( c('CNA mutation IDs ', paste(clone1$CNA_ID, collapse = ' ') ) )
#' message( c('Point mutation IDs ', paste(clone1$PointMut_ID, collapse = ' ') ) )
#' }
trial_mutagenesis <- function( clone1, num_mut, onco1 ) {

    # num_mut is a number of mutations in this NEW clone1
    # length of onco - the number of genes
    # onco1 is onco related to new clone1

    ### if (sum(mut1) == 0) stop("The mutation is zero, that is incorrect", call. = TRUE)    # We have known that mutation must occur
    type_events  =  sample( x = c('point_mut', 'del', 'dup' ), size = num_mut,
                            replace = TRUE, prob = (onco1$prob_1 + onco1$prob_2) / 2 )
    gm  =  modify_gene_map( clone1 , onco1 )
    # print( gm[[1]] )
    # print( gm[[2]] )
    for ( t in type_events ) {
        ### For each type of mutation to generate event of mutation
        if (t == 'point_mut') {
            pm  =  get_point_mutation( onco1, gm )
            prntl = unlist( pm[[1]] )
            gene  = unlist( pm[[2]] )
            pos   = unlist( pm[[3]] )
            Chr   = unlist( pm[[4]] )
            pnt0 = generate_pnt( prntl, gene, pos, onco1, Chr )
            if ( (clone1$PointMut_ID == 0)[1] ) {
                id   =  pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID
            } else  id   =  c( clone1$PointMut_ID, pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID )
            clone1$PointMut_ID  =  id   # pnt_clone is generated in the generate_pnt function

            if ( pnt0$MalfunctionedByPointMut ){
                clone1$gene[ which( pck.env$onco$name == gene ) ] = 1
            } else {
                clone1$pasgene[ which( pck.env$onco$name == gene ) ] = 1
            }
            gm[[ prntl ]]  =  add_pnt_mutation( gm = gm[[ prntl ]], pos_pnt = pos, Chr = Chr )
        }

        if (t == 'dup' | t == 'del') {
            cna_mut = get_cna_mutation( onco1, dupOrdel = t , gm_1_2 = gm)
            prntl =  unlist( cna_mut[[1]] )
            Chr   =  unlist( cna_mut[[2]] )
            genes =  unlist( cna_mut[[3]] )
            start_end   = unlist( cna_mut[[4]] )
            cna0 = generate_cna( prntl, genes, start_end, onco1, t )
            if ( (clone1$CNA_ID == 0)[1] ) {
                id   =  pck.env$cna_clones[[ length( pck.env$cna_clones ) ]]$CNA_ID
            } else  id   =  c( clone1$CNA_ID, pck.env$cna_clones[[ length( pck.env$cna_clones ) ]]$CNA_ID )
            clone1$CNA_ID  =  id

            if ( cna0$MalfunctionedByCNA ){
                clone1$gene[ sapply(genes, FUN = function(x) which(pck.env$onco$name == x) ) ] = 1
            } else {
                clone1$pasgene[ sapply(genes, FUN = function(x) which(pck.env$onco$name == x) ) ] = 1
            }

            ### Change the gene_map related chromosome
            ifelse( t == 'dup',
                    gm[[ prntl ]]  <-  add_duplication( gm = gm[[ prntl ]], Ref_start = start_end[1], Ref_end = start_end[2], Chr = Chr ),
                    gm[[ prntl ]]  <-     add_deletion( gm = gm[[ prntl ]], Ref_start = start_end[1], Ref_end = start_end[2], Chr = Chr )  )

            ### Check what point mutations match into the CNA
            sp   = FALSE
            sp_A = FALSE
            if ( clone1$PointMut_ID[ 1 ] != 0 ){
                sp = sapply( clone1$PointMut_ID , FUN = function( x )  {
                    chk_pnt_mut( pnt1  =  pck.env$pnt_clones[[ x ]], Ref_start = start_end[1],
                                 Ref_end = start_end[2], Chr = Chr, prntl  =  prntl )
                })
                ### to check original allele A do/don't match into CNA:
                prntl_inv  =  ifelse( prntl  ==  1, 2, 1 )
                sp_A = sapply( clone1$PointMut_ID , FUN = function( x )  {
                    chk_pnt_mut( pnt1  =  pck.env$pnt_clones[[ x ]], Ref_start = start_end[1],
                                 Ref_end = start_end[2], Chr = Chr, prntl  =  prntl_inv )
                })


            }

            if ( any( sp ) | any( sp_A ) ){
                ### Before changing point mutation for NEW clone
                ### we have to copy it and avoid any changes in OTHER (parental) clones

                ### Copy pnt mutations matched into CNA

                ### circle and function to copy pnt for a NEW clone
                for( q in clone1$PointMut_ID[ sp | sp_A ] ){
                    # pnt_clone is generated in the generate_to_copy_pnt function
                    pnt0 = generate_to_copy_pnt( pnt = pck.env$pnt_clones[[ q ]] )
                    x  =  which( clone1$PointMut_ID == q )
                    clone1$PointMut_ID[ x ]  =  pnt0$PointMut_ID
                }
                ### end of the new part of the code

                sapply( clone1$PointMut_ID[ sp ],
                        FUN = function( x ) change_pnt_by_cna( pnt1  =  pck.env$pnt_clones[[ x ]], start_end, t )  )
                sapply( clone1$PointMut_ID[ sp_A ],
                        FUN = function( x ) change_allele_A_by_cna( pnt1  =  pck.env$pnt_clones[[ x+1 ]], start_end, t )  )

            }

        }
    }

    onco_update( onco1, gm )
    if ( sum( clone1$gene ) > 0 ) clone1$primary = TRUE

}



# Functions related to point mutations ------------------------------------


#' Function to generate point mutations for initial clones
#'
#' @param clones List of objects of class 'Clone'
#' @param onco_clones List of objects of class 'OncoGene'
#'
#' @return NULL
#' @export
#'
#' @examples
#' clones = tugHall_dataset$clones
#' load_tugHall.Environment( tugHall_dataset )
#' onco = tugHall_dataset$onco
#' onco_clones = tugHall_dataset$onco_clones
#' copy_files_to_Input()
#' copy_files_to_Output()
#' define_gene_location()
#' pnt_clones = tugHall_dataset$pnt_clones
#' cna_clones = tugHall_dataset$cna_clones
#' mut_order = 234
#' \dontrun{
#' init_pnt_clones( clones, onco_clones )  # change pnt_clones for initialization
#' }
init_pnt_clones   <- function( clones, onco_clones ) {

    if ( is.null( clones ) )  return( NULL )

    for( i in 1:length( clones ) ){
        clone1  =  clones[[ i ]]
        onco1   =  onco_clones[[ i ]]

        genes  =  onco1$name[ clone1$gene == 1 ]
        gm  =  modify_gene_map( clone1 , onco1 )

        if ( length( genes ) == 0 ) {
            print('The clone of the normal cells is in simulation')  # stop('Length of mutated genes should be non-zero ') )
        } else {
            for( gene in genes ){

                pm  =   get_point_mutation_for_gene( onco1, gm_1_2 = gm, gene )   #  get_point_mutation( onco1, gm )
                prntl = unlist( pm[[1]] )
                # gene  = unlist( pm[[2]] )
                pos   = unlist( pm[[3]] )
                Chr   = unlist( pm[[4]] )
                pnt0 = generate_pnt( prntl, gene, pos, onco1, Chr, mutation = TRUE )

                ### Add pnt mutation ID to a clone:
                if ( (clone1$PointMut_ID == 0)[1] ) {
                    id   =  pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID
                } else  id   =  c( clone1$PointMut_ID, pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID )
                clone1$PointMut_ID  =  id
            }
        }
    }

}

#' Function to change the point mutation due to CNA
#'
#' @param pnt1 Object of class 'Point_Mutations'
#' @param start_end Vector with initial and final positions of CNA
#' @param t 'dup' or 'del' for duplication or deletion respectively
#'
#' @return NULL, but pnt1 data is updated due to CNA
#' @export
#'
#' @examples
#' pnt1 = tugHall_dataset$pnt_clones[[ 1 ]]
#' start_end = c( pnt1$Phys_pos - 50 , pnt1$Phys_pos + 50  )
#' message( pnt1$Copy_number )
#' change_pnt_by_cna( pnt1, start_end, t = 'dup' ) # View( pnt1 )
#' message( pnt1$Copy_number )
#' change_pnt_by_cna( pnt1, start_end, t = 'del' ) # View( pnt1 )
#' message( pnt1$Copy_number )
change_pnt_by_cna  <-  function( pnt1, start_end, t ) {

    ntrs  =  intersect( which( pnt1$Phys_pos >= start_end[1] ),
                        which( pnt1$Phys_pos <= start_end[2] )  )
    cf    =  length( ntrs )  ### coefficient for CN

    if ( pnt1$Copy_number != 0 )  {
        pn_cn  =  pnt1$Copy_number + cf * ifelse(t == 'dup', 1, -1 )
        pnt1$Copy_number  =  ifelse( pn_cn >= 0, pn_cn, 0 )
    }

    pos_pnt       =   pnt1$Phys_pos[ ntrs ]

    dlt  =  ifelse( t == 'dup', start_end[2]  -  start_end[1], start_end[1]  -  start_end[2] )
    pnt1$Phys_pos   =  c( pnt1$Phys_pos, pos_pnt + dlt )
    pnt1$Delta      =    pnt1$Phys_pos  -  pnt1$Ref_pos

}

#' Function to change copy number of the allele A of the point mutation at the allele B due to CNA
#'
#' @param pnt1 Object of class 'Point_Mutations'
#' @param start_end Vector with initial and final positions of CNA
#' @param t 'dup' or 'del' for duplication or deletion respectively
#'
#' @return NULL, but data of pnt1 is updated due to CNA
#' @export
#'
#' @examples
#' pnt1 = tugHall_dataset$pnt_clones[[ 2 ]]  # pnt of allele A
#' start_end = c( pnt1$Phys_pos - 50 , pnt1$Phys_pos + 50  )
#' message( pnt1$Copy_number )
#' change_allele_A_by_cna( pnt1, start_end, t = 'dup' )  # View( pnt1 )
#' message( pnt1$Copy_number )
#' change_allele_A_by_cna( pnt1, start_end, t = 'del' )  # View( pnt1 )
#' message( pnt1$Copy_number )
change_allele_A_by_cna  <-  function( pnt1, start_end, t ) {
    if ( pnt1$Copy_number != 0 )  {
        pn_cn  =  pnt1$Copy_number + ifelse(t == 'dup', 1, -1 )
        pnt1$Copy_number  =  ifelse( pn_cn >= 0, pn_cn, 0 )
    }
}

#' Function to update onco1 after mutation (for usage in trial_mutagenesis() function)
#'
#' @param onco1 Object of class 'OncoGene'
#' @param gm data.frame gene_map
#'
#' @return onco1 with updated info
#' @export
#'
#' @examples
#' onco1 = tugHall_dataset$onco_clones[[ 1 ]]
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' withr::with_environment( env = pck.env, code = onco_update( onco1, gm = list(gene_map, gene_map[1:42, ] ) )  )
#' # Check CDS length for TP53 gene
onco_update  <-  function( onco1, gm ){

    lst1  =  get_cds_rna( gm[[1]] )
    rd1   =  as.integer( sapply( onco1$name, FUN = function(x) which( x  ==  lst1[[1]]) ) )

    lst2  =  get_cds_rna( gm[[2]] )
    rd2   =  as.integer( sapply( onco1$name, FUN = function(x) which( x  ==  lst2[[1]]) ) )

    # change the onco1 related to new gene_map:
    onco1$cds_1  =  lst1$CDS[ rd1 ]
    onco1$cds_2  =  lst2$CDS[ rd2 ]

    onco1$rna_1  =  lst1$RNA[ rd1 ]
    onco1$rna_2  =  lst2$RNA[ rd2 ]

    onco1$prob_1 =  lst1$PROB
    onco1$prob_2 =  lst2$PROB

    onco1$sum_prob_1  =  lst1$SUM
    onco1$sum_prob_2  =  lst2$SUM

    onco1$p0_1   =  lst1$P0
    onco1$p0_2   =  lst2$P0

    return( onco1 )
}

#' Function to check point mutations match or don't match into duplication or deletion
#'
#' @param pnt1 Object of class 'Point_Mutations'
#' @param Ref_start Initial position of CNA
#' @param Ref_end Final position of CNA
#' @param Chr Chromosome name
#' @param prntl Parental chromosome 1 or 2
#'
#' @return Logical: TRUE if point mutation matches CNA, FALSE if it doesn't match
#' @export
#'
#' @examples
#' pnt1 = tugHall_dataset$pnt_clones[[ 5 ]]
#' pstn = pnt1$Phys_pos[1]
#' message( pstn )
#' prntl = pnt1$Parental_1or2
#' Chr = pnt1$Chr
#' chk_pnt_mut( pnt1 , Ref_start = pstn - 200, Ref_end = pstn + 200, Chr, prntl )
#' chk_pnt_mut( pnt1 , Ref_start = pstn - 200, Ref_end = pstn - 100, Chr, prntl )
chk_pnt_mut  <-  function( pnt1 , Ref_start, Ref_end, Chr, prntl ){

    for( X in pnt1$Phys_pos ){
        if ( pnt1$Chr == Chr &  X <= Ref_end & X >= Ref_start & pnt1$Parental_1or2 == prntl ) {
            return( TRUE )}
    }
    return( FALSE )
}


# Functions related to simulation -----------------------------------------

#' Function to make one copy for clone1 in clone_init function
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return New object of class 'Clone' with the same info and new ID
#' @export
#'
#' @examples
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' load_tugHall.Environment( tugHall_dataset )
#' withr::with_environment( pck.env, code = clone_copy(clone1) )
clone_copy <- function(clone1) {
    pck.env$env$last_id = pck.env$env$last_id + 1

    return(clone$new(id=pck.env$env$last_id, parent=clone1$id, c=clone1$c, d=clone1$d,
                     i=clone1$i, mutden=clone1$mutden, a=clone1$a,
                     k=clone1$k, E=clone1$E, Nmax=clone1$Nmax,
                     gene=clone1$gene, pasgene=clone1$pasgene,
                     PointMut_ID = clone1$PointMut_ID, CNA_ID = clone1$CNA_ID,
                     # posdriver=clone1$posdriver, pospasngr=clone1$pospasngr,
                     invasion=clone1$invasion, primary=clone1$primary,
                     s=clone1$s, birthday=pck.env$env$T)) #,
    # len = clone1$len ))
}

#' Function to make one copy for onco1 in init_onco_clones function
#'
#' @param onco1 Object of class 'OncoGene'
#'
#' @return New object of class 'OncoGene' with the same info
#' @export
#'
#' @examples
#' onco1  =  tugHall_dataset$onco_clones[[ 1 ]]
#' onco2 = onco_copy( onco1 )  # ID + 1
onco_copy <- function( onco1 ){

    onco2 = oncogene$new(id = onco1$id, name = onco1$name, onsp = onco1$onsp, len = onco1$len,
                         cds_1 = onco1$cds_1, cds_2 = onco1$cds_2,
                         rna_1 = onco1$rna_1, rna_2 = onco1$rna_2,
                         p0_1 = onco1$p0_1, p0_2 = onco1$p0_2,
                         prob_1 = onco1$prob_1, prob_2 = onco1$prob_2,
                         sum_prob_1 = onco1$sum_prob_1, sum_prob_2 = onco1$sum_prob_2  )
    return( onco2 )
}

#' Function to write log file
#'
#' @param genefile File name of initial OncoGene information
#' @param clonefile File name of info about initial clones
#' @param geneoutfile File name for output info about OncoGene information
#' @param cloneoutfile File name for output info with clone evolution data
#' @param logoutfile Name of log file with all the parameters
#' @param   E0 Parameter in the division probability, numeric type only
#' @param   F0 Parameter in the division probability, numeric type only
#' @param   m0 Mutation probability for point mutation, numeric type only
#' @param   uo Oncogene mutation probability, numeric type only
#' @param   us Suppressor mutation probability, numeric type only
#' @param   s0 Parameter in the sigmoid function, numeric type only
#' @param   k0 Environmental death probability, numeric type only
#' @param   d0 Initial probability to divide cells, numeric type only
#' @param   ctmax Hayflick limitation for cell division, integer type
#' @param   censor_cells_number Max cell number where the program forcibly stops, integer type only
#' @param   censor_time_step Max time where the program forcibly stops, integer type only
#' @param   real_time_stop Max time in seconds of running after that the program forcibly stops, integer type only
#' @param   n_repeat  Max number of repetition of the program until the NON-ZERO output will be getting, integer type only
#' @param  m_dup Mutation probability for duplication, numeric type only
#' @param  m_del Mutation probability for deletion, numeric type only
#' @param  lambda_dup  CNA duplication average length (of the geometrical distribution for the length), integer type only
#' @param  lambda_del  CNA deletion average length (of the geometrical distribution for the length), integer type only
#' @param  uo_dup Gene malfunction probability by CNA duplication for oncogene, numeric type only
#' @param  us_dup Gene malfunction probability by CNA duplication for suppressor, numeric type only
#' @param  uo_del Gene malfunction probability by CNA deletion    for oncogene, numeric type only
#' @param  us_del Gene malfunction probability by CNA deletion    for suppressor, numeric type only
#' @param  monitor The indicator to make monitor file during a simulation or do not make, logical type only
#' @param Compaction_factor Compaction factor, logical type only. True means 'to use', False means 'do not use' Compaction factor for hallmarks variables
#' @param model_name Name of the model to use. Can be  'proportional_metastatic' or 'threshold_metastatic' or 'simplified'
#' @param tumbler_for_metastasis_trial Logical parameter to turn on/off invasion/metastasis transformation trial
#' @param tumbler_for_apoptosis_trial Logical parameter to turn on/off the apoptosis trial
#' @param tumbler_for_immortalization_trial Logical parameter to turn on/off the immortalization trial
#' @param tumbler_for_angiogenesis_trial Logical parameter to turn on/off angiogenesis trial
#' @param tumbler_for_drug_intervention_trial Logical parameter to turn on/off drug intervention trial
#'
#' @return NULL, write log file to Output folder
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' define_files_names()
#' load_tugHall.Environment( tugHall_dataset )
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' \dontrun{
#' write_log(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile,
#' E0, F0, m0, uo, us, s0, k0, ctmax, m_dup, m_del, lambda_dup, lambda_del,
#' uo_dup, us_dup, uo_del, us_del, censor_cells_number, censor_time_step, d0,
#' Compaction_factor, model_name, real_time_stop, n_repeat, monitor )
#' }
write_log <- function(genefile, clonefile, geneoutfile, cloneoutfile, logoutfile,
                      E0, F0, m0, uo, us, s0, k0, ctmax,
                      m_dup, m_del, lambda_dup, lambda_del, # CNA parameters
                      uo_dup, us_dup, uo_del, us_del,       # CNA parameters
                      censor_cells_number, censor_time_step, d0, Compaction_factor, model_name,
                      real_time_stop, n_repeat, monitor,
                      tumbler_for_metastasis_trial, tumbler_for_apoptosis_trial,
                      tumbler_for_immortalization_trial, tumbler_for_angiogenesis_trial,
                      tumbler_for_drug_intervention_trial ) {
    data <- c('Working_folder', "genefile", "clonefile", "geneoutfile", "cloneoutfile", "logoutfile",
              "E0", "F0", "m0", "uo", "us", "s0", "k0", 'ctmax',
              "m_dup", "m_del", "lambda_dup", "lambda_del",
              "uo_dup", "us_dup", "uo_del", "us_del",
              "censor_cells_number", "censor_time_step", "d0", 'Compaction_factor', 'model_name',
              'real_time_stop', 'n_repeat', 'monitor',
              'tumbler_for_metastasis_trial', 'tumbler_for_apoptosis_trial',
              'tumbler_for_immortalization_trial', 'tumbler_for_angiogenesis_trial',
              'tumbler_for_drug_intervention_trial' )
    data <- rbind( data, c(pck.env$mainDir,
                           str_remove( genefile, pck.env$mainDir),
                           str_remove( clonefile, pck.env$mainDir),
                           str_remove( geneoutfile, pck.env$mainDir),
                           str_remove( cloneoutfile, pck.env$mainDir),
                           str_remove( logoutfile, pck.env$mainDir),
                           E0, F0, m0, uo, us, s0, k0, ctmax,
                           m_dup, m_del, lambda_dup, lambda_del, # CNA parameters
                           uo_dup, us_dup, uo_del, us_del,       # CNA parameters
                           censor_cells_number, censor_time_step, d0, Compaction_factor, model_name,
                           real_time_stop, n_repeat, monitor,
                           pck.env$tumbler_for_metastasis_trial, pck.env$tumbler_for_apoptosis_trial,
                           pck.env$tumbler_for_immortalization_trial, pck.env$tumbler_for_angiogenesis_trial,
                           pck.env$tumbler_for_drug_intervention_trial ) )
    write(data, logoutfile, ncolumns = 2, sep="\t")
}


#' Function to write info about HallMark data
#'
#' @param outfile File name for output info
#' @param hall Object of class "HallMark"
#' @param Compaction_factor Compaction factor, logical type only. True means 'to use', False means 'do not use' Compaction factor for hallmarks variables
#' @param CF Vector with values of compaction factor for each hallmark
#'
#' @return NULL, but data will save to a file
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' withr::with_environment( env = pck.env, code =  write_geneout(outfile = geneoutfile, hall, Compaction_factor, CF) )
write_geneout <- function(outfile, hall, Compaction_factor, CF) {
    gene <- c(pck.env$onco$name[hall$Ha], pck.env$onco$name[hall$Hi], pck.env$onco$name[hall$Hd],
              pck.env$onco$name[hall$Hb], pck.env$onco$name[hall$Him])
    hall_mark <- c(rep("apoptosis", length(pck.env$onco$name[hall$Ha])),
                   rep("immortalization", length(pck.env$onco$name[hall$Hi])),
                   rep("growth|anti-growth", length(pck.env$onco$name[hall$Hd])),
                   rep("angiogenesis", length(pck.env$onco$name[hall$Hb])),
                   rep("invasion", length(pck.env$onco$name[hall$Him])))
    weight_CF <- c(hall$Ha_w, hall$Hi_w, hall$Hd_w, hall$Hb_w, hall$Him_w)
    if ( Compaction_factor ) {
        weight <- c(hall$Ha_w / CF$Ha, hall$Hi_w / CF$Hi, hall$Hd_w / CF$Hd, hall$Hb_w / CF$Hb, hall$Him_w / CF$Him)
    }
    else {
        weight <- weight_CF
    }
    sup_or_onc <- c(pck.env$onco$onsp[hall$Ha], pck.env$onco$onsp[hall$Hi], pck.env$onco$onsp[hall$Hd],
                    pck.env$onco$onsp[hall$Hb], pck.env$onco$onsp[hall$Him] )
    df <- data.frame(gene=gene, hall_mark=hall_mark, weight=weight, weight_CF=weight_CF, sup_or_onc=sup_or_onc)
    write.table(df, outfile, quote=F, row.names=F, sep='\t')
}

#' Function to write the header to a file
#'
#' @param outfile File name for output info
#' @param env Object of class 'Environ'
#' @param onco Object of class "OncoGene"
#'
#' @return NULL, but the header will save to a file and delete old info
#' @export
#'
#' @examples
#' env = tugHall_dataset$env
#' onco = tugHall_dataset$onco
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' write_header(outfile='./Output/exmpl.txt', env, onco)
write_header <- function(outfile, env, onco) {
    header <- c('Time', 'AvgOrIndx', 'ID', 'N_cells', 'ParentID', 'Birth_time', 'c', 'd', 'i', 'im', 'a',
                'k', 'E', 'N_normal_intact',	'N_normal_speckled', 'Nmax', 'N_primary', 'N_metastatic', 'Ha', 'Him', 'Hi', 'Hd', 'Hb', 'type', 'mut_den',
                # paste("PosDriver:", onco$name, sep=""), paste("PosPasngr:", onco$name, sep="") )      #   , 'Clone number', 'Passengers Clone number', 'Mix Clone number')
                'driver_genes', 'passenger_genes',
                'PointMut_ID', 'CNA_ID' ,
                paste('Chr1_CDS_', onco$name, sep=''), paste('Chr1_Len_', onco$name, sep=''),
                'Chr1_p0', 'Chr1_prob_point_mut', 'Chr1_prob_del', 'Chr1_prob_dup' ,
                paste('Chr2_CDS_', onco$name, sep=''), paste('Chr2_Len_', onco$name, sep=''),
                'Chr2_p0', 'Chr2_prob_point_mut', 'Chr2_prob_del', 'Chr2_prob_dup' )
    write(header, outfile, append=FALSE, ncolumns = length(header), sep="\t")
}

#' Function to get type of the clone: normal, primary or metastatic
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return One of characters 'normal', 'primary' or 'metastatic'
#' @export
#'
#' @examples
#' clone1 = tugHall_dataset$clones[[1]]
#' get_type( clone1 )
#' clone1 = tugHall_dataset$clones[[56]]
#' get_type( clone1 )
get_type  <-  function( clone1 ){
    if ( clone1$invasion ) return( 'metastatic' )
    if ( !clone1$primary ) return( 'normal' )
    return( 'primary' )
}

#' Function to write data to cloneout file at a time step
#'
#' @param outfile File name for output info
#' @param env Object of class 'Environ'
#' @param clones List of objects of class 'Clone'
#' @param isFirst logical type = TRUE as default
#' @param onco_clones List of objects of class 'OncoGene'
#'
#' @return NULL, but add rows to output file with clone evolution data
#' @export
#'
#' @examples
#' env = tugHall_dataset$env
#' onco = tugHall_dataset$onco
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' clones = tugHall_dataset$clones
#' onco_clones = tugHall_dataset$onco_clones
#' write_header(outfile='./Output/exmpl.txt', env, onco)
#' write_cloneout( outfile = './Output/exmpl.txt', env, clones, isFirst = TRUE, onco_clones )
write_cloneout <- function( outfile, env, clones, isFirst, onco_clones ) {
    intact_normal  =  sum( sapply( clones, FUN = function( cl ) ifelse( cl$CNA_ID[ 1 ] == 0 & cl$PointMut_ID[ 1 ] == 0, cl$N_cells, 0 ) ) )
    data  =  c(env$T, 'avg', '-',  '-', '-', '-', env$c, env$d, env$i, env$im, env$a, env$k, env$E,
               intact_normal, env$N - intact_normal,  # env$N,
               env$Nmax, env$P, env$M, env$Ha, env$Him, env$Hi, env$Hd, env$Hb, '-', env$mutden,
               rep('-', 12 + 4*length( pck.env$onco$name ) )      # rep('-',17)
    )

    write(data, outfile, append=TRUE, ncolumns = length(data), sep="\t")

    if (length(clones) > 0 & isFirst) {
        for (i in 1:length( clones ) ) {
            clone1  =  clones[[i]]
            onco1   =  onco_clones[[i]]
            type    =  get_type( clone1 = clone1 )
            data    =  c(env$T, i, clone1$id, clone1$N_cells, clone1$parent, clone1$birthday, clone1$c, clone1$d,
                         clone1$i, clone1$im, clone1$a, clone1$k, clone1$E,
                         intact_normal, env$N - intact_normal, # env$N,
                         clone1$Nmax, env$P, env$M,
                         clone1$Ha, clone1$Him, clone1$Hi, clone1$Hd, clone1$Hb, type,  # ifelse(clone1$invasion,1,0),
                         clone1$mutden,
                         paste(clone1$gene, collapse =  ' '), paste(clone1$pasgene, collapse =  ' '),
                         paste(clone1$PointMut_ID, collapse = ', '), paste(clone1$CNA_ID, collapse = ', '),
                         # onco1$id,
                         onco1$cds_1, onco1$rna_1, onco1$p0_1, onco1$prob_1,
                         onco1$cds_2, onco1$rna_2, onco1$p0_2, onco1$prob_2      )
            write(data, outfile, append=TRUE, ncolumns = length(data), sep="\t")
        }
    }
}

#' Function to write a simulation monitoring data into the file_monitor
#'
#' @param outfile File name for output info
#' @param start Indicator to start from beginning (TRUE) or not (FALSE)
#' @param env Object of class 'Environ'
#' @param clones List of objects of class 'Clone'
#'
#' @return NULL, but info about current state of simulation will write to a file
#' @export
#'
#' @examples
#' env = tugHall_dataset$env
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' clones = tugHall_dataset$clones
#' onco_clones = tugHall_dataset$onco_clones
#' cna_clones = tugHall_dataset$cna_clones
#' pnt_clones = tugHall_dataset$pnt_clones
#' write_monitor( outfile = './Sim_monitoring.txt', start = TRUE , env, clones )
#' write_monitor( outfile = './Sim_monitoring.txt', start = FALSE , env, clones )
write_monitor  <- function( outfile, start = FALSE, env, clones ){

    if ( start ) {
        header <- c('Time', 'N_clones', 'N_normal_intact',  'N_normal_speckled', 'N_primary', 'N_metastatic',
                    'N_point_mutations', 'N_duplications',   'N_deletions' , 'TMB', 'TMB%5', 'TMB%10' )
        write( header, outfile, append = FALSE, ncolumns = length( header ), sep="\t" )
    } else {
        if ( length( clones )  >  0 ) {

            point_mut_list  =  sort( unique( as.numeric( unlist( sapply( X = 1:length( clones ), FUN = function( x ) clones[[ x ]]$PointMut_ID ) ) ) ) )
            if ( point_mut_list[ 1 ] == 0 ) l_pm  =  length( point_mut_list ) - 1 else l_pm  =  length( point_mut_list )
            cna_list  =  sort( unique( as.numeric( unlist( sapply( X = 1:length( clones ), FUN = function( x ) clones[[ x ]]$CNA_ID ) ) ) ) )
            if ( cna_list[ 1 ] == 0 ) cna_list  =  cna_list[ -1 ]
            dupdel  =  unlist( sapply( X = cna_list, FUN = function( x ) pck.env$cna_clones[[ x ]]$dupOrdel ) )
            l_dup   =  length( which( dupdel  ==  'dup' ) )
            l_del   =  length( which( dupdel  ==  'del' ) )
            TMB     =  c(0,0,0) # get_TMB( env = env, clones = clones, pnt_clones = pck.env$pnt_clones )   #      l_pm * 1E06 / sum( pck.env$onco$cds_1 ) / ( env$N + env$P + env$M )

            # Get intact and speckled normal cells:
            i_n  =  which( sapply( 1:length(clones), FUN = function(x) get_type( clones[[ x ]] ) ) == 'normal')
            if ( length( i_n ) > 0 ){

                int = sapply( i_n, FUN = function( x ) sum( clones[[ x ]]$pasgene ) )
                if ( length( which( int == 0 ) ) > 0 ){
                    N_intact  =  sum( sapply( which( int == 0 ), FUN = function( x ) clones[[ x ]]$N_cells ) )
                } else N_intact    =  0

                N_speckled  =  env$N - N_intact
            } else {
                N_intact    =  0
                N_speckled  =  0
            }

            data <- c( env$T, length( clones ), N_intact, N_speckled, env$P, env$M, l_pm, l_dup, l_del, TMB )

            write(data, outfile, append=TRUE, ncolumns = length(data), sep="\t")
        }


    }

}


#' @describeIn write_monitor  Function to get TMB - number of point mutations per 10^6 bps (per M bps)
#'
#' @param pnt_clones list of point mutations
#'
#' @return Vector of three numbers: TMB, TMB for point mutation with VAF > 5% and TMB for point mutation with VAF > 10%
#'
#' @export
#'
#' @examples
#' NULL
#'
get_VAF_primary_clones  <- function( env, clones, pnt_clones ){

    if ( is.null(clones) ) return( NULL )

    l   =   ( 1:(length(pnt_clones)/2) ) * 2 - 1
    pnt_mut_B  =  lapply( l, FUN = function( x ) pnt_clones[[ x ]] )
    pnt_mut_A  =  lapply( l+1, FUN = function( x ) pnt_clones[[ x ]] )

    # l do not choose invasion clones:
    l    =  which( !sapply( 1:length( clones ), FUN = function( x ) clones[[ x ]]$invasion ) )
    cl   =  lapply( l, FUN = function( x ) clones[[ x ]] )

    # indexes of point mutations for each clone:
    ids  =  lapply( 1:length( cl ), FUN = function( x ) cl[[ x ]]$PointMut_ID )
    nqu  =  sort( unique( unlist( ids ) ) )  # unique IDs of point mutations
    if ( nqu[1] == 0 ) nqu = nqu[ -1 ]       # exclude intact normal cells

    # start VAF
    VAF  =  NULL
    if ( length( nqu ) == 0 ) return( NULL )

    # get types for all the clones:
    clone_types  =  sapply( 1:length( cl ), FUN = function( x ) get_type( cl[[ x ]]) )
    # get number of cells for ll the clones:
    clone_numbers  =  sapply( 1:length( cl ), FUN = function( x ) cl[[ x ]]$N_cells )
    # iteration across all the unique point mutations:
    for( j in 1:length( nqu ) ){
        # wc - which clones have an ID of point mutation
        wc  =  which( sapply( X = 1:length( ids ), FUN = function( x ) is.element( nqu[j] , ids[[ x ]] ) ) )

        VAF_1  =  pnt_mut_B[[ (nqu[ j ] + 1) / 2 ]]$safe( )

        ### number of cells: speckled normal cells and primary tumor cells:
        if ( any( clone_types[ wc ] == 'normal' ) ){
            sm  =  sum( clone_numbers[ wc[ which( clone_types[ wc ]  ==  'normal' ) ] ] )
            if ( length( sm ) == 0 ) sm = 0  # if all the normal clones are intact normal
            VAF_1$N_speckled_normal =  sm
        } else VAF_1$N_speckled_normal =  0

        if ( any( clone_types[ wc ] == 'primary' ) ){
            VAF_1$N_primary  =  sum( clone_numbers[ wc[ which( clone_types[ wc ]  ==  'primary' ) ] ] )
        } else VAF_1$N_primary  =  0

        # add copy number of original allele A:
        # VAF_A  =  pnt_mut_A[ which( pnt_mut_A$PointMut_ID == nqu[ j ] ) , ]
        VAF_1$Copy_number_A  =  pnt_mut_A[[ ( nqu [ j ] + 1 ) / 2 ]]$Copy_number

        VAF  =  rbind( VAF, VAF_1)
    }

    i_n  =  which( sapply( 1:length(clones), FUN = function(x) get_type( clones[[ x ]] ) ) == 'normal')
    if ( length( i_n ) > 0 ){
        int = sapply( i_n, FUN = function( x ) sum( clones[[ x ]]$pasgene ) )
        if ( length( which( int == 0 ) ) > 0 ){
            N_intact  =  sum( sapply( which( int == 0 ), FUN = function( x ) clones[[ x ]]$N_cells ) )
        } else N_intact    =  0

        N_speckled  =  env$N - N_intact
    } else {
        N_intact    =  0
        N_speckled  =  0
    }

    VAF$N_speckled_normal_total  =  N_speckled
    VAF$N_primary_total          =  env$P
    VAF$N_metastatic   =    0
    VAF$N_metastatic_total   =   env$M

    return( VAF )
}


#' Function to write info about relationship between genes and hallmarks
#'
#' @param outfile File name for output info
#' @param hall Object of class 'HallMark'
#'
#' @return NULL, but info about relationship between genes and hallmarks will write to a file
#' @export
#'
#' @examples
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' hall = tugHall_dataset$hall
#' onco = tugHall_dataset$onco
#' write_weights(outfile = './Output/weights.txt', hall)
write_weights <- function(outfile, hall) {
    #data <- c("Hallmarks", "Designation", onco$name)
    data <- data.frame( "Gene" = pck.env$onco$name)
    data$Gene <-   as.character(data$Gene)

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Ha) ) > 0  ) w[j] = hall$Ha_w[which(j==hall$Ha)]  }
    data <- cbind(data, "Apoptosis ($H_a$)" = w)

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Hb) ) > 0  ) w[j] = hall$Hb_w[which(j==hall$Hb)]  }
    data <- cbind(data, "Angiogenesis ($H_b$)" = w)

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Hd) ) > 0  ) w[j] = hall$Hd_w[which(j==hall$Hd)]  }
    data <- cbind(data, "Growth / Anti-growth ($H_d$)" = w)

    w <- rep( 0.0, pck.env$onco$len )
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Hi) ) > 0  ) w[j] = hall$Hi_w[which(j==hall$Hi)]  }
    data <- cbind(data, "Immortalization ($H_i$)" = w)

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Him) ) > 0  ) w[j] = hall$Him_w[which(j==hall$Him)]  }
    data <- cbind(data, "Invasion / Metastasis ($H_{im}$)" = w)

    write.table(data, outfile, sep="\t", row.names = FALSE)
}


#' @describeIn write_weights Function to write info about relationship between genes and hallmarks in the framework of break points
#'
#'
#' @return \code{write_break_points} returns NULL,
#' but break points of weights between genes and hallmarks will write to a file
#'
#' @export
#'
#' @examples
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' hall = tugHall_dataset$hall
#' onco = tugHall_dataset$onco
#' write_break_points(outfile = './Output/break_points.txt', hall)
write_break_points <- function(outfile, hall) {
    #data <- c("Hallmarks", "Designation", onco$name)
    data <- data.frame( "Gene" = pck.env$onco$name)
    data$Gene <-   as.character( data$Gene )

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) {
        if ( length( which(j==hall$Ha) ) > 0  ) w[j] = hall$Ha_w[which(j==hall$Ha)]
    }
    data <- cbind(data, "Apoptosis" = w)

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Hb) ) > 0  ) w[j] = hall$Hb_w[which(j==hall$Hb)]  }
    data <- cbind(data, "Angiogenesis" = w)

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Hd) ) > 0  ) w[j] = hall$Hd_w[which(j==hall$Hd)]  }
    data <- cbind(data, "Growth_Antigrowth" = w)

    w <- rep( 0.0, pck.env$onco$len )
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Hi) ) > 0  ) w[j] = hall$Hi_w[which(j==hall$Hi)]  }
    data <- cbind(data, "Immortalization" = w)

    w <- rep(0.0, pck.env$onco$len)
    for (j in 1:pck.env$onco$len) { if ( length( which(j==hall$Him) ) > 0  ) w[j] = hall$Him_w[which(j==hall$Him)]  }
    data <- cbind(data, "Invasion_Metastasis" = w)

    # Correct normalization
    if ( pck.env$Compaction_factor ){

        data$Apoptosis     =  data$Apoptosis / pck.env$CF$Ha
        data$Angiogenesis  =  data$Angiogenesis / pck.env$CF$Hb
        data$Growth_Antigrowth  =  data$Growth_Antigrowth / pck.env$CF$Hd
        data$Immortalization    =  data$Immortalization   / pck.env$CF$Hi
        data$Invasion_Metastasis  =  data$Invasion_Metastasis / pck.env$CF$Him
    }

    data$break_points = paste0( 'break_point_', data$Gene )
    data = data[, -1]
    data = cbind( data[ , 6], data[, -6 ])
    names( data ) = c( 'break_points', names( data[ , -1]))
    for( i in 2:nrow( data ) ){
        data[ i, 2:6]  =  data[ i, 2:6] + data[ i-1, 2:6]
    }

    write.table(data, outfile, sep="\t", row.names = FALSE)
}

#' Function to write the point mutation info for all clones for all time steps, used at the last time step or after simulation
#'
#' @param pnt_clones List of objects of class 'Point_Mutations'
#' @param file_out File name to write
#'
#' @return NULL, but info will write to a file
#' @export
#'
#' @examples
#' pnt_clones = tugHall_dataset$pnt_clones
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' write_pnt_clones(pnt_clones, file_out = 'Output/point_mutations.txt')
write_pnt_clones <- function( pnt_clones, file_out = 'Output/point_mutations.txt' ){
    # if ( is.null(pnt_clones) ) return( print( 'Empty data.' ) )
    pn  <-  NULL
    if ( !is.null( pnt_clones ) ){
        for (i in 1:length(pnt_clones)) {
            pnt1 <-  unlist( pnt_clones[[i]] )
            pn1  <-  safe_pnt_mut( pnt1 )     ### pnt1$safe()
            pn   <-  rbind( pn, pn1)
        }
    }

    write.table( pn, file = file_out, append=FALSE, sep="\t", row.names = FALSE)
}

#' Function to read file with initial clones
#'
#' @param clonefile File to read
#' @param clone1 Object of class 'Clone'
#'
#' @return List of objects of class 'Clone
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' withr::with_environment( env = pck.env, code = init_clones(clonefile, clone1) )
init_clones <- function(clonefile, clone1) {
    mpos <- regexpr("\\.", clonefile)[1]
    if (mpos != -1) {
        name <- substr(clonefile, 1, mpos - 1)
    } else {
        name <- clonefile
    }
    clones = NULL
    n <- as.numeric(name)
    if (!is.na(n) && is.numeric(n)) {
        factor = n / sum(clone1$m * pck.env$onco$cds_1)
        f2 = 1.0
        while (TRUE) {
            if (sum(floor(clone1$m * pck.env$onco$cds_1*factor*f2 + 0.5)) >= n) {
                break
            }
            f2 = f2 + 0.1
        }
        nums = floor(clone1$m * pck.env$onco$cds_1*factor*f2 + 0.5)
        clones = NULL
        for (i in 1:n) {
            clones = c(clones, clone_copy(clone1))
        }
        pos = 0
        for (i in 1:length(nums)) {
            if (nums[i] > 0) {
                for (j in 1:nums[i]) {
                    if (pos + j <= n) {
                        clones[[pos + j]]$gene[i] = 1
                    }
                }
                pos = pos + nums[i]
            }
        }
    } else {
        data = read.table(clonefile, sep="\t", stringsAsFactors=FALSE )
        n <- nrow(data)

        for (i in 1:n) {
            clone2 = clone_copy(clone1)
            p <- match( pck.env$onco$name, str_trim(strsplit(as.character(data[i,2]),",")[[1]]))
            clone2$gene[seq(1,length( pck.env$onco$name))[!is.na(p)]] = 1
            clone2$N_cells = as.numeric(data[i,3])
            clones = c(clones, clone2)
        }
    }
    for (i in 1:n) {
        clones[[i]]$id = i
        clones[[i]]$parent = 0
        clones[[i]]$birthday = 0
        # clones[[i]]$posdriver = ifelse(clones[[i]]$gene == 1,
        #                               paste(ceiling(runif(onco$len)*onco$cds),"0",sep = ":"),
        #                               clones[[i]]$posdriver)
        clones[[i]]$calcMutden()
        clones[[i]]$calcApoptosis()
        if ( sum( clones[[ i ]]$gene ) > 0 ) clones[[ i ]]$primary = TRUE
    }
    pck.env$env$last_id = n
    return( as.list( clones ) )
}

#' Function to make list of objects of class 'OncoGene' and generate initial onco settings for all clones (onco_clones)
#'
#' @param onco1 Object of class 'OncoGene'
#' @param clones List of objects of class 'Clone'
#'
#' @return List of objects of class 'OncoGene'
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' clone1 = tugHall_dataset$clones[[ 1 ]]
#' withr::with_environment( env = pck.env, code = { clones = init_clones(clonefile, clone1)  } )
#' withr::with_environment( env = pck.env, code = { onco_clones = init_onco_clones( onco1 = onco, clones ) } )
init_onco_clones <- function( onco1, clones ) {
    # initialization of onco_clones in order to assign an onco to each clone
    onco_clones = NULL
    for ( i in 1:length( clones ) ) {
        clone1 = clones[[i]]
        onco_clone2 = onco_copy( onco1 )
        onco_clone2$id = clone1$id
        onco_clones = c( onco_clones, onco_clone2 )
    }

    return( as.list( onco_clones ) )
}

#' Function to calculate binomial distribution including BIG NUMBERS like 10^12 and more using approximation with normal distribution
#'
#' @param tr Length of vector with successes trials
#' @param n Number of independent Bernoulli trials
#' @param p Probability to get successes in trials
#'
#' @return Vector of integer numbers of successes trials
#' @export
#'
#' @examples
#' calc_binom(tr = 3, n = 40, p = 0.9)
#' calc_binom(tr = 3, n = 4E20, p = 9E-9)
calc_binom <- function(tr,n,p){
    if (n*p < 10^8){
        ou <- rbinom(tr,n,p)
    } else {
        m <- n * p
        s <- sqrt(  n*p* (1-p)  )
        ou <- rnorm( tr, mean = m, sd = s )
    }

    return(  round( ou )  )
}
