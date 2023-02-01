


#'  Function to emulate drug intervention to pool of tumor cells by killing a part of these cells
#'  and blocking a part of clones corresponding to malfunctioned gene
#'
#' @param kill_prob Probability of killing cancer cells corresponding to the malfunctioned gene
#' @param block_prob Probability of blocking cancer cells corresponding to the malfunctioned gene
#' @param gene Name of target gene to kill and block tumor cells by a drug
#' @param generate_mutations Logical to generate or not new mutations states with
#' the same positions but for passenger genes instead drivers
#'
#' @return  NULL changing clones and onco_clones objects in tugHall environment pck.env
#' @export
#'
#' @examples
#' NULL
drug_intervention  <-  function( kill_prob = 0, block_prob = 1, gene,
                                 generate_mutations = TRUE ){

    local_environment( env = pck.env )

    id_gene  =  which( pck.env$onco$name == gene )
    if ( length( id_gene ) == 0 ) stop( 'There is no defined gene in the pool of clones with genes of interest')

    slct  =  which ( unlist( sapply( X = 1:length( pck.env$clones ),
                        FUN = function( x ) pck.env$clones[[ x ]]$gene[ id_gene ] ==  1 ) ) )
    if (length( slct ) > 0 ){
        rest  =   ( 1 : length( pck.env$clones ))[ - slct ]

        # Kill the cells by calc_binom( 1, N_cells, (1 - kill_prob) ) function:
        sapply( X = slct, FUN = function( x) pck.env$clones[[ x ]]$N_cells  =
                    calc_binom( 1, pck.env$clones[[ x ]]$N_cells, (1 - kill_prob) ) )

        # Save copies of clones (splitting the blocked clones):
        clones_to_change  =   lapply( X = slct, FUN = function( x ) clone_copy( pck.env$clones[[ x ]] ) )
        onco_clones_to_change  =   lapply( X = slct, FUN = function( x ) onco_copy( pck.env$onco_clones[[ x ]] ) )
        sapply( X = 1:length( slct ), FUN = function( x ) onco_clones_to_change[[ x ]]$id  =  clones_to_change[[ x ]]$id )

        # Number of blocking cells for each clone:
        N_block  =  unlist( sapply( X = slct, FUN = function( x) calc_binom( 1, pck.env$clones[[ x ]]$N_cells , block_prob ) ) )

        # Change numbers of cells for splitted clones
        sapply( X = 1:length( slct ), FUN = function( x ) {
                    pck.env$clones[[ slct[x] ]]$N_cells  = pck.env$clones[[ slct[x] ]]$N_cells - N_block[ x ]
        } )
        sapply( X = 1:length( slct ), FUN = function( x ) clones_to_change[[ x ]]$N_cells  =  N_block[ x ] )
        # transform last of them into passenger genes:
        sapply( X = 1:length( slct ), FUN = function( x ) {
                        clones_to_change[[ x ]]$gene[ id_gene ]  =  0
                        clones_to_change[[ x ]]$pasgene[ id_gene ]  =  1
                        clones_to_change[[ x ]]$invasion  =  FALSE
                        clones_to_change[[ x ]]$primary  =  FALSE
                        } )

        if ( generate_mutations ){
            for( x in 1:length( slct ) ){

                PointMut_ID  =  clones_to_change[[ x ]]$PointMut_ID
                if ( PointMut_ID[1] != 0 ){
                    for( i in PointMut_ID ){
                        if ( pck.env$pnt_clones[[ i ]]$Gene_name == gene &
                             pck.env$pnt_clones[[ i ]]$MalfunctionedByPointMut ){

                            # Generate NEW point mutation instead old one:
                            # 1. delete old index of point mutation
                            PointMut_ID  =  PointMut_ID[ ! PointMut_ID  %in%  i ]

                            # 2. Generate new point mutation with passenger gene and the same info:
                            pnt0  =  generate_pnt(    prntl =  pck.env$pnt_clones[[ i ]]$Parental_1or2,
                                                      gene  =  gene,
                                                      pos   =  pck.env$pnt_clones[[ i ]]$Ref_pos,
                                                      onco1 =  pck.env$onco,
                                                      Chr   =  pck.env$pnt_clones[[ i ]]$Chr,
                                                      mutation = NA )
                            pnt0$MalfunctionedByPointMut  =  FALSE
                            pnt0$Phys_pos  =  pck.env$pnt_clones[[ i ]]$Phys_pos
                            pnt0$Delta     =  pck.env$pnt_clones[[ i ]]$Delta
                            pnt0$Copy_number  =  pck.env$pnt_clones[[ i ]]$Copy_number

                            # 3. Add new number to a clone:
                            if ( length( PointMut_ID ) == 0  ) {
                                    PointMut_ID   =  pnt0$PointMut_ID
                            } else  PointMut_ID   =  c( PointMut_ID, pnt0$PointMut_ID )
                        }
                    }
                    clones_to_change[[ x ]]$PointMut_ID  =  PointMut_ID
                }

                CNA_ID       =  clones_to_change[[ x ]]$CNA_ID
                if ( CNA_ID[1] != 0 ){
                    for( j in CNA_ID ){
                        if ( gene %in% pck.env$cna_clones[[ j ]]$Gene_names &
                                       pck.env$cna_clones[[ j ]]$MalfunctionedByCNA ){

                            # Generate NEW CNA mutation instead old one:

                            # 1. delete old index of CNA mutation:
                            CNA_ID  =  CNA_ID[ ! CNA_ID  %in%  j ]

                            # 2. Generate new CNA mutation with passenger gene and the same info:
                            cna0  =  generate_cna( prntl  =  pck.env$cna_clones[[ j ]]$Parental_1or2,
                                                  genes  =  pck.env$cna_clones[[ j ]]$Gene_names,
                                                  start_end  =  c( pck.env$cna_clones[[ j ]]$Ref_start,
                                                                   pck.env$cna_clones[[ j ]]$Ref_end    ),
                                                  onco1  =  pck.env$onco,
                                                  dupOrdel  =  pck.env$cna_clones[[ j ]]$dupOrdel )
                            cna0$MalfunctionedByCNA  =  FALSE

                            # 3. Add new ID to a clone:
                            if ( length( CNA_ID ) == 0  ) {
                                CNA_ID   =  cna0$CNA_ID
                            } else  CNA_ID   =  c( CNA_ID, cna0$CNA_ID )
                        }
                    }
                    clones_to_change[[ x ]]$CNA_ID  =  CNA_ID
                }
            }
        }

        # Combine new clones with previous clones:
        pck.env$clones      =  c( pck.env$clones,            clones_to_change )
        pck.env$onco_clones =  c( pck.env$onco_clones,  onco_clones_to_change )
    }
    return( NULL )
}

