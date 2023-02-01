# Functions related to Point Mutations:  --------------------------------------

#' Function to save 1 point mutation in a data frame
#'
#' @param pnt Object of class 'Point_Mutations'
#'
#' @return data frame with 1 row of point mutation info
#' @export
#'
#' @examples
#' pnt = tugHall_dataset$pnt_clones[[ 1 ]]
#' df = safe_pnt_mut( pnt ) # View( pnt )
safe_pnt_mut  <-  function(pnt){
    return( pnt$safe() )
}

#' Function to copy of pnt1 without mutation info for allele A
#' @param pnt1 Object of class 'Point_Mutations'
#'
#' @return Object of class 'Point_Mutations' for another chromosome
#' @export
#'
#' @examples
#' pnt = tugHall_dataset$pnt_clones[[ 1 ]]
#' copy_pnt_no_mutation( pnt )  # View( pnt )
copy_pnt_no_mutation  <-  function( pnt1 ){
    ### Function to copy of pnt1 without mutation info for allele A
    pnt2  =  Point_Mutations$new()
    pnt2$PointMut_ID  =  pnt1$PointMut_ID
    pnt2$Allele = 'A'
    pnt2$Parental_1or2 = ifelse( pnt1$Parental_1or2 == 2, as.integer(1), as.integer(2) )
    pnt2$Chr  =  pnt1$Chr
    pnt2$Ref_pos  =  pnt1$Ref_pos
    pnt2$Phys_pos  =  NA
    pnt2$Delta     =  NA
    pnt2$Copy_number  =  1
    pnt2$Gene_name    =  pnt1$Gene_name
    pnt2$MalfunctionedByPointMut  =  NA
    pnt2$mut_order    =  pnt1$mut_order

    return( pnt2 )
}

#' Function to copy of point mutation info
#' @param pnt1 Object of class 'Point_Mutations'
#'
#' @return The same object of class 'Point_Mutations' with the same ID
#' @export
#'
#' @examples
#' pnt = tugHall_dataset$pnt_clones[[ 1 ]]
#' copy_pnt( pnt ) # View( pnt )
copy_pnt  <-  function( pnt1 ){
    ### Function to copy of pnt1
    pnt2  =  Point_Mutations$new()
    pnt2$PointMut_ID  =  pnt1$PointMut_ID
    pnt2$Allele       = pnt1$Allele
    pnt2$Parental_1or2  =  pnt1$Parental_1or2
    pnt2$Chr  =  pnt1$Chr
    pnt2$Ref_pos  =  pnt1$Ref_pos
    pnt2$Phys_pos  =  pnt1$Phys_pos
    pnt2$Delta     =  pnt1$Delta
    pnt2$Copy_number  =  pnt1$Copy_number
    pnt2$Gene_name    =  pnt1$Gene_name
    pnt2$MalfunctionedByPointMut  =  pnt1$MalfunctionedByPointMut
    pnt2$mut_order    =  pnt1$mut_order

    return( pnt2 )
}

#' Function to copy CNA info
#'
#' @param CNA1 Object of class 'CNA_Mutations'
#'
#' @return The same object of class 'CNA_Mutations'
#' @export
#'
#' @examples
#' cna = tugHall_dataset$cna_clones[[ 1 ]]
#' cna2 = copy_CNA( cna )
#' cna$safe()
#' cna2$safe()
copy_CNA  <-  function( CNA1 ){
    ### Function to copy of CNA1
    CNA2  =  CNA_Mutations$new()
    CNA2$CNA_ID  =  CNA1$CNA_ID
    CNA2$Parental_1or2  =  CNA1$Parental_1or2
    CNA2$dupOrdel = CNA1$dupOrdel
    CNA2$Chr  =  CNA1$Chr
    CNA2$Ref_start  =  CNA1$Ref_start
    CNA2$Ref_end    =  CNA1$Ref_end
    CNA2$Gene_names =  CNA1$Gene_names
    CNA2$MalfunctionedByCNA  =  CNA1$MalfunctionedByCNA

    CNA2$mut_order  =  CNA1$mut_order

    return( CNA2 )
}

#' Generation point mutation info
#'
#' @param onco1 Object of class 'OncoGene'
#' @param gm_1_2 List of two data frames (for 1st and 2nd parental chromosomes) with genes' location information
#'
#' @return list of (prntl - 1 or 2 parental chromosome, gene - gene name, pos - position of point mutation, Chr - name of chromosome )
#' @export
#'
#' @examples
#' gm = tugHall_dataset$gene_map
#' gm_1_2 = list( gm, gm )
#' onco = tugHall_dataset$onco
#' get_point_mutation( onco, gm_1_2 )
get_point_mutation <- function( onco1, gm_1_2 ){
    ### Function to get position of point mutation related to gene_map info for each chr
    prntl  =   sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))
    gene   =   ifelse( prntl == 1,
                       sample( onco1$name, size = 1, replace = TRUE, prob = onco1$cds_1/sum(onco1$cds_1) ),
                       sample( onco1$name, size = 1, replace = TRUE, prob = onco1$cds_2/sum(onco1$cds_2) ))

    gm     =   gm_1_2[[ prntl ]]

    ### get position for  gene !!! + CNA later!!!
    sp  =  which( gm$Gene == gene)   ### get the name of related gene
    cds =  ifelse( prntl == 1,
                   onco1$cds_1[ which(onco1$name == gene ) ],
                   onco1$cds_2[ which(onco1$name == gene ) ] ) ### get CDS for the related gene
    p   =  sample( sp, size = 1, replace = TRUE,
                   prob = gm[sp,'Len'] / sum( gm[sp,'Len'] ) )   ### get the block in gene_map
    pos =  sample( gm$Start[p]:gm$End[p], 1, replace=FALSE)
    Chr =  gm$Chr[p]

    return( list( prntl, gene, pos, Chr ) )
}

#' Generation point mutation info for the particular gene
#'
#' @param onco1 Object of class 'OncoGene'
#' @param gm_1_2 List of two data frames (for 1st and 2nd parental chromosomes) with genes' location information
#' @param gene Gene's name where point mutation should be occured
#'
#' @return list of (prntl - 1 or 2 parental chromosome, gene - gene name, pos - position of point mutation, Chr - name of chromosome )
#' @export
#'
#' @examples
#' gm = tugHall_dataset$gene_map
#' gm_1_2 = list( gm, gm )
#' onco = tugHall_dataset$onco
#' get_point_mutation_for_gene( onco, gm_1_2, gene = 'APC')
#' get_point_mutation_for_gene( onco, gm_1_2, gene = 'KRAS')
get_point_mutation_for_gene <- function( onco1, gm_1_2, gene ){
    ### Function to get position of point mutation related to gene_map info for each chr
    prntl  =   sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))

    gm     =   gm_1_2[[ prntl ]]

    ### get position for  gene !!! + CNA later!!!
    sp  =  which( gm$Gene == gene)   ### get the name of related gene
    cds =  ifelse( prntl == 1,
                   onco1$cds_1[ which(onco1$name == gene ) ],
                   onco1$cds_2[ which(onco1$name == gene ) ] ) ### get CDS for the related gene
    p   =  sample( sp, size = 1, replace = TRUE,
                   prob = gm[sp,'Len'] / sum( gm[sp,'Len'] ) )   ### get the block in gene_map
    pos =  sample( gm$Start[p]:gm$End[p], 1, replace=FALSE)
    Chr =  gm$Chr[p]

    return( list( prntl, gene, pos, Chr ) )
}


#' Function to generate an object of class 'Point_Mutations'
#'
#' @param prntl Parental chromosome, could be 1 or 2
#' @param gene Gene name
#' @param pos Position of point mutation
#' @param onco1 Object of class 'OncoGene'
#' @param Chr Chromosome name
#' @param mutation If mutation is NOT NA then MalfunctionedByPointMut = TRUE, else it is defined by corresponding probabilities
#'
#' @return Object of class 'Point_Mutations'
#' @export
#'
#' @examples
#' gm = tugHall_dataset$gene_map
#' gm_1_2 = list( gm, gm )
#' onco = tugHall_dataset$onco
#' pnt_clones = tugHall_dataset$pnt_clones
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' mut_order = 234  # As an example
#' withr::with_environment( env = pck.env, code = generate_pnt( prntl = 1, gene = 'APC', pos = 112767192, onco, Chr = '5', mutation = NA ) )
generate_pnt  <-  function( prntl, gene, pos, onco1, Chr, mutation = NA ) {
    ### Function to generate object of point mutation pnt
    pnt0 = Point_Mutations$new()
    pnt0$Gene_name = gene
    pnt0$PointMut_ID =  ifelse( length( pck.env$pnt_clones ) == 0, 1,
                                pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID + 2 )
    pnt0$Allele = 'B'  # Mutation occurs on allele B
    pnt0$Parental_1or2  =  as.integer(prntl)
    pnt0$Chr = Chr
    pnt0$Ref_pos  = pos
    pnt0$Phys_pos = pos
    pnt0$Delta = 0
    pnt0$Copy_number = 1
    u = ifelse( onco1$onsp[ which(onco1$name == gene) ] == 'o', uo, us)
    if ( is.na( mutation ) ) {
        pnt0$MalfunctionedByPointMut  =  ifelse( (runif(1) < u), TRUE, FALSE )
    } else {
        pnt0$MalfunctionedByPointMut  =  TRUE
    }
    pck.env$mut_order  =  pck.env$mut_order  +  1
    pnt0$mut_order  =  pck.env$mut_order

    pck.env$pnt_clones = c( pck.env$pnt_clones, pnt0 )

    pnt1  =  copy_pnt_no_mutation( pnt0 )
    pck.env$pnt_clones = c( pck.env$pnt_clones, pnt1 )

    return( pnt0 )
}

#' Function to generate the same object of class 'Point_Mutations' with coping all information from input object
#'
#' @param pnt Object of class 'Point_Mutations'
#'
#' @return The same object of class 'Point_Mutations' with different ID
#' @export
#'
#' @examples
#' pnt = tugHall_dataset$pnt_clones[[ 1 ]]
#' pnt_clones = tugHall_dataset$pnt_clones
#' pnt2 = generate_to_copy_pnt( pnt )
generate_to_copy_pnt  <-  function( pnt ) {
    ### Function to generate object of point mutation pnt with data from input pnt
    pnt0 = Point_Mutations$new()
    pnt0$Gene_name = pnt$Gene_name
    pnt0$PointMut_ID =  ifelse( length( pck.env$pnt_clones ) == 0, 1,
                                pck.env$pnt_clones[[ length( pck.env$pnt_clones ) ]]$PointMut_ID + 2 )
    pnt0$Allele = pnt$Allele
    pnt0$Parental_1or2  =  pnt$Parental_1or2
    pnt0$Chr = pnt$Chr
    pnt0$Ref_pos  = pnt$Ref_pos
    pnt0$Phys_pos = pnt$Phys_pos
    pnt0$Delta = pnt$Delta
    pnt0$Copy_number = pnt$Copy_number
    pnt0$MalfunctionedByPointMut  =  pnt$MalfunctionedByPointMut

    pnt0$mut_order  =  pnt$mut_order

    pck.env$pnt_clones  =  c( pck.env$pnt_clones, pnt0 )

    pnt1  =  copy_pnt_no_mutation( pnt0 )
    pck.env$pnt_clones  =  c( pck.env$pnt_clones, pnt1 )

    return( pnt0 )
}



# Functions related to CNA ------------------------------------------------


#' Function to choose probability of CNA mutation for several genes
#'
#' @param genes Names of genes, vector of names
#' @param dupOrdel It could be 'dup' or 'del' to denote duplication or deletion
#'
#' @return Single value of maximal probability from probabilities for several genes
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' onco = tugHall_dataset$onco
#' withr::with_environment( env = pck.env, code = get_u_cna( genes = 'APC', dupOrdel = 'dup' ) )
#' withr::with_environment( env = pck.env, code = get_u_cna( genes = c('KRAS','APC'), dupOrdel = c('dup', 'del') ) )
get_u_cna <- function( genes, dupOrdel ){
    # input: genes, u_del or u_dup for oncogenes and suppressors
    u = NULL
    for (gene in genes ) {
        os = pck.env$onco$onsp[ which( pck.env$onco$name == gene) ]
        u1 = ifelse( dupOrdel == 'dup',
                     ifelse( os == 'o', uo_dup, us_dup ),   # for duplication
                     ifelse( os == 'o', uo_del, us_del ) )  # for deletion
        u = c( u, u1 )
    }

    return( max( u ) )
}


#' Generation CNA mutation info
#'
#' @param onco1 Object of class 'OncoGene'
#' @param dupOrdel It could be 'dup' or 'del' to denote duplication or deletion
#' @param gm_1_2 List of two data frames (for 1st and 2nd parental chromosomes) with genes' location information
#'
#' @return List of ( prntl - 1 or 2 parental chromosome, Chr - name of chromosome, genes - genes names, start_end - vector with start and end positions of CNA, w_cna - rows of CNA in gene_map data frame )
#'
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' onco = tugHall_dataset$onco
#' gm = tugHall_dataset$gene_map
#' withr::with_environment( env = pck.env, code = get_cna_mutation( onco1 = onco, dupOrdel = 'dup', gm_1_2 = list(gm, gm) ) )
#' withr::with_environment( env = pck.env, code = get_cna_mutation( onco1 = onco, dupOrdel = 'del', gm_1_2 = list(gm, gm) ) )
get_cna_mutation <- function( onco1, dupOrdel, gm_1_2 ){
    ### Function to get position of point mutation related to gm (gene_map) info for each chr
    prntl  =  sample( c(1,2), size = 1, replace = TRUE, prob = c( sum(onco1$cds_1), sum(onco1$cds_2) ))
    lambda =  ifelse( dupOrdel == 'dup', lambda_dup, lambda_del )
    l_cna  =  round( rexp(1, 1/lambda) ) +1 # rpois( n = 1, lambda ) + 1   # length of CNA

    gm     =  gm_1_2[[ prntl ]]

    w1  =  sample( 1:length( gm$Start ), size = 1, prob = gm$Len )  # choose the row in gene_map
    Chr = gm$Chr[ w1 ]
    w   =  which( gm$Chr == Chr )    ### Area of checking for CNA
    pos_st =  sample( gm$Start[w1]:gm$End[w1], 1, replace=FALSE)   ### position of start CNA
    pos_end  = pos_st  +  l_cna

    w3  =  which( gm$Start < pos_end & gm$Chr == Chr )
    w_cna  = w3[ which(w3 >= w1)]    ### rows of CNA in gene_map
    start_end  =  c(0,0)
    start_end[1]  =  pos_st
    start_end[2]  =  ifelse( pos_end < gm$End[ max(w_cna) ], pos_end, gm$End[ max(w_cna) ] )
    genes   =   unique( gm$Gene[w_cna] )   ### Genes of CNA at same Chr

    return( list( prntl, Chr, genes, start_end, w_cna ) )
}

#' Function to get order of mutation for all possible types
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return data.frame with fields order, type, ID
#' @export
#'
#' @examples
#' clone = tugHall_dataset$clones[[ 46 ]]
#' clone$PointMut_ID
#' clone$CNA_ID
#' pnt_clones = tugHall_dataset$pnt_clones
#' cna_clones = tugHall_dataset$cna_clones
#' mixed_mut_order( clone )
mixed_mut_order   <-   function( clone1 ) {
    # return the order, type and number of mutation (del, dup or point in order of appearance)
    order_tp_num  <-  data.frame( order = NULL, type = NULL, ID = NULL )
    i  =  0
    if ( length( clone1$PointMut_ID ) > 0  & clone1$PointMut_ID[ 1 ]  != 0 ){
        for (i in 1:length( clone1$PointMut_ID )) {
            order_tp_num[i,'type']   =  'pnt'
            order_tp_num[i,'ID']     =  as.numeric( clone1$PointMut_ID[i] )
            order_tp_num[i,'order']  =  pck.env$pnt_clones[[ order_tp_num[i,'ID'] ]]$mut_order   ### pn[order_tp_num[i,'ID'], 'mut_order']
        }
    }

    if ( length( clone1$CNA_ID ) > 0 & clone1$CNA_ID[ 1 ]  != 0 ){
        for (j in 1:length( clone1$CNA_ID ) ) {
            order_tp_num[j+i,'ID']      =   clone1$CNA_ID[ j ]
            cn1   =   pck.env$cna_clones[[ order_tp_num[ j+i, 'ID' ] ]]
            order_tp_num[j+i,'order']   =   cn1$mut_order ###  cn[order_tp_num[ j+i, 'ID' ], 'mut_order']
            order_tp_num[j+i,'type']    =   cn1$dupOrdel  ###  as.character( cn[order_tp_num[ j+i, 'ID' ], 'dupOrdel'] )
        }
    }

    if ( length(order_tp_num)  !=  0 ){
        order_tp_num  <-  order_tp_num[ order( order_tp_num$order), ]
        row.names(order_tp_num) <- 1:length( order_tp_num[,'type'])
    } else order_tp_num = NULL

    return( order_tp_num )
}

#'  Function to generate object of CNA mutation
#'
#'
#' @param prntl The 1st or 2nd parental chromosome
#' @param genes Genes names
#' @param start_end vector with start and final positions of CNA
#' @param onco1 Object of class 'OncoGene'
#' @param dupOrdel It could be 'dup' or 'del' to denote duplication or deletion
#'
#' @return Object of class 'CNA_Mutations'
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' load_tugHall.Environment( tugHall_dataset )
#' onco = tugHall_dataset$onco
#' gene_map = tugHall_dataset$gene_map
#' pnt_clones = tugHall_dataset$pnt_clones
#' cna_clones = tugHall_dataset$cna_clones
#' mut_order = 234
#' start_end = c(112775658, 112775716 )
#' withr::with_environment( env = pck.env, code = generate_cna( prntl = 1, genes = 'APC', start_end = start_end, onco1 = onco, dupOrdel = 'dup' ) )
generate_cna  <-  function( prntl, genes, start_end, onco1, dupOrdel ) {
    ### Function to generate object of point mutation cna
    cna0 = CNA_Mutations$new()
    cna0$CNA_ID =  ifelse( length( pck.env$cna_clones ) == 0, 1,
                           pck.env$cna_clones[[ length( pck.env$cna_clones ) ]]$CNA_ID + 1 )
    cna0$Parental_1or2  =  prntl
    cna0$dupOrdel = dupOrdel
    cna0$Chr = gene_map$Chr[ which( gene_map$Gene == genes[1] )[1] ]
    cna0$Ref_start  = start_end[1]
    cna0$Ref_end    = start_end[2]
    cna0$Gene_names = genes
    u = get_u_cna( genes, dupOrdel )
    cna0$MalfunctionedByCNA  =  ifelse( ( runif(1) < u ), TRUE, FALSE )

    pck.env$mut_order  =  pck.env$mut_order  +  1
    cna0$mut_order  =  pck.env$mut_order

    pck.env$cna_clones = c( pck.env$cna_clones, cna0 )

    return( cna0 )
}

#' Function to add the mutations to the data.frame gene_map
#'
#' @param clone1 Object of class 'Clone'
#' @param onco1 Object of class 'OncoGene'
#'
#' @return list( gm1, gm2 ), where gm1 and gm2 are data.frames gene_maps with mutation information
#' @export
#'
#' @examples
#' clone = tugHall_dataset$clones[[ 46 ]]
#' onco = tugHall_dataset$onco_clones[[ 46 ]]
#' gene_map = tugHall_dataset$gene_map
#' pnt_clones = tugHall_dataset$pnt_clones
#' cna_clones = tugHall_dataset$cna_clones
#' gene_map$pnts = ''
#' \dontrun{
#' gm_1_2 = modify_gene_map( clone , onco )  # View(gm_1_2)
#' }
modify_gene_map  <-  function( clone1 , onco1 ){
    # index    =   which( sapply( onco_clones , FUN = function(x)  x$id == clone1$id) )  # index of clone1 in clones list
    # onco1  =  onco_clones[[ index ]]
    gm1    =  gm2  =  gene_map        # locations for the first and second chromosomes
    gm1$pnts  =  gm2$pnts  = ''
    mixed_order  =  mixed_mut_order( clone1 = clone1 )
    if ( is.null(mixed_order) ) return( list( gm1, gm2 ) )

    for (l in 1:nrow( mixed_order ) ) {
        cn1  =  NULL
        if ( mixed_order[l, 'type'] == 'del' ) {

            cn1        =  pck.env$cna_clones[[  mixed_order$ID[l] ]]   # cn[ mixed_order$ID[l], ]
            Ref_start  =  cn1$Ref_start
            Ref_end    =  cn1$Ref_end
            Chr        =  cn1$Chr

            # change gene map only one of two chromosomes: 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 )

            gm  =  add_deletion( gm = gm, Ref_start = Ref_start,
                                 Ref_end = Ref_end, Chr = Chr )

            ### come back to gene map for certain chromosome 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
            rm( gm )
        }

        if ( mixed_order[l, 'type'] == 'dup' ) {

            cn1        =  pck.env$cna_clones[[  mixed_order$ID[l] ]]   ###  cn[ mixed_order$ID[l], ]
            Ref_start  =  cn1$Ref_start
            Ref_end    =  cn1$Ref_end
            Chr        =  cn1$Chr

            # change gene map only one of two chromosomes: 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 )

            gm  =  add_duplication( gm = gm, Ref_start = Ref_start,
                                    Ref_end = Ref_end, Chr = Chr )

            ### come back to gene map for the certain chromosome 1 or 2
            ifelse( cn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
            rm( gm )
        }

        if ( mixed_order[l, 'type'] == 'pnt' ) {

            pn1        =  pck.env$pnt_clones[[ mixed_order$ID[l] ]]   ###  pn[ mixed_order$ID[l], ]
            pos_pnt1   =  pn1$Ref_pos
            Chr        =  pn1$Chr

            # change gene map only one of two chromosomes: 1 or 2
            ifelse( pn1$Parental_1or2 == 1, gm <- gm1, gm <- gm2 )

            gm  =  add_pnt_mutation( gm = gm, pos_pnt = pos_pnt1, Chr = Chr )

            ### come back to gene map for the certain chromosome 1 or 2
            ifelse( pn1$Parental_1or2 == 1, gm1 <- gm, gm2 <- gm )
            rm( gm )

        }

    }

    return( list( gm1, gm2 ) )
}

#' Function to add point mutation to data.frame gene_map (chromosomal location data frame)
#'
#' @param gm Chromosomal location data frame
#' @param pos_pnt Position of point mutation
#' @param Chr Chromosome name
#'
#' @return Chromosomal location data frame with additional point mutation info
#' @export
#'
#' @examples
#' gene_map = tugHall_dataset$gene_map
#' gene_map$pnts  =  ''
#' gm2 = add_pnt_mutation( gm = gene_map, pos_pnt = 112775637 , Chr = '5' )
add_pnt_mutation   <-  function( gm = gm, pos_pnt, Chr = Chr ){

    w = which( gm$Chr == Chr & gm$Start <= pos_pnt  & gm$End >= pos_pnt )

    if ( length(w) != 1 ) {
        print( w )
        print( pos_pnt )
        print( gm )
        stop( 'Point mutation does not meet into location ranges' )
    } else {
        gm[w, 'pnts']  =  ifelse( gm[w, 'pnts'] ==  '',
                                  gm[w, 'pnts']  <-  as.character(pos_pnt),
                                  gm[w, 'pnts']  <-  paste(gm[w, 'pnts'], pos_pnt, sep = ',') )
    }

    return( gm )
}

#' Function to add deletion to gene map (chromosomal location data frame)
#'
#' @param gm Chromosomal location data frame
#' @param Ref_start Starting position of deletion
#' @param Ref_end Final position of deletion
#' @param Chr Chromosome name
#'
#' @return Chromosomal location data frame with additional deletion info
#' @export
#'
#' @examples
#' gene_map  =  tugHall_dataset$gene_map
#' gene_map$pnts = ''
#' gm_new = add_deletion( gm = gene_map, Ref_start = 112775658, Ref_end = 112775716, Chr = '5')
#' # to check add_del ... add_dupl, please, use diff library and:
#' # diff_data( gm_ref, gm, show_unchanged_columns = TRUE, always_show_order = TRUE )
add_deletion  <-  function( gm, Ref_start, Ref_end, Chr ){
    ### Change gene_map with CNA deletion:
    if (Ref_end < Ref_start)  stop( 'End should be larger or equal then Start for CNA' )
    if (Ref_start < min(gm$Start[gm$Chr == Chr]) ) stop( 'The CNA should be inside of the cromosomal locations' )
    w1  =  max( which( gm$Start <= Ref_start  &  gm$Chr == Chr ) )
    w2  =  max( which( gm$Start <= Ref_end    &  gm$Chr == Chr ) )

    if ( w1 == w2 ) {
        gm  =  rbind( gm[1:w1,], gm[ w1:nrow(gm), ] )  # duplicate the row w1=w2
        w2  =  w1 + 1
    }
    ### Correction of rows w1 and w2:
    if ( gm[ w1, 'Start']   <   Ref_start ){  # to delete part of row or leave as before
        gm[ w1, 'End']    =   ifelse( Ref_start <=  gm[ w1, 'End'], Ref_start - 1, gm[ w1, 'End']  )
        gm[ w1, 'Len']    =   gm[ w1, 'End']  -  gm[ w1, 'Start']    +  1
        gm[ w1, 'pnts']   =   check_pnts( gm[ w1,  ] )
    } else {
        w1  =  w1  -  1   # to delete whole row
    }

    if ( gm[ w2, 'End']   >   Ref_end ){   # to delete part of row or leave as before
        gm[ w2, 'Start']    =   Ref_end + 1
        gm[ w2, 'Len']      =   gm[ w2, 'End']  -  gm[ w2, 'Start']  +  1
        gm[ w2, 'pnts']     =   check_pnts( gm[ w2,  ] )
    } else {
        w2  =  w2  +  1    # to delete whole row
    }

    ### delete part of gm related to deletion
    if ( w2 - w1 > 1 )      gm   =    gm[ -c( (w1+1):(w2-1) ), ]

    # To subtract the delta from positions
    w    =  which( gm$Start > Ref_end & gm$Chr == Chr )
    dlt  =  Ref_end  -  Ref_start  +  1
    if ( length(w) > 0 ){
        gm[w, 'Start']    =  gm[w, 'Start']  -  dlt
        gm[w, 'End']      =  gm[w, 'End']    -  dlt
        gm[w, 'pnts']     =  unlist( sapply( w, FUN = function(x)  pnts_add_dlt( gm[ x, ], dlt  =  -dlt)  ) )
    }

    return( gm )
}

#' Function to add duplication to gene map (chromosomal location data frame)
#'
#' @param gm Chromosomal location data frame
#' @param Ref_start Starting position of duplication
#' @param Ref_end Final position of duplication
#' @param Chr Chromosome name
#'
#' @return Chromosomal location data frame with additional duplication info
#' @export
#'
#' @examples
#' gene_map  =  tugHall_dataset$gene_map
#' gene_map$pnts = ''
#' gm_new = add_duplication( gm = gene_map, Ref_start = 112775658, Ref_end = 112775716, Chr = '5')
#' # to check add_del ... add_dupl, please, use diff library and:
#' # diff_data( gm_ref, gm, show_unchanged_columns = TRUE, always_show_order = TRUE )
add_duplication  <-  function( gm, Ref_start, Ref_end, Chr ){

    if (Ref_end < Ref_start)  stop( 'End should be larger or equal then Start for CNA' )
    if (Ref_start < min(gm$Start[gm$Chr == Chr]) ) stop( 'The CNA should be inside of the cromosomal locations' )
    dlt        =  Ref_end  -  Ref_start  +  1  # delta for all next chromosomal positions

    ### Change gene_map with CNA duplication:
    w1  =  max( which( gm$Start <= Ref_start  &  gm$Chr == Chr ) )
    w2  =  max( which( gm$Start <= Ref_end    &  gm$Chr == Chr ) )

    # ADD rows for duplication
    if ( (w2+1) <= nrow(gm) ){
        gm  =  rbind( gm[1:w2, ], gm[w1:w2, ] , gm[(w2+1):nrow(gm), ])
    } else {
        gm  =  rbind( gm[1:w2, ], gm[w1:w2, ] )
    }

    w3  =  w2 + 1

    w5  =  max( which( gm$Chr  ==  Chr ) )
    # change the final position BEFORE duplication
    gm[w2, 'End']    =  ifelse( gm[w2, 'End'] > Ref_end, Ref_end, gm[w2, 'End'] )
    gm[w2, 'Len']    =  gm[w2, 'End'] - gm[w2, 'Start'] + 1
    gm[w2, 'pnts']   =  check_pnts( gm[ w2,  ] )

    # Add delta to all positions after duplication
    gm[w3:w5, 'Start']  =  gm[w3:w5, 'Start'] + dlt
    gm[w3:w5, 'End']    =  gm[w3:w5, 'End']   + dlt
    gm[w3:w5, 'pnts']   =  unlist( sapply( w3:w5, FUN = function(x) pnts_add_dlt( gm[ x, ], dlt )  ) )

    # Correct Start of w3 row if it's necessary
    if ( gm[w3, 'End']  >= (Ref_start + dlt) ){
        gm[w3, 'Start']  =  Ref_start + dlt
        gm[w3, 'Len']    =  gm[w3, 'End'] - gm[w3, 'Start'] + 1
        gm[w3, 'pnts']   =  check_pnts( gm[ w3,  ] )
    } else { # delete the row if duplication started from outside of gene positions
        gm  =  gm[-w3,]
    }

    return( gm )
}

## to check add_del ... add_dupl, please, use diff library and:
## diff_data( gm_ref, gm, show_unchanged_columns = TRUE, always_show_order = TRUE )

#' Function to subtract delta from position of point mutations
#' @param gm_w1 A row from data.frame gene_map
#' @param dlt Delta to subtract from positions of point mutations
#'
#' @return Return the pnts - dlt for one row of data.frame gene_map
#' @export
#'
#' @examples
#' gene_map = tugHall_dataset$gene_map
#' gene_map$pnts = ''
#' gene_map$pnts[6] = '112792451'
#' gm_w1 = gene_map[6,]
#' pnts_add_dlt( gm_w1 , dlt = 1000 )
#' pnts_add_dlt( gm_w1 , dlt = -1001 )
pnts_add_dlt  <-  function( gm_w1 , dlt ){
    ### Return the pnts - dlt for one row of gene map
    if ( is.null(gm_w1) )  stop( 'The input is null' )
    pnts = gm_w1$pnts
    if ( length(pnts)  ==  0 )  stop( 'The input is empty' )
    if ( pnts == '' ) return( '' )
    pnts  =  as.numeric( unlist( strsplit( pnts, split = ',') ) )
    if ( !is.numeric( pnts ) ) stop( 'Incorrect format of points mutation: should be numeric')
    pnts  =  pnts  +  dlt
    pnts  =  paste( as.character( pnts ) , collapse = ',')

    return( pnts )
}

#' Function to check what pnts do fall into the range?
#' @param gm_w1 A row from data.frame gene_map
#'
#' @return Return the point mutations which fall into the range
#' @export
#'
#' @examples
#' gene_map = tugHall_dataset$gene_map
#' gene_map$pnts = ''
#' gene_map$pnts[6] = '112792451, 112792442'
#' gm_w1 = gene_map[6,]
#' check_pnts( gm_w1 )
check_pnts  <-  function( gm_w1 ){
    ### Return the pnts which fall into the range
    if ( is.null(gm_w1) )  stop( 'The input is null' )
    pnts = gm_w1$pnts
    if ( length(pnts)  ==  0 )  stop( 'The input is empty' )
    if ( pnts == '' ) return( '' )
    pnts  =  as.numeric( unlist( strsplit( pnts, split = ',') ) )
    if ( !is.numeric( pnts ) ) stop( 'Incorrect format of points mutation: should be numeric')

    check_log  =  unlist( sapply( pnts, FUN = function( x ) ( x <= gm_w1$End & x >= gm_w1$Start ) ) )
    pnts  =  paste( as.character( pnts[ check_log ] ) , collapse = ',')
    pnts  =  ifelse( length(pnts) == 0, '', pnts )

    return( pnts )
}
## for tests:
## Prepare the input and then make test:
## data.frame( test = sapply( X = 1:12 , FUN = function( x ) check_pnts( gm_w1 = gm1[ x, ]) ) )
## data.frame( test = sapply(1:12 , FUN = function( x ) pnts_add_dlt( gm1[ x, ], dlt  =  -dlt ) ) )


#' Function to get length of CDS and of genes from data.frame gene_map and related probabilities
#' @param gm data.frame gene_map with info about genes' location
#'
#' @return list( names, CDS, RNA, PROB, SUM, P0 )
#' @export
#'
#' @examples
#' gene_map  =  tugHall_dataset$gene_map
#' load_tugHall.Environment( tugHall_dataset )
#' withr::with_environment( env = pck.env, code = get_cds_rna( gm = gene_map ) )
get_cds_rna  <-  function( gm ){

    name0  =  unique( gm$Gene )
    rna0  =  cds0  =  NULL
    for ( i in 1:length(name0) ) {
        w     =  which( gm$Gene == name0[ i ] )
        cds0  =  c( cds0, sum( gm$End[w]  -  gm$Start[w] + 1 ) )
        rna0  =  c( rna0, max( gm$End[w] ) - min( gm$Start[w]) + 1  )
    }

    sum_prob   =  sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
    prob       =  c(   m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob
    p0         =   (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)

    return( list( names = name0, CDS = cds0, RNA = rna0, PROB = prob, SUM = sum_prob, P0 = p0 ) )

}
