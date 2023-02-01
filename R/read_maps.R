
### The code to read maps of genes using CCDS code of genes  ###
# library(stringr)   # to use string data in input files

#' Function to make a gene_map data.frame with information of genes' locations
#'
#' @param f_out Name of file to save gene_map data.frame
#' @param ls List of IDs of genes corresponding CCDS database https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/
#' @param f_in Name of file to input downloaded from CCDS database
#'
#' @return gene_map data.frame with information of genes' locations for genes of interest
#' @export
#'
#' @examples
#' url = 'https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt'
#' download.file( url = url, destfile = 'CCDS.current.txt')
#' ls   =  c( 'CCDS4107.1', 'CCDS8702.1', 'CCDS43171.1', 'CCDS11118.1' )
#' gene_map = make_map(f_out    =  'map.txt', ls = ls, f_in =  'CCDS.current.txt' )
make_map  <-  function(f_out    =  'Input/map.txt',
                       ls   =  c( 'CCDS4107.1', 'CCDS8702.1',
                                  'CCDS43171.1', 'CCDS11118.1' ),
                       f_in =  'Input/CCDS.current.txt'
                       ){
    ### ls - the list of CCDS ids of genes for simulation

    ### Read of the file with genome information (get dataset)
    df  <-  read.delim( file = f_in, sep = '\t', header = TRUE)

    ### Extract the row ids of dataset
    ids  <-  NULL
    for (i in 1:length( ls ) )        ids  <-  c( ids, which( df$ccds_id == ls[i] ) )

    ### df[ ids, ]

    ###  Get the only needed data

    dt  <- df[ ids, ]
    rm( df, i, ids )

    ### get a vector of strings:
    dfALL  <-  NULL
    for (i in 1:length(dt[,1]) ) {
        lct  <-   as.character( dt[ i, 'cds_locations'] )

        lct  <-  str_sub( lct, 2, str_length( lct )-1 )
        lct  <-   strsplit( lct , split = ',')
        lct  <-  lct[[1]]

        lct  <-  strsplit( lct[1:length(lct)] , split = '-')


        df <- data.frame( matrix( unlist(lct), nrow=length(lct), byrow=TRUE ),
                         stringsAsFactors=FALSE )
        df[ , 1] <- as.numeric( df[ , 1] )
        df[ , 2] <- as.numeric( df[ , 2] )
        names( df )  <-  c( 'Start', 'End' )

        dfALL  <-  rbind( dfALL,
                         data.frame( Chr = as.character( dt$X.chromosome[i] ),
                                     CCDS_ID = as.character( dt$ccds_id[i] ),
                                     Gene = as.character( dt$gene[i] ),
                                     df, stringsAsFactors = FALSE )
                        )
    }

    dfALL$Len  <-  dfALL$End  -  dfALL$Start + 1

    write.table( dfALL, file = f_out, quote = FALSE, sep = '\t',
                 row.names = FALSE, col.names = TRUE )

    return( dfALL )
}


#' Function to order info in gene_map data.frame with information of genes' locations
#'
#' @param gene_map data.frame with information of genes' locations
#'
#' @return The same data.frame gene_map with ordered positions for each gene and each chromosome
#' @export
#'
#' @examples
#' gene_map = tugHall_dataset$gene_map
#' gene_map = order_gene_map( gene_map )
order_gene_map  <-  function( gene_map ){
    gm = NULL
    for( t in   unique(gene_map$Chr) ) {

        w = which( gene_map$Chr == t )
        gm2  =  gene_map[ w, ]
        gmt  =  gm2[ order( gm2$Start ), ]
        gm   =  rbind( gm, gmt )

    }
    rownames(gm) <- 1:length(gm$Chr)
    return( gm )
}

#' Function to get length of CDS and whole gene from gene_map data.frame
#'
#' @param gene_map data.frame with info about genes' locations
#'
#' @return list of ( Name, CDS, LEN_Genes ) where Name is a vector of genes' names, CDS is a vector of CDS lengths, LEN_Genes is a vector of length of whole genes including introns and exons
#' @export
#'
#' @examples
#' gene_map = tugHall_dataset$gene_map
#' onco = tugHall_dataset$onco
#' get_len_cds_rna( gene_map)
get_len_cds_rna  <-  function( gene_map){

    name0  =  unique( gene_map$Gene )
    cds0   =  NULL
    rna0   =  NULL
    # Get length of CDS and RNA from gene_map:
    for ( i in 1:length(name0) ) {
        w     =  which( gene_map$Gene == name0[ i ] )
        cds0  =  c( cds0, sum( gene_map$End[w]  -  gene_map$Start[w] + 1 ) )
        rna0  =  c( rna0, max( gene_map$End[w] ) - min( gene_map$Start[w]) + 1  )
    }

    w = unlist( sapply( name0, FUN = function(x) which( pck.env$onco$name  ==  x ) ) )

    name0  =  name0[w]
    rna0   =  rna0[w]
    cds0   =  cds0[w]

    return( list( Name = name0, CDS = cds0, LEN_Genes = rna0 ) )
}



