
#' Function to get data about last simulation from cloneoutfile
#'
#' @param cloneoutfile Name of file to read data about clone evolution
#' @param genefile Name of file with hallmarks values
#' @param mainDir Working directory, by default mainDir = getwd()
#' @param sbdr_Output Directory for output data getting from mainDir
#'
#' @return list of data.frames like onco, hall, data_last (data of last time step), \cr
#' data_avg (average data for all time steps), data_flow (data without average rows), time_max (max time step), \cr
#' pnt_mut and pnt_mut_B (data.frame of point mutations for both alleles and for allele B only ) \cr
#' and cna_mut (data.frame of CNA mutations)
#' @export
#' @examples
#' copy_files_to_Input()
#' load_tugHall.Environment( results = tugHall_dataset )
#' copy_files_to_Output()
#' withr::with_environment( env = pck.env, code = { dataset = get_flow_data(cloneoutfile, genefile, mainDir = getwd(), sbdr_Output = '/Output' ) } )
#' # View(dataset)
get_flow_data <- function(cloneoutfile, genefile, mainDir = getwd(), sbdr_Output = '/Output' ) {

    # Get data about onco and hallmarks
    onco   =  pck.env$onco  #  oncogene$new()        # make the vector onco about the hallmarks
    # onco$read( genefile )          # read the input info to the onco from genefile - 'gene_cds2.txt'
    hall   =  pck.env$hall    #  hallmark$new()        # make a vector hall with hallmarks parameters
    #  hall$read( genefile, onco$name )     # read from the genefile - 'gene_cds2.txt'

    ## Get data from output file
    data_out  <-  read_file( file_name = cloneoutfile )   # read.table( cloneoutfile, sep="\t", header = TRUE )
    data_out[is.na(data_out)]  <-  ""

    # average data
    data_avg  =  data_out[ which( data_out$AvgOrIndx  ==  "avg"), ]

    # data without averaging - flow data
    data_flow  =  data_out[ which( !data_out$AvgOrIndx == "avg" ), ]

    data_flow[ , c( 1:22, 24, 30:41 ) ]  =  sapply( X = c( 1:22, 24, 30:41 ),
              FUN = function( x ) {
                  as.numeric( data_flow[ , x ] )
    } )


    # the data of the last time step
    time_max   =  max( data_flow$Time )
    data_last  =  data_flow[ which( data_flow$Time  ==  time_max ), ]
    rm( data_out )
    dr  =  file.path( mainDir, sbdr_Output )
    cna_mut  =  read_file( file_name = paste0( dr, '/CNA_mutations.txt' ) )
      #  read.table( file = 'Output/CNA_mutations.txt', sep = '\t', header = TRUE )
    pnt_mut  =  read_file( file_name = paste0( dr, '/point_mutations.txt' ) )
      #  read.table( file = 'Output/point_mutations.txt', sep = '\t', header = TRUE )

    return( list( onco = onco, hall = hall, data_avg = data_avg,
                  data_flow = data_flow, time_max = time_max, data_last = data_last,
                  cna_mut = cna_mut, pnt_mut = pnt_mut,
                  pnt_mut_B = pnt_mut[ which( !is.na( pnt_mut$MalfunctionedByPointMut ) ) , ]
                  ) )
}

#' Function to plot main data from data.frame with average data
#'
#' @param data_avg data.frame with average values from cloneout.txt file
#' @param time_max Maximal time step in a simulation
#'
#' @return NULL, draw many plot with average data
#' @export
#'
#' @examples
#' data_avg = tugHall_dataset$data_avg
#' time_max = tugHall_dataset$time_max
#' plot_average_simulation_data( data_avg , time_max = time_max )
plot_average_simulation_data  <-  function( data_avg, time_max ){
    # plot number of cells:
    sapply( X = c('grDevices', 'graphics' ), FUN = function( x ) check_pkg( x ) )

    flnm =  readline(prompt=" Press Enter to skip saving to a file or enter name of file to save data ")
    if ( flnm != '' ){
        safe_pdf = TRUE
    } else safe_pdf = FALSE

    g_range_y  =  range( 0, data_avg$N_normal + 1 , data_avg$N_primary + 1, data_avg$N_metastatic + 1 )
    g_range_x  =  range( min( data_avg$Time ), max( data_avg$Time, round( time_max/10 +0.5 ) * 10  ) )
    plot_2D_lines( x = data_avg$Time, data_avg,
                   nl = c('N_normal', 'N_primary', 'N_metastatic'), xr = g_range_x, yr = g_range_y,
                   legend_names = c('Normal', 'Primary tumor', 'Metastatic'),
                   names = c( 'Time step', 'Number of cells' ) ,
                   safe_pdf  =  safe_pdf,
                   filename = paste0( flnm, '_NM.pdf' )  )
    rl =  readline(prompt="This is a plot for Numbers of primary tumor and Metastasis cells - Press Enter  ")


    g_range_y  =  range( 0, 1 )
    plot_2D_lines( x = data_avg$Time, data_avg,
                   nl = c('a', 'd', 'i', 'im', 'k'), xr = g_range_x, yr = g_range_y,
                   legend_names = '',
                   names = c( 'Time step', 'Probabilities' ) ,
                   safe_pdf  =  safe_pdf,
                   filename = paste0( flnm, '_probabilities.pdf' )  )
    rl =  readline(prompt="This is a plot for average probabilities - Press Enter  ")


    g_range_y  =  range( 0, max( data_avg[c('Ha', 'Hd', 'Hi', 'Him', 'Hb')] ) + 0.04 )
    plot_2D_lines( x = data_avg$Time, data_avg,
                   nl = c('Ha', 'Hd', 'Hi', 'Him', 'Hb'), xr = g_range_x, yr = g_range_y,
                   legend_names = '',
                   names = c( 'Time step', 'Hallmarks' ) ,
                   safe_pdf  =  safe_pdf,
                   filename = paste0( flnm, '_hallmarks.pdf' )  )

    rl =  readline(prompt="This is a plot for average values of Hallmarks - Press Enter")


    g_range_y  =  range( 0, max( data_avg[c('mut_den')] ) + 0.05 )
    plot_2D_lines( x = data_avg$Time, data_avg,
                   nl = c('mut_den'), xr = g_range_x, yr = g_range_y,
                   legend_names = '',
                   names = c( 'Time step', 'Mutation rate' ) ,
                   safe_pdf  =  safe_pdf,
                   filename = paste0( flnm, '_mutation_rate.pdf') ,
                   draw_key = FALSE )

    rl =  readline(prompt="This is a plot for average value of mutation rate - Press Enter")
}


#' Function to get order of genes' dysfunction
#'
#' @param pnt_mut data.frame with info about all the point mutations
#' @param file_name Name of file to save data
#' @param data_last data.frame with data of simulation at the last time step
#' @param cna_mut data.frame with info about all the CNA mutations
#'
#' @return data.frame of genes' dysfunction and save it in a file
#' @export
#'
#' @examples
#' load_tugHall.Environment( results = tugHall_dataset )
#' copy_files_to_Output()
#' dtst = get_flow_data( pck.env$cloneoutfile, pck.env$genefile )
#' pnt_mut   =  dtst$pnt_mut
#' data_last  =  dtst$data_last
#' cna_mut = dtst$cna_mut
#' file_name = './Output/order_genes_dysfunction.txt'
#' rdr = get_order_of_genes_dysfunction( pnt_mut, data_last, cna_mut, file_name = file_name )
get_order_of_genes_dysfunction  <-  function( pnt_mut, data_last, cna_mut, file_name = './Output/order_genes_dysfunction.txt' ){

    sapply( X = c('grDevices', 'graphics' ), FUN = function( x ) check_pkg( x ) )

    # as.numeric( unlist( str_split( data_last[ , 'PointMut_ID' ], pattern = ',' ) ) )

    # sum( as.numeric( unlist( str_split( data_last[ 1, 'driver_genes' ], pattern = ' ' ) ) ) ) != 0

    ch = sapply(X = 1:nrow( data_last ), FUN = function( x )  sum( as.numeric( unlist( str_split( data_last[ x , 'driver_genes' ], pattern = ' ' ) ) ) ) != 0 )
    genes_dysfunction  =  data_last[ ch, c('N_cells', 'ID', 'ParentID', 'Birth_time',
                                           'type', 'mut_den', 'driver_genes', 'passenger_genes',
                                           'PointMut_ID', 'CNA_ID' )  ]
    if ( nrow( genes_dysfunction ) > 0  ) {
        for (i in 1:nrow( genes_dysfunction ) ){

            drivers  =  as.numeric( unlist( str_split( genes_dysfunction$driver_genes[ i ] , pattern = ' ') )  )
            # nm       =  onco$name[ as.logical(drivers) ]

            ### Point mutations:
            id_pnt_mut    =  as.numeric( unlist( str_split( genes_dysfunction$PointMut_ID[ i ] , pattern = ',') ) )
            if ( id_pnt_mut[1] != 0 ){
                r = as.vector( sapply( X = id_pnt_mut, FUN = function( x ) which( pnt_mut$PointMut_ID ==   x  ) ) )

                l = pnt_mut[ r , 'MalfunctionedByPointMut' ]
                l[ is.na( l ) ] = FALSE
                nm_pnt  =  pnt_mut[ r, 'Gene_name'][ l ]
                order_pnt  =  pnt_mut[ r, 'mut_order'][ l ]
            } else{
                nm_pnt  =  NULL
                order_pnt  =  NULL
            }

            ### CNA mutations:
            id_cna   =  as.numeric( unlist( str_split( genes_dysfunction$CNA_ID[ i ] , pattern = ',') ) )
            if ( id_cna[1] != 0 ){
                r_cna    =  as.vector( sapply( X = id_cna, FUN = function( x ) which( cna_mut$CNA_ID ==   x  ) ) )
                l_cna    =  cna_mut[ r_cna , 'MalfunctionedByCNA' ]
                l_cna[ is.na( l_cna ) ] = FALSE
                nm_cna   =  cna_mut[ r_cna, 'Gene_names' ][ l_cna ]
                order_cna  =  cna_mut[ r_cna, 'mut_order'][ l_cna ]
            } else {
                nm_cna    =  NULL
                order_cna =  NULL
            }

            # Combine together:
            nm  =  c( nm_pnt, nm_cna )
            rdr = c( order_pnt, order_cna )
            genes_dysfunction[ i, 'order' ]  =  paste( nm[ order( rdr ) ], collapse = ' -> ')
        }
    } else { return( NULL ) }

    write.table( genes_dysfunction, file = file_name,
               sep = '\t', row.names = FALSE, col.names = TRUE, append = FALSE )
    cat( paste0('Order of genes dysfunction saved to the file ', file_name ,' \n') )

    return( genes_dysfunction )
}



#' Function to get data about Variant allele frequencies (VAF)
#'
#' @param file_name Name of file to save data
#' @param pnt_mut data.frame with point mutation info
#' @param data_last data.frame with data of simulation at the last time step
#'
#' @return data.frame with info about Variant allele frequencies
#' @export
#'
#' @examples
#' pnt_mut  =  tugHall_dataset$pnt_mut
#' data_last  =  tugHall_dataset$data_last
#' vf = get_VAF( pnt_mut, data_last, file_name = 'Output/VAF_data.txt')
get_VAF  <-  function( pnt_mut, data_last, file_name = 'Output/VAF_data.txt'){

    pnt_mut_A  =  pnt_mut[ which(  is.na( pnt_mut$MalfunctionedByPointMut ) ) , ]
    pnt_mut_B  =  pnt_mut[ which( !is.na( pnt_mut$MalfunctionedByPointMut ) ) , ]

    if ( nrow(data_last) == 0 ) return( NULL )

    ids  =  str_split( data_last$PointMut_ID, pattern = ',' )
    ids  =  lapply( X = 1:nrow( data_last), FUN = function(x) as.numeric( ids[[ x ]] ) )

    nqu  =  sort( unique( unlist( ids ) ) )  # unique IDs of point mutations
    if ( nqu[1] == 0 ) nqu = nqu[ -1 ]       # exclude intact normal cells

    # start VAF
    VAF  =  NULL

    if ( length( nqu ) == 0 ) return( NULL )

    for( j in 1:length( nqu ) ){
        # wc - which clones have an ID of point mutation
        wc  =  sapply( X = 1:length( ids ), FUN = function( x ) is.element( nqu[j] , ids[[ x ]] ) )

        VAF_1  =  pnt_mut_B[ which( pnt_mut_B$PointMut_ID == nqu[ j ] ) , ]

        ### number of cells: speckled normal cells,
        ###     primary tumor cells and metastatic cells:
        if ( any( data_last[ which( wc ), 'type' ] == 'normal' ) ){
            sm  =  sum( as.numeric( data_last[ which( wc & data_last$type ==  'normal' ),  'N_cells'  ] ) )
            if ( length( sm ) == 0 ) sm = 0
            VAF_1$N_speckled_normal =  sm
        } else VAF_1$N_speckled_normal =  0

        if ( any( data_last[ which( wc ), 'type' ] == 'primary' ) ){
            VAF_1$N_primary  =  sum( as.numeric( data_last[ which( wc & data_last$type ==  'primary' ),  'N_cells'  ] ) )
        } else VAF_1$N_primary  =  0

        if ( any( data_last[ which( wc ), 'type' ] == 'metastatic' ) ){
            VAF_1$N_metastatic =  sum( as.numeric( data_last[ which( wc & data_last$type ==  'metastatic' ),  'N_cells'  ] ) )
        } else VAF_1$N_metastatic =  0

        # add copy number of original allele A:
        VAF_A  =  pnt_mut_A[ which( pnt_mut_A$PointMut_ID == nqu[ j ] ) , ]

        VAF_1$Copy_number_A  =  VAF_A$Copy_number

        VAF  =  rbind( VAF, VAF_1)
    }

    # Total number of speckled normal cells (no driver mutations but at least one passenger )
    if ( any( data_last[ , 'type' ] == 'normal' & ( data_last[ , 'CNA_ID'] != '0' | data_last[ , 'PointMut_ID'] != '0' ) ) ) {
      w = which( data_last$type == 'normal' & ( data_last[ , 'CNA_ID'] != '0' | data_last[ , 'PointMut_ID'] != '0' ) )
      VAF$N_speckled_normal_total = sum( as.numeric( data_last[ w, 'N_cells' ] ) )
    } else  VAF$N_speckled_normal_total  =  0

    # Total number of primary tumor cells (at least one driver)
    if ( any( data_last[ , 'type' ] == 'primary' ) ) {
        VAF$N_primary_total = sum( as.numeric( data_last[ which( data_last$type == 'primary' ), 'N_cells' ] ) )
    } else  VAF$N_primary_total  =  0

    # Total number of metastatic cells
    if ( any( data_last[ , 'type' ] == 'metastatic' ) ) {
        VAF$N_metastatic_total = sum( as.numeric( data_last[ which( data_last$type == 'metastatic' ), 'N_cells' ] ) )
    } else  VAF$N_metastatic_total  =  0

    ### Save VAF to the file:
    write.table( VAF, file = file_name, append = FALSE, row.names = FALSE, sep="\t" )
    print( paste0("VAF data for allele B and A is saved to the file ", file_name , ' \n ' ) )

    return( VAF )
}


#' Function to get Variant allele frequencies (VAF) based on rho input parameters
#'
#' @param vf data.frame getting from get_VAF() function
#' @param rho Vector of rho parameter in the range (0,1)
#' @param file_name Name of file to save VAF
#'
#' @return VAF for different rho with separation for metastatic cells and (primary tumor + speckled normal) cells
#' @export
#'
#' @examples
#' pnt_mut  =  tugHall_dataset$pnt_mut
#' data_last  =  tugHall_dataset$data_last
#' if ( !dir.exists('./Output') ) dir.create('./Output')
#' vf = get_VAF( pnt_mut, data_last, file_name = 'Output/VAF_data.txt')
#' VAF = get_rho_VAF( vf = vf, rho = c( 0.0, 0.1, 0.5 ) , file_name = './Output/VAF.txt' )
get_rho_VAF  <-  function( vf = NULL, rho = c( 0.0, 0.1, 0.5 ) , file_name = './Output/VAF.txt' ){
    # rho is an admixture rate of intact normal cells, it can be vector of numbers
    # vf is a VAF data getting from get_VAF function
    # file_name is  file name for VAF file

    if ( min(rho) < 0 | max(rho) >1 ) return( NULL )
    nq_i  =  unique( vf$Ref_pos )
    if ( length(nq_i) < 1 ) return( NULL )

    N_speckled_normal_total  =  vf$N_speckled_normal_total[1]
    N_primary_total          =  vf$N_primary_total[1]
    N_metastatic_total       =  vf$N_metastatic_total[1]

    VAF  =  NULL
    for( k in 1:length( rho ) ){
        # Scale for admixture rate of intact normal cells rho[ k ]:
        k_scale  =  ( 1 - rho[ k ] )   # rho[ k ] * ( N_primary_total + N_speckled_normal_total ) / N_primary_total

        for( i in nq_i ){
            w  =  which( vf$Ref_pos == i )
                                        # for primary tumor and speckled normal cells
            numenator_N    =  sum( vf[ w , c( 'N_speckled_normal', 'N_primary' ) ] * vf[ w, 'Copy_number'] ) / ( N_primary_total + N_speckled_normal_total )
            denominator_N  =  sum( vf[ w , c( 'N_speckled_normal', 'N_primary' ) ] * ( vf[ w, 'Copy_number'] + vf[ w, 'Copy_number_A'] ) ) / ( N_primary_total + N_speckled_normal_total )
                                        # for metastatic cells
            numenator_M    =  sum( vf[ w , 'N_metastatic'] * vf[ w, 'Copy_number'] ) / N_metastatic_total
            denominator_M  =  sum( vf[ w , 'N_metastatic'] * ( vf[ w, 'Copy_number'] + vf[ w, 'Copy_number_A'] ) ) / N_metastatic_total
                                        # VAF calculations:
            VAF_N_rho  =  k_scale * numenator_N / ( 2*( 1 - k_scale ) + k_scale * denominator_N )
            VAF_M_rho  =  k_scale * numenator_M / ( 2*( 1 - k_scale ) + k_scale * denominator_M )
                                        # save to data.frame:
            VAF_1 = data.frame(site = i,
                               Chr = vf[ w[1], 'Chr' ] ,
                               gene = vf[ w[1], 'Gene_name' ] ,
                               rho = rho[ k ],
                               VAF_primary = VAF_N_rho,
                               VAF_primary_numerator   =  k_scale * numenator_N,
                               VAF_primary_denominator    =  ( 2*( 1 - k_scale ) + k_scale * denominator_N ),
                               VAF_metastatic = VAF_M_rho,
                               VAF_metastatic_numerator = k_scale * numenator_M,
                               VAF_metastatic_denominator =  ( 2*( 1 - k_scale ) + k_scale * denominator_M )
                    )
            VAF_1[ is.na.data.frame( VAF_1 ) ]  =  0  #  division by 0 if rho = 1
            VAF  =  rbind( VAF, VAF_1 )
        }
    }

    Stop_reason = if ( sum(pck.env$data_last$N_cells) > pck.env$censor_cells_number ){
        Stop_reason  =  'Cells number'
    } else {
        if ( pck.env$data_last$Time[1] >= pck.env$censor_time_step ){
            Stop_reason  =  'Time step'
        } else {
            Stop_reason  =  'Real time'
        }
    }

    VAF$Stop_reason  =  Stop_reason

    write.table( VAF, file = file_name, append = FALSE, sep = '\t',
                    row.names = FALSE, col.names = TRUE )
    cat( paste0( ' VAF is saved in the file ', file_name ) )
    return( VAF )
}

