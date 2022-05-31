### Define the FOLDERS and files' names ---------------------------------------------------
## Create folders:  /Input, /Output and /Figures

#' Function to define all the files names
#'
#' @param mainDir Working directory
#' @param sbdr_Input Sub directory for input files
#' @param sbdr_Output Sub directory for output files
#'
#' @return NULL, but all file names are defined in GLOBAL environment
#' @export
#'
#' @examples
#' define_files_names()
define_files_names  <-  function( mainDir = getwd(), sbdr_Input = '/Input', sbdr_Output = '/Output' ){

    if (! file.exists( paste0( mainDir, sbdr_Output ) ) ){  dir.create( file.path( mainDir, sbdr_Output ) ) }

    if (! file.exists( paste0( mainDir, sbdr_Input ) ) ){  dir.create( file.path( mainDir, sbdr_Input ) ) }

    if (! file.exists( paste0( mainDir, '/Figures' ) ) ){  dir.create( file.path( mainDir, 'Figures' ) ) }

    ### Files to output and input data
    genefile       <<-   paste0( mainDir, sbdr_Input, '/gene_hallmarks.txt' )    # gene file
    clonefile      <<-   paste0( mainDir, sbdr_Input, '/cloneinit.txt'      )    # initial Cells

    ### Output files
    geneoutfile    <<-   paste0( mainDir, sbdr_Output, '/geneout.txt'       )    # Gene Out file with Hallmarks
    cloneoutfile   <<-   paste0( mainDir, sbdr_Output, '/cloneout.txt'      )    # output information of simulation
    logoutfile     <<-   paste0( mainDir, sbdr_Output, '/log.txt'     )    # log file to save the input information of simulation - "log.txt"
    ### Output/Weights.txt               # file with gene weights for hallmarks
    file_monitor   <<-   './Sim_monitoring.txt'
}
### Define the gene map - chromosomal locations --------------------------



#' Define genes' location in chromosome
#'
#' @param file_input is a name of file to input where the information about genes location is defined. That is loaded from CCDS database https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/
#' @param genes_list is a list of genes' names like CCDS4107.1 in the CCDS database.
#'
#' @return Function returns the table of genes' locations in DNA
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' define_gene_location()
#' file_input  =  'Input/CCDS.current.txt'
#' genes_list  =  c( 'CCDS4107.1', 'CCDS8702.1', 'CCDS43171.1', 'CCDS11118.1' )
#' define_gene_location( file_input = file_input,  genes_list = genes_list )
define_gene_location  <-  function( file_input  =  'Input/CCDS.current.txt',
                                    genes_list  =  c( 'CCDS4107.1', 'CCDS8702.1',
                                                      'CCDS43171.1', 'CCDS11118.1' ) ){

    ### Make a map of genes with sorting of start position for each chromosome:
    gene_map  <<-   make_map(f_out    =  'Input/gene_map.txt',
                             ls   =  genes_list,
                             f_in =  file_input )
    gene_map  <<-  order_gene_map( gene_map )  ### We have to be sure in the sorting of start position for each chromosome

    write.table(gene_map, file = 'Output/gene_MAP.txt', col.names = TRUE,
                sep = "\t", row.names = FALSE)
}

### Define the PARAMETERS ------------------------------------------------

#' Define all the parameters for a simulation
#'
#' @param   E0 Parameter in the division probability, numeric type only
#' @param   F0 Parameter in the division probability, numeric type only
#' @param   m0 Mutation probability for point mutation, numeric type only
#' @param   uo Oncogene mutation probability, numeric type only
#' @param   us Suppressor mutation probability, numeric type only
#' @param   s0 Parameter in the sigmoid function, numeric type only
#' @param   k0 Environmental death probability, numeric type only
#' @param   d0 Initial probability to divide cells, numeric type only
#' @param   censore_n Max cell number where the program forcibly stops, integer type only
#' @param   censore_t Max time where the program forcibly stops, integer type only
#' @param   time_stop Max time in seconds of running after that the program forcibly stops, integer type only
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
#' @param CF Compaction factor, logical type only. True means 'to use', False means 'do not use' Compaction factor for hallmarks variables
#' @param model Name of the model to use. Can be  'proportional_metastatic' or 'threshold_metastatic' or 'simplified'
#' @param read_fl Indicator to read file or not, logical type only
#' @param file_name File name to rad all the parameters, it is used only if read_fl == TRUE
#'
#' @return Values of all the parameters
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' define_parameters( read_fl = TRUE , file_name = './Input/parameters.txt' )
#' define_parameters( read_fl = FALSE )
define_parameters  <-  function( E0 =  1E-4, F0 =  10, m0 =  1E-7, uo =  0.9, us =  0.9,
                                 s0 =  10, k0 =  0.12, d0 =  0.4, censore_n = 10^5,
                                 censore_t = 50, m_dup  = 1E-8, m_del  = 1E-8,
                                 lambda_dup  = 5000, lambda_del  = 7000,
                                 uo_dup  = 0.8, us_dup  = 0.5, uo_del  = 0, us_del  = 0.8,
                                 CF  =  TRUE,
                                 model  =  c( 'proportional_metastatic', 'threshold_metastatic', 'simplified' )[ 1 ],
                                 time_stop = 120,
                                 read_fl = FALSE, file_name ='./Input/parameters.txt',
                                 n_repeat = 1000, monitor  =  TRUE ){
    if ( read_fl ){
        data_log  =  read.table( file = file_name, sep = '\t', stringsAsFactors = FALSE )
        names( data_log )  =  c( 'var', 'value' )
        # Model definition
        Compaction_factor  <<-  as.logical( data_log[ which( data_log$var == 'Compaction_factor' ), 2 ] )
        model_name         <<-  data_log[ which( data_log$var == 'model_name' ), 2 ]
        time_stop          <<-  as.numeric( data_log[ which( data_log$var == 'time_stop' ), 2 ] )  # max time in seconds
        n_repeat           <<-  as.numeric( data_log[ which( data_log$var == 'n_repeat' ), 2 ] )  # max number of repetitions
        # Parameters:
        E0 <<-  as.numeric( data_log[ which( data_log$var == 'E' ), 2 ] )       # parameter in the division probability
        F0 <<-  as.numeric( data_log[ which( data_log$var == 'F' ), 2 ] )       # parameter in the division probability
        m0 <<-  as.numeric( data_log[ which( data_log$var == 'm0' ), 2 ] )     # mutation probability
        uo <<-  as.numeric( data_log[ which( data_log$var == 'uo' ), 2 ] )        # oncogene mutation probability
        us <<-  as.numeric( data_log[ which( data_log$var == 'us' ), 2 ] )        # suppressor mutation probability
        s0 <<-  as.numeric( data_log[ which( data_log$var == 's' ), 2 ] )         # parameter in the sigmoid function
        d0 <<-  as.numeric( data_log[ which( data_log$var == 'd0' ), 2 ] )      # Initial probability to divide cells
        k0 <- as.character( data_log[ which( data_log$var == 'k0' ), 2 ] )   # Environmental death probability
        if ( is.na( k0 ) ) {
            k0 <<- 1 - (1 + d0) ^ (-1)
        }
        else {
            k0 <<-  as.numeric( k0 )
        }
        ### Additional parameters of simulation
        censore_n <<- as.numeric( data_log[ which( data_log$var == 'censore_n' ), 2 ] )       # Max cell number where the program forcibly stops
        censore_t <<- as.numeric( data_log[ which( data_log$var == 'censore_t' ), 2 ] )       # Max time where the program forcibly stops
        ### New parameters for CNA:
        m_dup  <<- as.numeric( data_log[ which( data_log$var == 'm_dup' ), 2 ] ) # mutation probability for duplication
        m_del  <<- as.numeric( data_log[ which( data_log$var == 'm_del' ), 2 ] ) # mutation probability for deletion
        lambda_dup  <<- as.numeric( data_log[ which( data_log$var == 'lambda_dup' ), 2 ] )  # CNA duplication average length (of the geometrical distribution for the length)
        lambda_del  <<- as.numeric( data_log[ which( data_log$var == 'lambda_del' ), 2 ] )  # CNA deletion average length
        uo_dup  <<- as.numeric( data_log[ which( data_log$var == 'uo_dup' ), 2 ] ) # Gene malfunction probability by CNA duplication for oncogene
        us_dup  <<- as.numeric( data_log[ which( data_log$var == 'us_dup' ), 2 ] )   # Gene malfunction probability by CNA duplication for suppressor
        uo_del  <<- as.numeric( data_log[ which( data_log$var == 'uo_del' ), 2 ] )   # Gene malfunction probability by CNA deletion    for oncogene
        us_del  <<- as.numeric( data_log[ which( data_log$var == 'us_del' ), 2 ] ) # Gene malfunction probability by CNA deletion    for suppressor
        monitor <<- as.logical( data_log[ which( data_log$var == 'monitor' ), 2 ] )
    } else {

        # Model definition:
        Compaction_factor  <<-  CF
        model_name         <<-  model
        # Parameters:
        E0 <<-  E0       # parameter in the division probability
        F0 <<-  F0         # parameter in the division probability
        m0 <<-  m0      # mutation probability
        uo <<-  uo        # oncogene mutation probability
        us <<-  us        # suppressor mutation probability
        s0 <<-  s0         # parameter in the sigmoid function
        k0 <<-  k0        # Environmental death probability
        d0 <<-  d0       # Initial probability to divide cells
        ### Additional parameters of simulation
        censore_n <<- censore_n       # Max cell number where the program forcibly stops
        censore_t <<- censore_t         # Max time where the program forcibly stops
        time_stop  <<-  time_stop     # Max time in seconds of running after that the program forcibly stops
        n_repeat     <<-   n_repeat     # Max number of repetition of the program until the NON-ZERO output will be getting
        ### New parameters for CNA:
        m_dup  <<- m_dup # mutation probability for duplication
        m_del  <<- m_del # mutation probability for deletion
        lambda_dup  <<- lambda_dup  # CNA duplication average length (of the geometrical distribution for the length)
        lambda_del  <<- lambda_del  # CNA deletion average length
        uo_dup  <<- uo_dup # Gene malfunction probability by CNA duplication for oncogene
        us_dup  <<- us_dup   # Gene malfunction probability by CNA duplication for suppressor
        uo_del  <<- uo_del   # Gene malfunction probability by CNA deletion    for oncogene
        us_del  <<- us_del # Gene malfunction probability by CNA deletion    for suppressor
        monitor <<- monitor  # The indicator to make monitor file during a simulation or do not make
    }
}

#' Function to print GLOBAL parameters
#'
#' @return Message with values of all the GLOBAL parameters
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' define_parameters( read_fl = FALSE )
#' define_compaction_factor()
#' print_parameters()
print_parameters  <-  function(){

    msg  =  c(
        'Model definition:  \n ' ,
        'Compaction_factor = ', Compaction_factor, '\n',
        'model_name  =  ', model_name, '\n',
        'Parameters:  \n',
        'parameter of the division probability E0 =  ', E0, '\n',
        'another parameter of the division probability F0  = ',  F0, '\n',
        'mutation probability m0 =  ', m0, '\n',
        'oncogene mutation probability uo = ', uo, '\n',
        'suppressor mutation probability  us  =  ', us, '\n',
        'parameter in the sigmoid function  s0  =  ', s0, '\n',
        'Environmental death probability  k0 =  ',  k0, '\n',
        'Initial probability to divide cells  d0  =  ',  d0, '\n',
        'Additional parameters of simulation  \n ',
        'Max cell number where the program forcibly stops  censore_n  = ',  censore_n,  '\n',
        'Max time steps where the program forcibly stops  censore_t  = ',  censore_t,  '\n',
        'Max time (in seconds) where the program forcibly stops time_stop  =  ',  time_stop,  '\n',
        'Max number of repetition of the program until the NON-ZERO output will be getting, n_repeat  =  ', n_repeat ,   '\n',
        'New parameters for CNA:  \n',
        'mutation probability for duplication  m_dup  =  ', m_dup ,  '\n',
        'mutation probability for deletion',  m_del,  '\n',
        'CNA duplication average length (of the geometrical distribution for the length)  lambda_dup  =  ', lambda_dup ,  '\n',
        'CNA deletion average length  lambda_del  = ', lambda_del ,  '\n',
        'Gene malfunction probability by CNA duplication for oncogene  uo_dup  =  ', uo_dup ,  '\n',
        'Gene malfunction probability by CNA duplication for suppressor  us_dup  =  ', us_dup ,  '\n',
        'Gene malfunction probability by CNA deletion for oncogene  uo_del  = ', uo_del ,  '\n',
        'Gene malfunction probability by CNA deletion for suppressor  us_del  = ', us_del,  '\n',
        'Compaction factor is applied if variable Compaction_factor ==  TRUE \n',
        'Compaction factor for apoptosis hallmark CF$Ha = ', CF$Ha, ' \n',
        'Compaction factor for angiogenesis hallmark CF$Hb = ', CF$Hb, ' \n',
        'Compaction factor for growth/antigrowth hallmark CF$Hd = ', CF$Hd, ' \n',
        'Compaction factor for immortalization hallmark CF$Hi = ', CF$Hi, ' \n',
        'Compaction factor for invasion/metastasis hallmark CF$Him = ', CF$Him,
        '\n Monitoring: \n indicator monitor  =  ', monitor, '\n \n '
    )

    cat( paste0( msg, collapse = ' ' ) )

    if ( model_name  ==  'simplified' ) {
        msg  =  c(
            'The model is simplified that means next: \n ',
            '    1) all the hallmarks are defined but do not affect \n ',
            '       excepting hallmark of growth/antigrowth Hd \n',
            '    2) apoptosis trial is deleted \n',
            '    3) all the cells are metastatic \n',
            '    4) Hayflic limitation is deleted \n',
            '    5) Only exponential growth is simulated'
        )
        message( paste0( msg, collapse = ' ' ) )
    }
}

#' Define compaction factor
#'
#' @param cf Data frame with compaction factors for all the hallmarks, for example, data.frame( Ha = 1, Hb = 1, Hd = 1, Hi = 1, Him = 1 )
#' @param read_fl Indicator to read file or not, logical type only
#' @param file_name File name to rad all the parameters, it is used only if read_fl == TRUE
#'
#' @return Data frame with with compaction factors for all the hallmarks
#' @export
#'
#' @examples
#' copy_files_to_Input()
#' define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' )
#' CF1 = CF
#' cf = data.frame( Ha = 0.1, Hb = 0.2, Hd = 0.7, Hi = 1, Him = 0.5 )
#' define_compaction_factor( cf = cf, read_fl = FALSE )  # View( c( CF, CF1 ) ) to compare
define_compaction_factor  <-  function( cf = data.frame( Ha = 1, Hb = 1, Hd = 1,
                                                         Hi = 1, Him = 1 ),
                                        read_fl = TRUE , file_name = './Input/CF.txt' ){
    if ( read_fl ){
        data_log  =  read.table( file = file_name, sep = '\t', stringsAsFactors = FALSE )
        names( data_log )  =  c( 'var', 'value' )

        cf$Ha   =  data_log$value[ data_log$var == 'apoptosis' ]
        cf$Hb   =  data_log$value[ data_log$var == 'angiogenesis' ]
        cf$Hd   =  data_log$value[ data_log$var == 'growth' ]
        cf$Hi   =  data_log$value[ data_log$var == 'immortalization' ]
        cf$Him  =  data_log$value[ data_log$var == 'invasion' ]
    }

    CF  <<-  cf
}

