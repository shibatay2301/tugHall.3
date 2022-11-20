### Pipeline to perform Approximate Bayesian Computation for tugHall dataset:

### Here is just an example, so, please, be careful to use this pipeline for your data.


# 1. GET DATASET -------------------------------------------------------------

library( 'utils' )

setwd( '../ABC' )
if ( !dir.exists( './Data') ) dir.create( './Data' )

### GET dataset from repository https://data.mendeley.com:
addrss  =  'https://data.mendeley.com/public-files/datasets/spszxd8r3z/files/ae40c72a-4b93-40c1-a74b-77fdc59651d6/file_downloaded'
setwd('./Data/')
download.file  (url = addrss, destfile = './Archive.zip' )
unzip(zipfile = './Archive.zip')
setwd('../')
unlink('./Data/Archive.zip')

### To understand the data, be kind read Description.pdf file:
###     https://data.mendeley.com/datasets/spszxd8r3z/1

# Here we use
# 'Initial_parameters_Discrete_ALL.txt' as dataset of parameters:
data.param  =  read.csv2(file = "./Data/Initial_parameters_Discrete_ALL.txt",
                         header = TRUE, sep = "\t")

### For summary statistics of simulations:
data.statistics  =  read.csv2(file = "./Data/data_base_MODELS_PRIMARY.txt",
                         header = TRUE, sep = "\t")

### That is data of VAF for primary tumor cells only (metastatic cells are removed).

### NEXT, we use:\
##  - name_model    =  'STRONG' with condition ( im' = 1 ) to be metastatic
##  - name_weights  =  'Discrete' with discrete set of weight values
##  - name_init     =  'Mutated_cell_in_Thousand_cells'

### Because these condition yield with maximal number of non-zero output 34059 of 34602 rows:
cndt  =  which( data.statistics$name_weights   ==  'Discrete' &
                    data.statistics$name_model == 'STRONG'    &
                    data.statistics$name_init  == 'Mutated_cell_in_Thousand_cells' )

data.statistics  =  data.statistics[ cndt, ]

cndt  =  data.param$ID_Simulation  %in% data.statistics$ID_Simulation
data.param   =   data.param[ cndt, ]
rm ( cndt, addrss )


# 2.  GET OBSERVATION POINT -----------------------------------------------

### Here we use OUR prepared file 'VAF_observations.csv'
###  TARGET DATA OR OBSERVATION DATA

############## Function to get the observation data for different patients
read_obs_data  =  function(filename = "VAF_observations.csv") {
    W <- read.csv2(file = filename, header = TRUE)

    W$Tumor_Sample_Barcode <- as.character( W$Tumor_Sample_Barcode )
    ########################## to define how many patients with 1 mutation for each gene

    x <- W[which(W$FILTER == "PASS"),]
    y <- x[which(x$Hugo_Symbol == "APC" | x$Hugo_Symbol == "KRAS" | x$Hugo_Symbol == "TP53" | x$Hugo_Symbol == "PIK3CA" ),]
    z <- table(y$Tumor_Sample_Barcode)
    s <- names(z[which(z < 5)])   # s - id of patients with 1 mutation / gene (probably)

    ### check it and exclude others:

    v_TR <- rep(TRUE,length(s))
    for (u in 1:length(s)) {
        id_patient <- s[u]
        gns <- as.character( y$Hugo_Symbol[ which(y$Tumor_Sample_Barcode == id_patient) ] )
        if ( length(unique(gns)) != length(gns)   ) v_TR[u] <- FALSE   # exclude double mutations
    }

    # names of patients with only 1 mutation for 1 gene or less
    s <- s[v_TR]

    ##############################

    W <- y

    data.target <- as.data.frame(matrix(data = 0, nrow = length(s), ncol = 5 ))
    names(data.target) <- c("APC", "KRAS", "TP53", "PIK3CA", "Name")

    for (ID in 1:length(s)) {

        id_patient <- s[ID]

        W1 <- W[(W$Tumor_Sample_Barcode == id_patient),]

        data.target[ID, 1] <- ifelse(length(W1$t_VAF[which(W1$Hugo_Symbol == "APC")]) == 0, 0, W1$t_VAF[which(W1$Hugo_Symbol == "APC")] )
        data.target[ID, 2] <- ifelse(length(W1$t_VAF[which(W1$Hugo_Symbol == "KRAS")]) == 0, 0, W1$t_VAF[which(W1$Hugo_Symbol == "KRAS")] )
        data.target[ID, 3] <- ifelse(length(W1$t_VAF[which(W1$Hugo_Symbol == "TP53")]) == 0, 0, W1$t_VAF[which(W1$Hugo_Symbol == "TP53")] )
        data.target[ID, 4] <- ifelse(length(W1$t_VAF[which(W1$Hugo_Symbol == "PIK3CA")]) == 0, 0, W1$t_VAF[which(W1$Hugo_Symbol == "PIK3CA")] )
        data.target[ID, 5] <- id_patient
    }
    data.target[is.na(data.target)] <- 0

    return(data.target)
}

data.target_ALL  =  read_obs_data(filename = "VAF_observations.csv")

# Get histogram over all the patients
hist( data.target_ALL$APC, breaks = 18)

### GET ONE observation:
ID = 40         # Use any number to choose patient ID
data.target  =  data.target_ALL[ ID, ]
patient      =  data.target_ALL[ ID, "Name"]
data.target  =  data.target[ , -5 ]
rm ( ID, data.target_ALL )


# 3. GET A Posterior Distribution based on given dataset ------------------------------------------------

library('abc')

# Restrict number of columns for the summary statistics:
data.sum_stat  =  data.statistics[ , c( 'APC_max_1', 'KRAS_max_1',  'TP53_max_1',  'PIK3CA_max_1' ) ]

# Transform characters to numeric formats:
for ( i in 1:ncol(data.sum_stat ) ){
    data.sum_stat[ , i ]  = as.numeric( data.sum_stat[ , i ] )
}

data.param$Mutated_Gene  =  as.integer( as.factor( data.param$Mutated_Gene ) )
for ( i in 1:ncol(data.param ) ){
    data.param[ , i ]  = as.numeric( data.param[ , i ] )
}

### ABC
### Choose tolerance:
tolerance  =  0.001
rejection  =  abc( target = data.target,
                   param  = data.param,
                   sumstat = data.sum_stat,
                   tol = tolerance,
                   method = 'rejection')


hist( x = rejection, unadj = FALSE, true = NULL, file = NULL,
          postscript = FALSE, onefile = TRUE,
          col.hist = "grey", col.true = "red",
          caption = NULL, breaks = 12 )


# 4.GET ACCURATE MAP based on SAMPLING ------------------------------------

library( 'EasyABC' )



toy_model  =  function(x){ c( 100 * exp( - (x[1] - 30) ** 2 / 32 ),
                              100 * exp( - (x[2] - 55) ** 2 / 32 ) ) }


toy_prior  =  list( c( "unif", 0, 100 ), c( "unif", 0, 100 ) )
sum_stat_obs  =  c( 100, 100 )
set.seed(1)


############# REJECTION ABC
n=1000

ABC_rej  =  ABC_rejection( model = toy_model, prior = toy_prior,
                           nb_simul = n,
                           summary_stat_target = sum_stat_obs,
                           tol = 0.008,
                           progress_bar = TRUE )
ABC_rej$param
ABC_rej$stats

hist( ABC_rej$param[ , 1] )
hist( ABC_rej$param[ , 2] )


############# Adaptive ABC

### Ref:
# Beaumont, M. A., Cornuet, J., Marin, J., and Robert, C. P. (2009)
# Adaptive approximate Bayesian computation. Biometrika, 96, 983–990.

tolerance  =  c( 1E-1, 1E-2 )
n = 20

ABC_Beaumont  =  ABC_sequential( method = "Beaumont",
                                 model  = toy_model,
                                 prior  = toy_prior,
                                 nb_simul = n,
                                 summary_stat_target = sum_stat_obs,
                                 tolerance_tab = tolerance,
                                 verbose = TRUE )


ABC_Beaumont$weights
ABC_Beaumont$param

hist( ABC_Beaumont$param[, 1])
hist( ABC_Beaumont$param[, 2])


