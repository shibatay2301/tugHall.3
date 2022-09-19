###   I)  Define CLONE'S CLASSES ----------------------------------------------------------

#' Class 'Clone' for clones
#'
#' @field id numeric. ID of a clone
#' @field parent numeric. Parent ID (for first - 0)
#' @field N_cells numeric. Number of cells in clone
#' @field c numeric. Split counter as average value for all cells in clone
#' @field d numeric. Probability of division
#' @field i numeric. Probability of Hayflick limit
#' @field m numeric. Probability that gene normal function is destroyed due to epigenome abnormality / mutation rate
#' @field a numeric. Probability of apoptosis for a cell in the clone
#' @field s numeric. Coefficient in the sigmoid function of the mutation density
#' @field k numeric. Probability of cell death by environment
#' @field E numeric. Coefficient of friction term against to the split probability.
#' @field Nmax numeric. Coefficient for determination the max number of cells \cr
#' that can exist in the primary tumor (Nmax = 1/E) \cr
#' @field im numeric. Probability of the invasion/ metastatic transformation
#' @field Ha numeric. Apoptosis hallmark value
#' @field Him numeric. Invasion/ metastasis hallmark
#' @field Hi numeric. Mitotic restriction hallmark (immortalization hallmark)
#' @field Hd numeric. Growth/antigrowth hallmark (division rate hallmark)
#' @field Hb numeric. Angiogenesis hallmark
#' @field gene numeric. Vector of flags for each genes if they have driver mutation
#' @field pasgene numeric. Vector of flags for each genes if they have passenger mutation
#' @field PointMut_ID numeric. ID of point mutation in list of objects of class 'Point_Mutations'
#' @field CNA_ID numeric. ID of CNA mutation in list of objects of class 'CNA_Mutations'
#' @field mutden numeric. Gene mutation density
#' @field invasion logical. Indicator that clone is metastatic (invasion/metastatic transformation occured or not)
#' @field primary logical. Logical variable is clone primary tumor or not (normal)
#' @field birthday numeric. Time step of birth of clone
#'
#' @return NULL
#'
#' @export
#' @examples
#' clone = tugHall_dataset$clones[[ 1 ]]
#' print(clone$Ha)
#' print(clone$N_cells)
#' clone$calcApoptosis()    # to calculate apoptosis death probability based on mutation density
clone <- setRefClass(
    # name of the class
    Class = "Clone",
    # Field
    fields = list(
        id = "numeric",          # identificator
        parent = "numeric",      # parent ID (for first - 0)
        N_cells = "numeric",     # Number of cells in clone
        c = "numeric",           # Split counter
        d = "numeric",           # Probability of division
        i = "numeric",           # Probability of Hayflick limit
        m = "numeric",           # Probability that gene normal function is destroyed due to epigenome abnormality
        a = "numeric",           # Probability of apoptosis
        s = "numeric",           # Coefficient in the of apoptosis
        k = "numeric",           # Probability of cell death by environment
        E = "numeric",           # Coefficient of friction term against to the split probability,
        # Coefficient for determination the max number of cells that can exist in the primary tumor (Nmax = 1/E)
        Nmax = "numeric",        # The max number of cells that can exist in the primary tumor (Nmax = 1/E)
        im = "numeric",          # Invasion/ metastasis probability
        Ha = "numeric",          # Apoptosis probability difference (apoptosis gene x weight)
        Him = "numeric",         # Invasion/ metastasis probability difference (invasion/ metastasis gene x weight)
        Hi = "numeric",          # Mitotic restriction probability difference (immortalized gene x weight)
        Hd = "numeric",          # Divide probability difference (cancer gene x weight)
        Hb = "numeric",          # Friction coefficient (angiogenic gene x weight)
        gene = "numeric",        # Flag for cancer gene function deficit (for drivers)
        pasgene = "numeric",     # Flag for cancer gene as passenger dysfunction
        # deleted: posdriver = "character", # position of cancer gene damage (function deficit)
        # deleted: pospasngr = "character", # position of cancer gene damage (maintenance of function)
        PointMut_ID  =  "numeric", # ID of point mutation in object PointMut / pnt
        CNA_ID       =  "numeric", # ID of CNA mutation in object CNA
        mutden = "numeric",      # Gene mutation density
        invasion = "logical",    # Wetting/Displacement flag:    TRUE: Wetting/Displacement      FALSE: Limited pattern
        primary  = "logical",    # Logical variable is clone primary tumor or not
        birthday = "numeric"    # Time step of birth of cell
        # lenCDS      = "numeric"     # length of all oncogenes of interest
        # lenRNA      = "numeric"      # length of all oncogenes of interest including introns
    ),

    # Method
    methods = list(
        # Initialize
        initialize = function(gene_size, id=1, parent=0, c=0, d=d0, i=1, m=m0,  N_cells = 1,
                              mutden=0, a=0, k=k0, E=E0, Nmax=0, gene=NULL, pasgene=NULL,
                              PointMut_ID = 0, CNA_ID = 0,
                              # deleted: posdriver=NULL, pospasngr=NULL,
                              invasion=FALSE, primary=FALSE, s=s0, birthday=0) {
            id <<- id
            parent <<- parent
            N_cells <<- N_cells
            c <<- c
            d <<- d
            i <<- i
            m <<- m
            s <<- s
            birthday <<- birthday
            mutden <<- mutden
            if (is.null(a)) {
                a <<- 1/(1+exp(-s*(mutden - 0.5)))
            } else {
                a <<- a
            }
            k <<- k
            E <<- E
            Nmax <<- 1.0 / E
            im <<- 0
            Ha <<- 0
            Him <<- 0
            Hi <<- 0
            Hd <<- 0
            Hb <<- 0

            if (is.null(gene)) {
                gene <<- rep(0, gene_size)
            } else {
                gene <<- gene
            }
            if (is.null(pasgene)) {
                pasgene <<- rep(0, gene_size)
            } else {
                pasgene <<- pasgene
            }
            # if (is.null(posdriver)) {
            #    posdriver <<- rep("", gene_size)
            # } else {
            #    posdriver <<- posdriver
            # }
            # if (is.null(pospasngr)) {
            #     pospasngr <<- rep("", gene_size)
            # } else {
            #     pospasngr <<- pospasngr
            # }
            PointMut_ID <<-  PointMut_ID
            CNA_ID     <<-  CNA_ID
            invasion   <<-  invasion
            primary    <<-  primary
            # lenCDS     <<-  sum( onco$cds )
            # lenRNA     <<-  integer(0)
        },
        # Apoptosis
        calcApoptosis = function() {
            #            if (mutden <= sum(gene)/length(gene)) {
            #                a1 = 1/(1+exp(s*(mutden - 0.5)))
            #                mutden <<- sum(gene)/length(gene)
            #                a2 = 1/(1+exp(s*(mutden - 0.5)))
            #                a <<- a - (a1 - a2)

            mutden <<- sum(gene) / length(gene)
            a <<- 1 / ( 1 + exp( -1 * s * ( mutden - 0.5 ) ) )
            if ( a < 0 ) {
                a <<- 0
            }

        },
        # Aggregate
        calcMutden = function() {
            mutden <<- sum( gene ) / length( gene )
        }
    )
)


#' Class 'Environ'
#'
#' @field T numeric. Time counter
#' @field N numeric. Number of normal cells
#' @field P numeric. Number of primary tumor cells
#' @field M numeric. Number of metastatic cells
#' @field F numeric. Coefficient that determines the maximal number of cells in pool of primary tumor cells
#' @field c numeric. Average number of divisions in pool of clones
#' @field d numeric. Mean value of splitting probability
#' @field i numeric. Average value of immortalization probability
#' @field a numeric. Average value of apoptosis probability
#' @field k numeric. Average probability of cell death via environment death
#' @field E numeric. Average value of coefficients of friction term
#' @field Nmax numeric. Maximal number of primary tumor cells that can exist in pool of clones
#' @field im numeric. Average value of invasion/metastasis probability
#' @field Ha numeric. Average value of apoptosis hallmark Ha
#' @field Him numeric. Average value of invasion/metastasis hallmark Him
#' @field Hi numeric. Average value of immortalization hallmark Hi
#' @field Hd numeric. Average value of growth/antigrowth hallmark Hd
#' @field Hb numeric. Average value of angiogenesis hallmark Hb
#' @field type numeric. Invasion / metastatic ratio
#' @field gene numeric. Cancer gene damage rate
#' @field mutden numeric. Average density of gene malfunction
#' @field last_id numeric. Maximal ID in the pool of clones.
#'
#' @return NULL
#' @export
#'
#' @examples
#' env = tugHall_dataset$env
#' print( env )
#' env$initFields()
environ <- setRefClass(
    # the class name
    Class = "Environ",

    # Fields
    fields = list(
        T = "numeric",           # time counter
        N = "numeric",           # number of normal cells
        P = "numeric",           # number of primary tumor cells
        M = "numeric",           # number of infiltrating / metastatic cells
        F = "numeric",           # a coefficient (Nmax = F/E) that determines the maximal number of cells
        # that can exist in the primary tumor when the hallmark is engraved
        c = "numeric",           # average number of divisions
        d = "numeric",           # mean value of splitting probability
        i = "numeric",           # average value of splitting probability
        a = "numeric",           # average value of apoptosis probability
        k = "numeric",           # average probability of cell death
        E = "numeric",           # average value of coefficients of friction term proportional to N, for splitting probability
        Nmax = "numeric",        # Maximal number of cells that can exist
        im = "numeric",          # average value of invasion / metastasis probability
        Ha = "numeric",
        Him = "numeric",
        Hi = "numeric",
        Hd = "numeric",
        Hb = "numeric",
        type = "numeric",        # invasion / metastatic ratio
        gene = "numeric",        # cancer gene damage rate
        # posdriver = "character", # cancer gene damage position (function deficit)
        # pospasngr = "character", # cancer gene damage position (maintaince of function)
        mutden = "numeric",      # density of gene malfunction
        last_id = "numeric"
    ),

    # Methods
    methods = list(
        # Initialize
        initialize = function(F0) {
            T <<- 0
            N <<- 0
            M <<- 0
            P <<- 0
            F <<- F0
        }
    )
)


#' Class 'OncoGene'
#'
#' @field id numeric. ID is same as in clone (key for clones)
#' @field name character. Onco genes' names list
#' @field onsp character. Oncogene/suppressor indicator for each gene in list of names
#' @field len numeric. Lengths of onco genes
#' @field cds_1 numeric. Onco genes' CDS base lengths for parental chr 1
#' @field cds_2 numeric. Onco genes' CDS base lengths for parental chr 2
#' @field rna_1 numeric. Onco genes RNA base number length for parental chr 1 (exons+introns)
#' @field rna_2 numeric. Onco genes RNA base number length for parental chr 2 (exons+introns)
#' @field p0_1 numeric. Probability of absent of mutations for parental chr 1
#' @field p0_2 numeric. Probability of absent of mutations for parental chr 2
#' @field prob_1 numeric. Vector of relative probabilities for point mutation, deletion and duplication: \cr
#' prob = c( m0 x sumCDS,  m_del x sumRNA,  m_dup x sumRNA ) / sum( m0 x sumCDS,  m_del x sumRNA,  m_dup x sumRNA )
#' @field prob_2 numeric.
#' @field sum_prob_1 numeric.
#' @field sum_prob_2 numeric.
#'
#' @return NULL
#' @export
#'
#' @examples
#' onco = tugHall_dataset$onco
#' onco$copy()
oncogene <- setRefClass(
    # name of the class
    Class = "OncoGene",

    # Fields
    fields = list(
        id = "numeric",       # identificator is same as in clone (key for clones)
        name  = "character",   # Cancer gene name list
        onsp  = "character",   # oncogene/suppressor indicator
        len   = "numeric",      # lengths of cancer genes
        cds_1 = "numeric",      # cancer gene CDS base length for parental chr 1
        cds_2 = "numeric",      # cancer gene CDS base length for parental chr 2
        rna_1 = "numeric",      # cancer gene RNA base number length for parental chr 1
        rna_2 = "numeric",      # cancer gene RNA base number length for parental chr 2
        p0_1  = "numeric",      # the probability of absent of mutations for parental chr 1
        p0_2  = "numeric",      # the probability of absent of mutations for parental chr 2
        prob_1 = "numeric",     # prob = c(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA) / sum(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA)
        prob_2 = "numeric",     # same for parental chr 2
        sum_prob_1 = "numeric",  # sum_prob = sum(m0*sumCDS,  m_del*sumRNA,  m_dup*sumRNA)
        sum_prob_2 = "numeric"  # same for parental chr 2
    ),

    # Methods
    methods = list(
        # read the configuration file
        read = function(file) {
            data = read.table(file, sep="\t", stringsAsFactors=FALSE )
            name0 = NULL
            onsp0 = NULL
            cds0 = NULL
            rna0 = NULL
            for (i in 1:nrow(data)) {
                name <<- as.character(data[i, 1])
                if (!is.element(name, name0)) {
                    name0 = c(name0, name)
                    type = as.character(data[i, 3])
                    if (type == "?") {
                        if (runif(1) > 0.5) {
                            type = "o"
                        } else {
                            type = "s"
                        }
                    }
                    onsp0 = c(onsp0, type)
                    # changed the length CDS, now get from gene_map file
                    # cds0 = c(cds0, as.numeric(as.character(data[i, 2])))
                }
            }

            # Get length of CDS and gene's length from gene_map:
            for ( i in 1:length(name0) ) {
                w     =  which( pck.env$gene_map$Gene == name0[ i ] )
                cds0  =  c( cds0, sum( pck.env$gene_map$End[w]  -  pck.env$gene_map$Start[w] + 1 ) )
                rna0  =  c( rna0, max( pck.env$gene_map$End[w] ) - min( pck.env$gene_map$Start[w]) + 1  )
            }

            id   <<- 1
            name <<- name0
            onsp <<- onsp0
            cds_1  <<- cds0
            cds_2  <<- cds0
            len  <<- length(name0)
            rna_1  <<- rna0
            rna_2  <<- rna0
            sum_prob_1 <<- sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
            sum_prob_2 <<- sum( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) )
            prob_1     <<- c( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob_1
            prob_2     <<- c( m0*sum(cds0),  m_del*sum(rna0),  m_dup*sum(rna0) ) / sum_prob_2
            p0_1   <<- (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)
            p0_2   <<- (1-m0)**sum(cds0) * (1-m_dup)**sum(rna0) * (1-m_del)**sum(rna0)
        }
    )
)


#' Class 'HallMark'
#'
#' @field Ha numeric. Apoptosis hallmark indexes of genes in onco$name, \cr where onco is object of class OncoGene
#' @field Hi numeric. Immortalization hallmark indexes of genes in onco$name, \cr where onco is object of class OncoGene
#' @field Hd numeric. Growth/antigrowth hallmark indexes of genes in onco$name, \cr where onco is object of class OncoGene
#' @field Hb numeric. Angiogenesis hallmark indexes of genes in onco$name, \cr where onco is object of class OncoGene
#' @field Him numeric. Invasion/metastatic transformation hallmark indexes of genes in onco$name, \cr where onco is object of class OncoGene
#' @field Ha_w numeric. Apoptosis hallmark weights of genes
#' @field Hi_w numeric. Immortalization hallmark weights of genes
#' @field Hd_w numeric. Growth/antigrowth hallmark weights of genes
#' @field Hb_w numeric. Angiogenesis hallmark weights of genes
#' @field Him_w numeric. Invasion/metastatic transformation hallmark weights of genes
#' @field notHa numeric. Indexes of genes which are not in apoptosis hallmark
#'
#' @return NULL
#' @export
#'
#' @examples
#' hall = tugHall_dataset$hall
#' print( hall )
#' hall$copy()
#' hall$show()
hallmark <- setRefClass(
    #
    Class = "HallMark",

    #
    fields = list(
        Ha = "numeric",       # (Evading apoptosis)
        Hi = "numeric",       # (immortalization limit)
        Hd = "numeric",       # (Insensitivity to anti-growth signals || Self-sufficiency in growth signals)
        Hb = "numeric",       # (Sustained angiogenesis)
        Him = "numeric",      # (Tissue invasion & metastasis)
        Ha_w = "numeric",     #
        Hi_w = "numeric",     #
        Hd_w = "numeric",     #
        Hb_w = "numeric",     #
        Him_w = "numeric",    #
        notHa = "numeric"
    ),

    #
    methods = list(
        #
        read = function(file, names, normalization  =  TRUE ) {
            # normalization is an indicator to normalize all Hallmarks values
            data <- read.table(file, sep="\t", stringsAsFactors=FALSE )
            Ha0 = NULL
            Hi0 = NULL
            Hd0 = NULL
            Hb0 = NULL
            Him0 = NULL
            Ha0_w = NULL
            Hi0_w = NULL
            Hd0_w = NULL
            Hb0_w = NULL
            Him0_w = NULL
            w_flg = FALSE
            if (ncol(data) >= 4) {
                w_flg = TRUE
            }
            # Acquire gene name and Hallmark coefficient by function from definition file
            for (i in 1:nrow(data)) {
                if (data[i, 2] == "apoptosis") {
                    Ha0 = c(Ha0, as.character(data[i, 1]))
                    if (w_flg) {
                        Ha0_w = c(Ha0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "immortalization") {
                    Hi0 = c(Hi0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hi0_w = c(Hi0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "anti-growth" | data[i, 2] == "growth") {
                    Hd0 = c(Hd0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hd0_w = c(Hd0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "angiogenesis") {
                    Hb0 = c(Hb0, as.character(data[i, 1]))
                    if (w_flg) {
                        Hb0_w = c(Hb0_w, as.numeric(as.character(data[i, 4])))
                    }
                } else if (data[i, 2] == "invasion") {
                    Him0 = c(Him0, as.character(data[i, 1]))
                    if (w_flg) {
                        Him0_w = c(Him0_w, as.numeric(as.character(data[i, 4])))
                    }
                }
            }
            Ha <<- match(Ha0, names)
            notHa <<- setdiff(seq(1,length(names)),Ha)
            Hi <<- match(Hi0, names)
            Hd <<- match(Hd0, names)
            Hb <<- match(Hb0, names)
            Him <<- match(Him0, names)

            # if there is no Hallmark coefficient then generate a Hallmark coefficient as a random number - beta distribution
            if (!w_flg) {
                if (length(Ha) > 0) {
                    Ha_rnd = 1:length(Ha)
                } else {
                    Ha_rnd = c()
                }
                total0 = length(Ha)
                if (length(Hi) > 0) {
                    Hi_rnd = (total0 + 1):(total0 + length(Hi))
                } else {
                    Hi_rnd = c()
                }
                total0 = total0 + length(Hi)
                if (length(Hd) > 0) {
                    Hd_rnd = (total0 + 1):(total0 + length(Hd))
                } else {
                    Hd_rnd = c()
                }
                total0 = total0 + length(Hd)
                if (length(Hb) > 0) {
                    Hb_rnd = (total0 + 1):(total0 + length(Hb))
                } else {
                    Hb_rnd = c()
                }
                total0 = total0 + length(Hb)
                if (length(Him) > 0) {
                    Him_rnd = (total0 + 1):(total0 + length(Him))
                } else {
                    Him_rnd = c()
                }
                total = total0 + length(Him)
                # random = runif(total)
                random = rbeta(total, 0.01, 1)
                Ha0_w = random[Ha_rnd]
                Hi0_w = random[Hi_rnd]
                Hd0_w = random[Hd_rnd]
                Hb0_w = random[Hb_rnd]
                Him0_w = random[Him_rnd]
            }
            # Total by genetic mode
            if ( normalization ){
                Ha_sum = sum(Ha0_w)
                Hi_sum = sum(Hi0_w)
                Hd_sum = sum(Hd0_w)
                Hb_sum = sum(Hb0_w)
                Him_sum = sum(Him0_w)
            } else {
                Ha_sum = 1
                Hi_sum = 1
                Hd_sum = 1
                Hb_sum = 1
                Him_sum = 1
            }

            Ha_w <<- Ha0_w / Ha_sum
            Hi_w <<- Hi0_w / Hi_sum
            Hd_w <<- Hd0_w / Hd_sum
            Hb_w <<- Hb0_w / Hb_sum
            Him_w <<- Him0_w / Him_sum

            if ( Compaction_factor ){
                Ha_w  <<-  pck.env$CF$Ha  * Ha_w
                Hi_w  <<-  pck.env$CF$Hi  * Hi_w
                Hd_w  <<-  pck.env$CF$Hd  * Hd_w
                Hb_w  <<-  pck.env$CF$Hb  * Hb_w
                Him_w <<-  pck.env$CF$Him * Him_w

            }

        },
        # Change the cell variables
        # mode = 2 Corresponding (Hallmark) Gene Mode
        updateClone = function(clone1, F) {
            # Apoptosis
            clone1$calcApoptosis()
            clone1$Ha  =  sum( clone1$gene[Ha] * Ha_w )
            clone1$a   =  clone1$a - clone1$Ha
            if ( clone1$a < 0 ) {
                clone1$a = 0
            }
            # Not dead - Immortalized
            clone1$Hi = sum( clone1$gene[Hi] * Hi_w )
            clone1$i = 1 - clone1$Hi
            if ( clone1$i < 0 ) {
                clone1$i = 0
            }
            # Angiogenesis
            clone1$Hb = sum( clone1$gene[Hb] * Hb_w )

            clone1$E     =   ifelse( !pck.env$tumbler_for_angiogenesis_trial, E0, ( E0  / (1 + F * clone1$Hb) ) )
            clone1$Nmax  =  1.0 / clone1$E

            # Cancer gene, tumor suppressor gene
            clone1$Hd  = sum( clone1$gene[Hd] * Hd_w )
            clone1$Him = sum( clone1$gene[Him] * Him_w )

            clone1$d = clone1$Hd + d0     # d0 is defined in parameters
            if ( clone1$d > 1 ) { clone1$d = 1 }
            if ( !clone1$invasion ) {
                clone1$d = clone1$d - clone1$E * pck.env$env$P
                if ( clone1$d < 0 ) { clone1$d = 0 }
            }

            # Invasion metastasis
            if ( !clone1$invasion ) {
                clone1$im = clone1$Him
            } else {
                clone1$im = 1
            }
        },
        # Change the environment variables
        updateEnviron = function(env, clones) {
            sum_cell(env, clones)
            # env$M = sum()  # ceiling(length(clones) * env$type)
            # env$N = sum()  # length(clones)  - env$M
        }
    )
)


#' Function to update Hallmark and variable after division or under initialization
#'
#' @param clone1 Object of class 'Clone'
#'
#' @return The same object of class 'Clone' with updated fields
#' @export
#'
#' @examples
#' clone = tugHall_dataset$clones[[ 1 ]]
#' load_tugHall.Environment( tugHall_dataset )
#' withr::with_environment( env = pck.env, code = update_Hallmarks( clone ) )
update_Hallmarks <- function(clone1) {
    # Hallmark
    pck.env$hall$updateClone(clone1, pck.env$env$F)
}


###  II) CNA and Point Mutations: CLASSES -------------------------------------------------

#' Class 'Point_Mutations'
#'
#' @field PointMut_ID numeric. ID of point mutation
#' @field Allele character. A or B allele
#' @field Parental_1or2 numeric. Parental chromosome, could be 1 or 2
#' @field Chr character. Chromosome name
#' @field Ref_pos numeric. Reference position
#' @field Phys_pos vector. Physical positions
#' @field Delta vector. Delta of positions
#' @field Copy_number numeric. Copy number of allele
#' @field Gene_name character. Gene's name
#' @field MalfunctionedByPointMut logical. True for driver mutation and False for passenger mutation
#' @field mut_order numeric. Number in order of mutation to reproduce the gene_map data.frame
#'
#' @return NULL
#' @export
#'
#' @examples
#' pnt = tugHall_dataset$pnt_clones[[ 1 ]]
#' print( pnt )
#' pnt$copy()
#' pnt$show()
#' pnt$initialize()
#' pnt$show()
#' pnt = tugHall_dataset$pnt_clones[[ 3 ]]
#' pnt$safe()   # save as row of data.frame
Point_Mutations <- setRefClass(
    #
    Class = "Point_Mutations",
    #
    fields = list(
        PointMut_ID	  = "numeric",         # ID of point mutation
        Allele        = "character",        # A or B allele
        Parental_1or2	= "numeric",         # 1 or 2
        Chr	          = "character",       # Chromosome name
        Ref_pos	      = "numeric",         # Reference position
        Phys_pos	    = "vector",          # Physical positions
        Delta	        = "vector",          # Delta of positions
        Copy_number	  = "numeric",         # Copy number of allele
        Gene_name	    = "character",       # Gene's name
        MalfunctionedByPointMut  = "logical",       # True or False
        mut_order     = "numeric"          # order of mutation to reproduce the map
    ),

    #
    methods = list(
        # Initialization with default values
        initialize = function( ID = 1, ... ){
            PointMut_ID	  <<- ID       # ID of point mutation
            Allele        <<- ""       # A or B allele
            Parental_1or2	<<-  0       # 1 or 2
            Chr	          <<-  ""       # Chromosome name
            Ref_pos	      <<-  0       # Reference position
            Phys_pos	    <<-  0       # Physical positions [from:to]
            Delta	        <<-  0       # Delta of positions [from:to]
            Copy_number	  <<-  0       # Copy number of allele
            Gene_name	    <<-  ""       # Gene's name
            MalfunctionedByPointMut  <<- NA       # True of False
            mut_order     <<-  0       # order of mutation
        },

        # Function to safe data to data frame df1
        safe = function() {
            df1 = data.frame( PointMut_ID = PointMut_ID,
                              Parental_1or2 = Parental_1or2,
                              Chr = Chr,
                              Ref_pos = Ref_pos,
                              Phys_pos = paste0( '[', paste(Phys_pos, collapse = ', ' ),  ']'),
                              Delta = paste0( '[', paste(Delta, collapse = ', ' ),  ']'),
                              Copy_number = Copy_number,
                              Gene_name = Gene_name,
                              MalfunctionedByPointMut = MalfunctionedByPointMut,
                              mut_order  =  mut_order,
                              stringsAsFactors = FALSE)
            # df = rbind( df, df1 )
            return( as.data.frame( df1 ) )
        }

    )
)



#' Class 'CNA_Mutations'
#'
#' @field CNA_ID numeric. ID of CNA mutation
#' @field Parental_1or2 numeric. Parental chromosome, could be 1 or 2
#' @field dupOrdel character. dup for duplication or del for deletion
#' @field Chr character. Chromosome name
#' @field Ref_start numeric. Reference start position
#' @field Ref_end numeric. Reference final position
#' @field Gene_names character. Names of genes involved in CNA
#' @field MalfunctionedByCNA logical. True for driver mutation and False for passenger mutation
#' @field mut_order numeric. Order of mutations in the lists of point mutations and CNA mutations
#'
#' @return NULL
#' @export
#'
#' @examples
#' cna = tugHall_dataset$cna_clones[[ 1 ]]
#' cna$safe()   # to save as row of data.frame
#' cna$copy()
#' cna$initialize()
#' cna$show()   # After initialization
CNA_Mutations <- setRefClass(
    #
    Class = "CNA_Mutations",
    #
    fields = list(
        CNA_ID	  = "numeric",         # ID of point mutation
        Parental_1or2	= "numeric",         # 1 or 2
        dupOrdel      = "character",         # duplication OR deletion indicator
        Chr	          = "character",       # Chromosome name
        Ref_start	    = "numeric",         # Reference start position
        Ref_end	      = "numeric",         # Reference end position
        # Phys_pos	    = "vector",          # Physical positions [from:to]
        # Delta	        = "vector",          # Delta of positions [from:to]
        # Copy_number	  = "numeric",         # Copy number of allele
        Gene_names	    = "character",       # Genes' name
        MalfunctionedByCNA  = "logical",     # True of False
        mut_order       = "numeric"          # Order of mutation
    ),

    #
    methods = list(
        # Initialization with default values
        initialize = function( ID = 1, ... ){
            CNA_ID	      <<-  ID       # ID of point mutation
            Parental_1or2	<<-  0       # 1 or 2
            dupOrdel      <<-  ""
            Chr	          <<-  ""       # Chromosome name
            Ref_start	    <<-  0       # Reference position
            Ref_end 	    <<-  0       #
            Gene_names	  <<- ""       # Gene's name
            MalfunctionedByCNA  <<- NA       # True of False
            mut_order     <<-  0
        },

        # Function to safe data to data frame df1
        safe = function() {
            df1 = data.frame( CNA_ID = CNA_ID,
                              Parental_1or2 = Parental_1or2,
                              dupOrdel = dupOrdel,
                              Chr = Chr,
                              Ref_start = Ref_start,
                              Ref_end = Ref_end,
                              Gene_names = Gene_names,
                              MalfunctionedByCNA = MalfunctionedByCNA,
                              mut_order  =  mut_order,
                              stringsAsFactors = FALSE)
            # df = rbind( df, df1 )
            return( as.data.frame( df1 ) )
        }

    )
)


