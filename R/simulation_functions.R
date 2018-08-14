

## Function to set up base population

setup_population <- function(scenario) {

    if (! scenario %in% c("f2_diverse", "f2_inbred",
                          "f8_diverse", "f8_inbred")) {
        stop("Unknown scenario")
    }
    
    if (scenario == "f2_diverse" | scenario == "f8_diverse") {
        founderpop <- runMacs(segSites = 5000,
                              nInd = 8,
                              nChr = 10,
                              species = "CATTLE",
                              split = 10000,
                              inbred = FALSE)
        
    } else {
        founderpop <- runMacs(segSites = 5000,
                              nInd = 2,
                              nChr = 10,
                              species = "CATTLE",
                              split = 10000,
                              inbred = TRUE)
    }

    SIMPARAM <- SimParam$new(founderpop)
    SIMPARAM$setGender("yes_sys")
    SIMPARAM$restrSegSites(maxQtl = 1000,
                           maxSnp = 1000,
                           snpQtlOverlap = FALSE)
    SIMPARAM$addTraitA(nQtlPerChr = 1,
                       mean = 100,
                       var = 10)
    SIMPARAM$addSnpChip(1000)
    SIMPARAM$addSnpChip(100)

    population <- newPop(founderpop,
                         simParam = SIMPARAM) %>%
        setPheno(varE = 20,
                 simParam = SIMPARAM)


    ## Inbred F2
    if (scenario == "f2_inbred") {
        print("A")
        f1 <- randCross(population,
                        nCrosses = 100,
                        simParam = SIMPARAM)
        print("B")
        f2 <- randCross(f1,
                        nCrosses = 100,
                        nProgeny = 10,
                        simParam = SIMPARAM) %>%
            setPheno(varE = 20,
                     simParam = SIMPARAM)
        mapping_population <- f2
    }


    ## Inbred F8
    
    if (scenario == "f8_inbred") {        
        f1 <- randCross(population,
                        nCrosses = 100,
                        simParam = SIMPARAM)

        generations <- vector(mode = "list",
                              length = 8)
        generations[[1]] <- f1
        
        for (i in 2:7) {
            generations[[i]] <- randCross(generations[[i - 1]],
                                          nCrosses = 10,
                                          nProgeny = 10,
                                          simParam = SIMPARAM)
        }
        generations[[8]] <- randCross(generations[[7]],
                                      nCrosses = 100,
                                      nProgeny = 10,
                                      simParam = SIMPARAM) %>%
            setPheno(varE = 20,
                     simParam = SIMPARAM)
        mapping_population <- generations[[8]]
    }

    
    ## Diverse F8
    
    if (scenario == "f8_diverse") {
        pop1 <- population[1:4]
        pop2 <- population[5:8]
        
        f1 <- c(randCross2(females = pop1[pop1@gender == "F"],
                           males = pop2[pop2@gender == "M"],
                           nCrosses = 50),
                randCross2(females = pop2[pop2@gender == "F"],
                           males = pop1[pop1@gender == "M"],
                           nCrosses = 50))
        
        generations <- vector(mode = "list",
                              length = 8)
        generations[[1]] <- f1
        
        for (i in 2:7) {
            generations[[i]] <- randCross(generations[[i - 1]],
                                          nCrosses = 10,
                                          nProgeny = 10)
        }
        generations[[8]] <- randCross(generations[[7]],
                                      nCrosses = 100,
                                      nProgeny = 10) %>%
            setPheno(varE = 20,
                     simParam = SIMPARAM)
        mapping_population <- generations[[8]]
    }

    ## Diverse F2
    if (scenario == "f2_diverse") {
        pop1 <- population[1:4]
        pop2 <- population[5:8]
        
        f1 <- c(randCross2(females = pop1[pop1@gender == "F"],
                           males = pop2[pop2@gender == "M"],
                           nCrosses = 50,
                           simParam = SIMPARAM),
                randCross2(females = pop2[pop2@gender == "F"],
                           males = pop1[pop1@gender == "M"],
                           nCrosses = 50,
                           simParam = SIMPARAM))
        
        f2 <- randCross(f1,
                        nCrosses = 100,
                        nProgeny = 10,
                        simParam = SIMPARAM) %>%
            setPheno(varE = 20,
                     simParam = SIMPARAM)

        mapping_population <- f2
    }


    list(SIMPARAM = SIMPARAM,
         founder_population = population,
         mapping_population = mapping_population)
}



## Functions to find QTL locations

pull_snp_map <- function(simparam, chip_no = 1) {
    n_chr <- simparam$nChr
    maps <- vector(length = n_chr, mode = "list")
    for (chr_ix in 1:n_chr) {
        chr_map <- simparam$genMaps[chr_ix,][[1]]
        start_ix <- sum(simparam$snpChips[[chip_no]]@lociPerChr[1:chr_ix]) -
            simparam$snpChips[[chip_no]]@lociPerChr[chr_ix] + 1
        end_ix <- sum(simparam$snpChips[[chip_no]]@lociPerChr[1:chr_ix])
        snp_ix <- simparam$snpChips[[chip_no]]@lociLoc[start_ix:end_ix]
        maps[[chr_ix]] <- data.frame(chr = chr_ix,
                                     loc = chr_map[snp_ix])
    }
    Reduce(rbind, maps)
}

pull_qtl_map <- function(simparam) {
    n_chr <- simparam$nChr
    maps <- vector(length = n_chr, mode = "list")
    for (chr_ix in 1:n_chr) {
        chr_map <- simparam$genMaps[chr_ix,][[1]]
        start_ix <- sum(simparam$traits[[1]]@lociPerChr[1:chr_ix]) -
            simparam$traits[[1]]@lociPerChr[chr_ix] + 1
        end_ix <- sum(simparam$traits[[1]]@lociPerChr[1:chr_ix])
        snp_ix <- simparam$traits[[1]]@lociLoc[start_ix:end_ix]
        maps[[chr_ix]] <- data.frame(chr = chr_ix,
                                     loc = chr_map[snp_ix])
    }
    Reduce(rbind, maps)
}


window_locations <- function(map, window_size) {
    n_snps <- nrow(map)
    n_windows <- ceiling(n_snps/window_size)
    
    window_start <- seq(from = 1,
                        to = n_windows * window_size,
                        by = window_size)
    window_end <- window_start + window_size
    if (window_end[n_windows] > n_snps) {
        window_end[n_windows] <- n_snps
    }
    
    data.frame(chr = map$chr[window_start],
               start = map$loc[window_start],
               end = map$loc[window_end])
}



qtl_in_windows <- function(window_locations,
                           qtl_locations) {

    qtl_locations$id <- 1:nrow(qtl_locations)

    window_locations$contains_qtl <- FALSE
    for (window_ix in 1:nrow(window_locations)) {
        qtl_in_window <- subset(qtl_locations,
                                chr == window_locations$chr[window_ix] &
                                loc >= window_locations$start[window_ix] &
                                loc < window_locations$end[window_ix])
        if (nrow(qtl_in_window) > 0) {
            print(qtl_in_window)
            window_locations$contains_qtl[window_ix] <- TRUE
        }
    }

    window_locations
}


## Recode genotypes for R/qtl

recode_genotypes <- function(snps,
                             parental_snps) {

    codeA <- parental_snps[1,] == 0
    codeB <- parental_snps[2,] == 0
    assert_that(all(codeA != codeB))
    
    recoded <- snps + 1
    for (snp_ix  in 1:ncol(snps)) {
        if (codeB[snp_ix]) {
            change_to_3 <- which(recoded[, snp_ix] == 1)
            change_to_1 <- which(recoded[, snp_ix] == 3)
            recoded[change_to_3, snp_ix] <- 1
            recoded[change_to_1, snp_ix] <- 3
        }
    }
    recoded
}


recode_genotypes_diverse <- function(snps,
                                     parental_snps) {
    codeA <- logical(ncol(parental_snps))
    codeB <- logical(ncol(parental_snps))
    for (snp_ix in 1:ncol(parental_snps)) {
        codeA[snp_ix] <- codeB[snp_ix] <- NA
        if (all(parental_snps[1:4, snp_ix] == 0)) {
            codeA[snp_ix] <- TRUE
        }
        if (all(parental_snps[5:8, snp_ix] == 0)) {
            codeB[snp_ix] <- TRUE
        }
        if (all(parental_snps[1:4, snp_ix] == 2)) {
            codeA[snp_ix] <- FALSE
        }
        if (all(parental_snps[5:8, snp_ix] == 2)) {
            codeB[snp_ix] <- FALSE
        }
    }
    assert_that(all(na.exclude(codeA != codeB)))

    recoded <- snps + 1
    for (snp_ix  in 1:ncol(snps)) {
        if (!is.na(codeA[snp_ix]) &
            !is.na(codeB[snp_ix]) &
            codeB[snp_ix]) {
            change_to_3 <- which(recoded[, snp_ix] == 1)
            change_to_1 <- which(recoded[, snp_ix] == 3)
            recoded[change_to_3, snp_ix] <- 1
            recoded[change_to_1, snp_ix] <- 3
        }
        if (is.na(codeA[snp_ix]) | is.na(codeB[snp_ix])) {
            recoded[, snp_ix] <- NA
        }
    }
    recoded[,!is.na(colSums(recoded))]
}


## Calculate LD with a target marker (needs genetics package)

get_target_ld <- function(geno,
                          focal_snp) {
    
    focal_geno <- as.genotype.allele.count(geno[, focal_snp])

    r2 <- numeric(ncol(geno))
    for(col_ix in 1:ncol(geno)) {
        r2[col_ix] <- LD(focal_geno,
                         as.genotype.allele.count(geno[, col_ix]))$"R^2"
    }
    names(r2) <- colnames(geno)
    
    r2
}


## Check linkage results

check_linkage_results <- function(hit_intervals,
                                  SIMPARAM) {
    
    qtl_locations <- pull_qtl_map(SIMPARAM)
    qtl_locations$cM <- qtl_locations$loc * 100
    qtl_locations$detected_linkage <- FALSE
    false_positives_linkage <- 0

    if (length(hit_intervals) > 0) {
        for (hit_ix in 1:length(hit_intervals)) {
            chr <- unique(as.numeric(as.character(hit_intervals[[hit_ix]]$chr)))
            detected <- qtl_locations$chr == chr &
                qtl_locations$cM >= hit_intervals[[hit_ix]]$pos[1] &
                qtl_locations$cM <= hit_intervals[[hit_ix]]$pos[3]
            qtl_locations$detected_linkage[which(detected)] <-  TRUE
            if (all(detected == FALSE)) {
                false_positives_linkage <- false_positives_linkage + 1
            }
        }
    }
    list(qtl_locations = qtl_locations,
         false_positives_linkage = false_positives_linkage)
}


## Check gwas results (building on the same data frame)

check_gwas_results <- function(gwas_hits,
                               qtl_locations,
                               snps,
                               map) {
    qtl_locations$detected_gwas <- FALSE
    ld_markers <- vector(mode = "list", length = nrow(gwas_hits))
    false_positives_gwas <- 0
    
    false_positives_list <- vector(mode = "list", length = nrow(gwas_hits))
    false_positive_counter <- 1
    if (nrow(gwas_hits) > 0) {
        for (hit_ix in 1:nrow(gwas_hits)) {
            
            chr_markers <- subset(map, chr == gwas_hits$chr[hit_ix])
            chr_geno <- snps[, colnames(snps) %in% chr_markers$marker_name]
            
            r2 <- get_target_ld(chr_geno, gwas_hits$marker_name[hit_ix])
            
            ld_markers[[hit_ix]] <- subset(chr_markers, marker_name %in% names(which(r2 > 0.2)))
            
            detected <- qtl_locations$chr == gwas_hits$chr[hit_ix] &
                qtl_locations$cM >= min(ld_markers[[hit_ix]]$cM) &
                qtl_locations$cM <= max(ld_markers[[hit_ix]]$cM)
            qtl_locations$detected_gwas[which(detected)] <-  TRUE
            if (all(detected == FALSE)) {
                ##                false_positives_gwas <- false_positives_gwas + 1
                false_positives_list[[false_positive_counter]] <-
                    gwas_hits[hit_ix,]
                false_positive_counter <- false_positive_counter + 1
            }
        }
    }
    false_positives_df <- Reduce(rbind, false_positives_list)
    false_positives_gwas <- length(unique(false_positives_df$chr))

    list(qtl_locations = qtl_locations,
         false_positives_gwas = false_positives_gwas,
         ld_markers = ld_markers)
}
