
## Simulate a simple F2 cross and try linkage mapping and GWAS on it

library(AlphaSimR)
library(assertthat)
library(genetics)
library(magrittr)
library(qtl)
library(rrBLUP)

source("R/simulation_functions.R")




system("mkdir simulations")

base_dir <- getwd()

for (i in 1:20) {

    system(paste("mkdir simulations/f2_inbred", i, sep = ""))
    setwd(paste(base_dir,
                "/simulations/f2_inbred",
                i, sep = ""))
    
    mapping_population <- setup_population("f2_inbred")
    SIMPARAM <- mapping_population$SIMPARAM
    population <- mapping_population$founder_population
    f2 <- mapping_population$mapping_population    
    
    snps_dense <- pullSnpGeno(f2, 1, simParam = SIMPARAM)
    snps_sparse <- pullSnpGeno(f2, 2, simParam = SIMPARAM)
    
    map_dense <- pull_snp_map(SIMPARAM, 1)
    map_dense$cM <- map_dense$loc * 100
    map_dense$marker_name <- colnames(snps_dense)
    
    map_sparse <- pull_snp_map(SIMPARAM, 2)
    map_sparse$cM <- map_sparse$loc * 100
    map_sparse$marker_name <- colnames(snps_sparse)
    
    
    ## R/qtl
    
    parental_snps_sparse <- pullSnpGeno(population, 2, simParam = SIMPARAM)
    
    recoded_sparse <- recode_genotypes(snps_sparse, parental_snps_sparse)
    
    pheno_geno_sparse <- data.frame(id = 1:f2@nInd, trait = f2@pheno[,1], recoded_sparse)
    
    
    rqtl_file <-  file("rqtl_sparse.csv", "w")
    writeLines(paste(colnames(pheno_geno_sparse), collapse = ","), rqtl_file)
    writeLines(paste(c("","", map_sparse$chr), collapse = ","), rqtl_file)
    writeLines(paste(c("","", map_sparse$cM), collapse = ","), rqtl_file)
    write.table(pheno_geno_sparse,
                sep = ",",
                file = rqtl_file,
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE)
    close(rqtl_file)
    
    
    cross <- read.cross(file = "rqtl_sparse.csv",
                        format = "csv",
                        genotypes = c(1, 2, 3, 4, 5))
    
    scan <- scanone(cross, pheno.col = 2)
    perm <- scanone(cross, pheno.col = 2, n.perm = 500, n.cluster = 7)
    threshold <- summary(perm)[1,1]
    hits <- summary(scan, threshold = threshold)
    hit_intervals <- lapply(hits$chr, lodint, results = scan, drop = 1.8, expandtomarkers = TRUE)
    
    
    linkage_check <- check_linkage_results(hit_intervals,
                                           SIMPARAM)
    qtl_locations <- linkage_check$qtl_locations
    false_positives_linkage <- linkage_check$false_positives_linkage    
    
    
    
    ## rrBLUP
    
    pheno <- data.frame(gid = paste("ind", 1:f2@nInd, sep = ""),
                        trait = f2@pheno[,1],
                        stringsAsFactors = FALSE)
    
    geno <- data.frame(marker_name = colnames(snps_sparse),
                       map_sparse[, c("chr", "cM")],
                       t(snps_sparse - 1),
                       stringsAsFactors = FALSE)
    colnames(geno)[4:ncol(geno)] <- pheno$gid
    
    
    gwas <- GWAS(pheno, geno, P3D = TRUE, n.core = 7)
    gwas$p <- 10^-gwas$trait
    gwas$p_fdr <- p.adjust(gwas$p, "fdr")
    
    
    gwas_hits <- subset(gwas, p_fdr < 0.05)
    
    gwas_check <- check_gwas_results(gwas_hits,
                                     qtl_locations,
                                     snps_sparse,
                                     map_sparse)
    qtl_locations <- gwas_check$qtl_locations
    false_positives_gwas <- gwas_check$false_positives_gwas
    ld_markers <- gwas_check$ld_markers
    
    ## Summarise true and false positives
    
    true_false_positives <- data.frame(detected_linkage = sum(qtl_locations$detected_linkage),
                                       detected_gwas = sum(qtl_locations$detected_gwas),
                                       false_positives_linkage,
                                       false_positives_gwas)
    
    
    
    saveRDS(qtl_locations, file = "qtl_locations.Rds")
    saveRDS(hits, file = "linkage_hits.Rds")
    saveRDS(hit_intervals, file = "linkage_intervals.Rds")
    saveRDS(gwas_hits, file = "gwas_hits.Rds")
    saveRDS(ld_markers, file = "ld_markers.Rds")
    saveRDS(true_false_positives, file = "true_false_positives.Rds")


    setwd(base_dir)
}
