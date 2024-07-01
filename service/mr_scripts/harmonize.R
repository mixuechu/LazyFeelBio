
    library(VariantAnnotation)
    library(TwoSampleMR)
    exposure_dat <- read.csv("exposure.F.csv")
    vcfRT <- readVcf(outcomeFile)
    outcomeData <- gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="outcome")
    outcomeTab <- merge(exposure_dat, outcomeData, by.x="SNP", by.y="SNP")
    write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv", row.names=FALSE)
    outcome_dat <- read_outcome_data(snps=exposure_dat$SNP,
                     filename="outcome.csv", sep = ",",
                     snp_col = "SNP",
                     beta_col = "beta.outcome",
                     se_col = "se.outcome",
                     effect_allele_col = "effect_allele.outcome",
                     other_allele_col = "other_allele.outcome",
                     pval_col = "pval.outcome",
                     eaf_col = "eaf.outcome")
    exposure_dat$exposure <- "BMI"
    outcome_dat$outcome <- "Coronary heart disease"
    dat <- harmonise_data(exposure_dat=exposure_dat, outcome_dat=outcome_dat)
    outTab <- dat[dat$mr_keep=="TRUE",]
    write.csv(outTab, file="table.SNP.csv", row.names=FALSE)
    