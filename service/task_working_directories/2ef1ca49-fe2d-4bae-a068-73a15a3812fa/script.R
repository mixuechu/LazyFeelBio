
    library(VariantAnnotation)
    library(gwasglue)
    library(TwoSampleMR)

    library(MendelianRandomization)
    library(dplyr)
    library(tidyr)
    library(CMplot)

    options(download.file.method = "curl")
    options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
    setwd("./task_working_directories/2ef1ca49-fe2d-4bae-a068-73a15a3812fa")

    inputFile <- "ieu-a-2.vcf.gz"

    vcfRT <- readVcf(inputFile)
    data <- gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")
    outTab <- subset(data, pval.exposure < 5e-08)
    exposure_dat_clumped <- clump_data(outTab, clump_kb=10000, clump_r2=0.001)

    Ffilter <- 10
    N <- exposure_dat_clumped[1, "samplesize.exposure"]
    dat <- transform(exposure_dat_clumped, R2 = 2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
    dat <- transform(dat, F = (N-2)*R2/(1-R2))
    outTab <- dat[dat$F > Ffilter,]
    # 此处存疑
    write.csv(outTab, "ieu-a-2.csv", row.names=FALSE)
    