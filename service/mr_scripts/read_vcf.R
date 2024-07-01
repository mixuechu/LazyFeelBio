
    vcfRT <- readVcf(inputFile)
    data <- gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")
    outTab <- subset(data, pval.exposure < 5e-08)
    write.csv(outTab, file="exposure.pvalue.csv", row.names=FALSE)
