
    data <- read.csv("exposure.pvalue.csv")
    exposure_dat_clumped <- clump_data(data, clump_kb=10000, clump_r2=0.001)
    write.csv(exposure_dat_clumped, file="exposure.LD.csv", row.names=FALSE)

