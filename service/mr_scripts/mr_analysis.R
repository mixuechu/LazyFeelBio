
    dat <- read.csv("table.SNP.csv")
    presso <- run_mr_presso(dat)
    write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file="table.MR-PRESSO.csv")
    mrResult <- mr(dat)
    mrTab <- generate_odds_ratios(mrResult)
    write.csv(mrTab, file="table.MRresult.csv", row.names=FALSE)
    heterTab <- mr_heterogeneity(dat)
    write.csv(heterTab, file="table.heterogeneity.csv", row.names=FALSE)
    pleioTab <- mr_pleiotropy_test(dat)
    write.csv(pleioTab, file="table.pleiotropy.csv", row.names=FALSE)

    pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
    mr_scatter_plot(mrResult, dat)
    dev.off()

    res_single <- mr_singlesnp(dat)
    pdf(file="pic.forest.pdf", width=7, height=6.5)
    mr_forest_plot(res_single)
    dev.off()

    pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
    mr_funnel_plot(singlesnp_results = res_single)
    dev.off()

    pdf(file="pic.leaveoneout.pdf", width=7, height=6.5)
    mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    dev.off()
    