# simple_mr_processor.py
import subprocess
import argparse
import shutil
import time
from utils import *


def process(input_param):
    exposure_id = input_param.get('exposure_id')
    outcome_id = input_param.get('outcome_id')

    exposure_name = get_data_entity_name(exposure_id)
    outcome_name = get_data_entity_name(outcome_id)

    working_directory = f"./task_working_directories/{input_param.get('task_id')}"
    script_content = f"""
    library(VariantAnnotation)
    library(gwasglue)
    library(TwoSampleMR)

    library(MendelianRandomization)
    library(dplyr)
    library(tidyr)
    library(CMplot)

    options(download.file.method = "curl")
    options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
    setwd("{working_directory}")
    exposure_dat <- read.csv("{exposure_id}_purified.csv")
    # vcfRT <- readVcf("{outcome_id}.vcf.gz")
    # outcomeData <- gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="outcome")
    
    
    # outcomeData <- read.csv("{outcome_id}_original.csv")
    # outcomeTab <- merge(exposure_dat, outcomeData, by.x="SNP", by.y="SNP")
    # write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv", row.names=FALSE)
    outcome_dat <- read_outcome_data(snps=exposure_dat$SNP,
                     filename="outcome.csv", sep = ",",
                     snp_col = "SNP",
                     beta_col = "beta.outcome",
                     se_col = "se.outcome",
                     effect_allele_col = "effect_allele.outcome",
                     other_allele_col = "other_allele.outcome",
                     pval_col = "pval.outcome",
                     eaf_col = "eaf.outcome")
    exposure_dat$exposure <- "{exposure_name}"
    outcome_dat$outcome <- "{outcome_name}"
    dat <- harmonise_data(exposure_dat=exposure_dat, outcome_dat=outcome_dat)
    outTab <- dat[dat$mr_keep=="TRUE",]
    write.csv(outTab, file="table.SNP.csv", row.names=FALSE)
    
    dat <- outTab
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
    
    """

    complete_script_path = f"{working_directory}/script.R"

    with open(complete_script_path, 'w', encoding='utf-8') as file:
        file.write(script_content)

    start_time = time.time()

    shutil.copy(os.path.join(working_directory, f"{outcome_id}_original.csv"),
                os.path.join("./purified_data", f"{outcome_id}_original.csv"))
    shutil.copy(os.path.join(working_directory, f"{exposure_id}_purified.csv"),
                os.path.join("./purified_data", f"{exposure_id}_purified.csv"))

    print(f"搬运文件花费的时间: {time.time() - start_time:.2f} 秒")
    start_time = time.time()
    # modify_to_outcome_csv(os.path.join(working_directory, f"{outcome_id}_original.csv"))
    print(f"转表头花费时间: {time.time() - start_time:.2f} 秒")

    start_time = time.time()
    subprocess.run(['Rscript', complete_script_path], check=True)
    print(f"simple mr cost: {time.time() - start_time:.2f} 秒")
    # shutil.copy(os.path.join(working_directory, f"{gwas_id}.csv"), os.path.join("./purified_data", f"{gwas_id}.csv"))
    # import time
    # time.sleep(30)

    # if check_purified(gwas_id):
    #     print(f"PURIFICATION Data purified, result saved to ./purified_data/{gwas_id}.csv")
    #     mark_data_status(gwas_id, DataStatus.PURIFIED)
    # else:
    #     mark_data_status(gwas_id, DataStatus.FAILED)
    #     raise


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(
    #     description="使用提纯过的数据进行快速孟德尔随机化分析")
    #
    # parser.add_argument('--exposure_id', required=True, type=str, help='暴露实体的GWAS ID')
    # parser.add_argument('--outcome_id', required=True, type=str, help='结局实体的GWAS ID')
    # parser.add_argument('--task_id', required=True, type=str, help='任务ID')
    # args = parser.parse_args()
    #
    # input_param = {
    #     "exposure_id": args.exposure_id,
    #     "outcome_id": args.outcome_id,
    #     "task_id": args.task_id
    # }

    input_param = {
        "exposure_id": "ieu-a-2",
        "outcome_id": "ieu-a-7",
        "task_id": "demo"
    }

    process(input_param)
