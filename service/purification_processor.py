import subprocess
import argparse
import shutil
from utils import *


def process(input_param):
    gwas_id = input_param.get('gwas_id')
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

    inputFile <- "{gwas_id}.vcf.gz"

    vcfRT <- readVcf(inputFile)
    outTab <- gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")
    write.csv(outTab, "{gwas_id}_original.csv", row.names=FALSE)
    outTab <- subset(outTab, pval.exposure < 5e-08)
    exposure_dat_clumped <- clump_data(outTab, clump_kb=10000, clump_r2=0.001)

    Ffilter <- 10
    N <- exposure_dat_clumped[1, "samplesize.exposure"]
    dat <- transform(exposure_dat_clumped, R2 = 2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
    dat <- transform(dat, F = (N-2)*R2/(1-R2))
    outTab <- dat[dat$F > Ffilter,]
    # 此处存疑
    write.csv(outTab, "{gwas_id}_purified.csv", row.names=FALSE)
    """

    complete_script_path = f"{working_directory}/script.R"

    mark_data_status(gwas_id, DataStatus.PURIFYING)

    with open(complete_script_path, 'w', encoding='utf-8') as file:
        file.write(script_content)

    subprocess.run(['Rscript', complete_script_path], check=True)
    shutil.copy(os.path.join(working_directory, f"{gwas_id}_original.csv"),
                os.path.join("./purified_data", f"{gwas_id}_original.csv"))
    shutil.copy(os.path.join(working_directory, f"{gwas_id}_purified.csv"),
                os.path.join("./purified_data", f"{gwas_id}_purified.csv"))
    # import time
    # time.sleep(30)

    if check_purified(gwas_id):
        print(f"PURIFICATION Data purified, result saved to ./purified_data/{gwas_id}_purified.csv")
        mark_data_status(gwas_id, DataStatus.PURIFIED)
    else:
        mark_data_status(gwas_id, DataStatus.FAILED)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="purify exposure data for future use")
    parser.add_argument('--gwas_id', required=True, type=str, help='Gwas ID of the exposure.')
    parser.add_argument('--id', required=True, type=str, help='Task id.')
    parser.add_argument('--entity_name', required=True, type=str, help='实体名.')
    args = parser.parse_args()

    input_param = {
        "gwas_id": args.gwas_id,
        "task_id": args.id,
        "entity_name": args.entity_name
    }

    process(input_param)
