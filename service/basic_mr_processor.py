# basic_mr_processor.py

import os.path
from enum import Enum
import subprocess
import shutil

import argparse

from utils import *


class MRProcess(Enum):
    READ_VCF = "read_vcf.R"
    CLUMP = "clump.R"
    FILTER_BY_F = "filter_by_f.R"
    HARMONIZE = "harmonize.R"
    MR_ANALYSIS = "mr_analysis.R"
    PLOT_MANHATTAN = "plot_manhattan.R"


def create_and_run_r_script(input_param):
    working_directory = f"./task_working_directories/{input_param.get('id')}"
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
    inputFile <- "{input_param.get('exposure_file')}"
    outcomeFile <- "{input_param.get('outcome_file')}"
    exposureName <- "{input_param.get('exposure_name')}"
    outcomeName <- "{input_param.get('outcome_name')}"
    """

    for process in input_param.get('processes'):
        script_path = os.path.join("mr_scripts", MRProcess[process].value)
        module_content = read_r_script(script_path)
        script_content += "\n" + module_content
    complete_script_path = f"{working_directory}/script.R"
    with open(complete_script_path, 'w', encoding='utf-8') as file:
        file.write(script_content)

    subprocess.run(['Rscript', complete_script_path], check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="do the most basic mendelian randomization analysis.")

    parser.add_argument('--exposure_id', required=True, type=str, help='暴露数据ID')
    parser.add_argument('--outcome_id', required=True, type=str, help='结局数据ID')
    parser.add_argument('--task_id', required=True, type=str, help='任务ID')
    args = parser.parse_args()

    input_param = {
        "exposure_id": args.exposure_id,
        "outcome_id": args.outcome_id,
        "task_id": args.task_id
    }

    create_and_run_r_script(input_param)