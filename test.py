import subprocess


def download_exposure(exposure_name):
    pass

def download_outcome(outcome_name):
    pass

def run_process(input_file, outcome_file):
    subprocess.run(["Rscript", "whole_process.r", input_file, outcome_file])


def test_single_func():
    subprocess.run(["Rscript", "test.R"])

if __name__ == '__main__':
    download_exposure("exposure")
    download_outcome("outcome")

    test_single_func()

    # run_process("ieu-a-2.vcf.gz", "ieu-a-7.vcf.gz")



