import json
import uuid
import shutil
import os
import subprocess
import argparse
from pathlib import Path
import logging
from enum import Enum

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
BASE_DIR = "./task_working_directories"


class ScriptMap(Enum):
    BASIC_MR = "./basic_mr_processor.py"
    ZIP_DEMO = "./zip_processor.py"


def create_task_folder(task_data, base_dir):
    """
    根据输入的任务数据创建以 UUID 命名的任务文件夹，并复制需要的文件。

    Args:
    task_data (dict): 输入的任务数据，包含需要的数据文件。
    base_dir (str): 基础目录，所有任务文件夹将创建在此目录下。

    Returns:
    str: 创建的任务文件夹路径。
    """
    task_uuid = str(uuid.uuid4())
    task_folder = Path(base_dir) / task_uuid
    task_folder.mkdir(parents=True, exist_ok=True)

    return task_folder, task_uuid


def execute_task(task_id, task_data):
    task_type = task_data.get("task_type", "BASIC_MR")
    processor_script = ScriptMap[task_type].value
    task_folder = os.path.join("./task_working_directories", task_id)
    if task_type == ScriptMap.ZIP_DEMO.name:
        # 复制文件到任务文件夹
        for file_path in task_data['files']:
            file_path = os.path.join("gwas_data", file_path)
            shutil.copy(file_path, task_folder)

        # 获取任务文件夹中的所有 zip 文件
        zip_files = [str(p) for p in Path(task_folder).glob("*.zip")]
        output_zip = Path(task_folder) / "output.zip"

        # 调用 zip_handler.py 脚本
        subprocess.run(["python", processor_script, *zip_files, "-o", str(output_zip)], check=True)
        logging.info(f'Task execution completed, output saved to: {task_folder / "output.zip"}')

    if task_type == ScriptMap.BASIC_MR.name:
        outcome_file_name, exposure_file_name = task_data.get("outcome_file"), task_data.get("exposure_file")

        file_path = os.path.join("gwas_data", outcome_file_name)
        shutil.copy(file_path, task_folder)
        file_path = os.path.join("gwas_data", exposure_file_name)
        shutil.copy(file_path, task_folder)

        subprocess.run(["python", processor_script, "--outcome_name", task_data.get("outcome_name"), "--outcome_file",
                        outcome_file_name, "--exposure_name", task_data.get("exposure_name"), "--id", task_id,
                        "--exposure_file", exposure_file_name], check=True)
        logging.info(f'Basic MR completed, output saved in folder {task_folder}')


def create_task(task_data, base_dir=BASE_DIR):
    logging.info('Creating task folder based on input JSON data')
    task_folder, task_id = create_task_folder(task_data, base_dir)
    logging.info(f'Task folder created: {task_folder}')

    logging.info('Executing task')
    execute_task(task_id, task_data)
    return task_id

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Handle task based on input JSON and execute zip handling script.")
#     parser.add_argument('--input-json', type=str, help='Path to the input JSON file.')
#     parser.add_argument('--base-dir', type=str, required=False, help='Base directory to create task folders.')
#     parser.add_argument('--zip-handler-script', type=str, required=False, help='Path to the zip_handler.py script.')
#     args = parser.parse_args()
#
#     args.zip_handler_script = "./demo_processor.py"
#     args.input_json = "input_config.json"
#     args.base_dir = "./task_working_directories"
#
#     main(args.input_json, args.base_dir, args.zip_handler_script)