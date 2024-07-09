# task_handler.py

import json
import uuid
import shutil
import os
import subprocess
import argparse
from pathlib import Path
import logging
from enum import Enum
from mr_db import *
from utils import *

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
BASE_DIR = "./task_working_directories"


class TASK_STATUS(Enum):
    READY = "READY"
    PROCESSING = "PROCESSING"
    FINISHED = "FINISHED"
    FAILED = "FAILED"


class SCRIPT_MAP(Enum):
    BASIC_MR = "./simple_mr_processor.py"
    ZIP_DEMO = "./zip_processor.py"
    DATA_PURIFICATION = "./purification_processor.py"


def create_task_folder(base_dir):
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
    processor_script = SCRIPT_MAP[task_type].value
    task_folder = os.path.join("./task_working_directories", task_id)
    if task_type == SCRIPT_MAP.ZIP_DEMO.name:
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

    if task_type == SCRIPT_MAP.DATA_PURIFICATION.name:
        data_id = task_data.get("data_id")
        file_path = os.path.join("gwas_data", data_id + ".vcf.gz")
        shutil.copy(file_path, task_folder)

        # 库中插入
        new_task = MTask(
            id=task_id,
            name=f"{data_id} 净化",
            state=TASK_STATUS.PROCESSING.value,
            type=SCRIPT_MAP.DATA_PURIFICATION.name,
            input_params=task_data,
            description=f"对于 {task_data.get('data_id')} - {task_data.get('entity_name')} 所进行的净化任务"
        )
        mr_db.session.add(new_task)
        mr_db.session.commit()

        subprocess.run(["python", processor_script, "--gwas_id", task_data.get("data_id"), "--id",
                        task_id, "--entity_name", task_data.get("entity_name")], check=True)
        logging.info(
            f'Purification complete, output saved in folder {task_folder} and f"./purified_data/{data_id}.csv"')

    if task_type == SCRIPT_MAP.BASIC_MR.name:
        # 创建任务只需要传入两个gwas_id
        # 选择数据需要从数据列表进行选择

        exposure_id, outcome_id = task_data.get("exposure_id"), task_data.get("outcome_id")

        file_path = os.path.join("purified_data", exposure_id + '_exposure.csv')
        shutil.copy(file_path, task_folder)
        file_path = os.path.join("purified_data", outcome_id + '_outcome.csv')
        shutil.copy(file_path, task_folder)

        # 库中插入
        new_task = MTask(
            id=task_id,  # Ensure the correct field name 'id' is used
            name=task_data.get('name', f'{exposure_id} -> {outcome_id} 的简单MR分析'),
            state=TASK_STATUS.PROCESSING.value,
            type=SCRIPT_MAP.BASIC_MR.name,
            input_params=task_data,
            description=task_data.get('description',
                                      f"对于 {task_data.get('data_id')} - {task_data.get('entity_name')} 所进行的基础MR任务")
        )

        mr_db.session.add(new_task)
        mr_db.session.commit()

        # modify_to_outcome_csv(os.path.join(task_folder, f"{outcome_id}.csv"))
        subprocess.run(["python", processor_script, "--exposure_id", exposure_id, "--outcome_id",
                        outcome_id, "--task_id", task_id], check=True)

        logging.info(f'Basic MR completed, output saved in folder {task_folder}')


def create_task(task_data, base_dir=BASE_DIR):
    logging.info('Creating task folder based on input JSON data')
    task_folder, task_id = create_task_folder(base_dir)
    logging.info(f'Task folder created: {task_folder}')

    logging.info('Executing task')
    try:
        execute_task(task_id, task_data)
        mark_task_status(task_id, TASK_STATUS.FINISHED)
    except Exception as e:
        save_failed_task(task_id)
        # mark_task_status(task_id, TASK_STATUS.FAILED)
        # logging.error("task failed")

    return task_id

# if __name__ == "__main__":
#     # parser = argparse.ArgumentParser(description="Handle task based on input JSON and execute zip handling script.")
#     # parser.add_argument('--input-json', type=str, help='Path to the input JSON file.')
#     # parser.add_argument('--base-dir', type=str, required=False, help='Base directory to create task folders.')
#     # parser.add_argument('--zip-handler-script', type=str, required=False, help='Path to the zip_handler.py script.')
#     # args = parser.parse_args()
#     #
#     # args.zip_handler_script = "./demo_processor.py"
#     # args.input_json = "input_config.json"
#     # args.base_dir = "./task_working_directories"
#     #
#     # main(args.input_json, args.base_dir, args.zip_handler_script)
#     task_req = {"data_id": "ieu-a-2", "task_type": "DATA_PURIFICATION", "entity_name": "BMI"}
#     create_task(task_req)
