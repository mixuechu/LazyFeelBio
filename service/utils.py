# utils.py

import os
from flask import jsonify
from sqlalchemy.exc import SQLAlchemyError
from mr_db import *
from sqlalchemy import text
from enum import Enum
import pandas as pd


class GWASDataExistsError(Exception):
    pass


class DataStatus(Enum):
    RAW = "RAW"
    PURIFIED = "PURIFIED"
    PURIFYING = "PURIFYING"
    MISSING = "MISSING"
    FAILED = "FAILED"


def modify_to_outcome_csv(file_path):
    df = pd.read_csv(file_path)
    original_headers = df.columns
    new_headers = [header.replace('exposure', 'outcome') for header in original_headers]
    df.columns = new_headers

    # 将修改后的文件保存
    df.to_csv(file_path, index=False)

    print(f"{file_path} 文件转为Outcome格式...")


def read_r_script(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return file.read()


def check_purified(gwas_id):
    data = MGwasData.query.get(gwas_id)
    file_path = f"./purified_data/{gwas_id}_purified.csv"
    if not data:
        print("PURIFICATION: Data not found")
        return False
    if not os.path.exists(file_path):
        print(f"PURIFICATION Error: File {file_path} not found")
        return False
    return True


def mark_data_status(gwas_id, data_status):
    data = MGwasData.query.get(gwas_id)
    data.state = data_status.value
    mr_db.session.add(data)
    mr_db.session.commit()


def mark_task_status(task_id, task_status):
    task = MTask.query.get(task_id)
    task.state = task_status.value
    mr_db.session.add(task)
    mr_db.session.commit()


def get_data_entity_name(gwas_id):
    data = MGwasData.query.get(gwas_id)
    return data.name


def delete_data_from_db(data_id):
    try:
        # 查找要删除的记录
        data = MGwasData.query.get(data_id)
        if not data:
            return jsonify({"message": "Data not found"}), 404

        # 删除记录
        mr_db.session.delete(data)
        mr_db.session.commit()
        return jsonify({"message": "Data deleted successfully"}), 200
    except Exception as e:
        mr_db.session.rollback()
        return jsonify({"message": "An error occurred", "error": str(e)}), 500


def delete_task_from_db(task_id):
    try:
        # 查找要删除的记录
        task = MTask.query.get(task_id)
        if not task:
            return jsonify({"message": "Task not found"}), 404

        # 删除记录
        mr_db.session.delete(task)
        mr_db.session.commit()
        return jsonify({"message": "Task deleted successfully"}), 200
    except Exception as e:
        mr_db.session.rollback()
        return jsonify({"message": "An error occurred", "error": str(e)}), 500


@app.route('/tasks', methods=['GET'])
def get_task_list():
    try:
        # 查询数据库获取所有任务
        task_list = MTask.query.all()

        # 将数据转换为字典列表
        data_list = [{
            'id': task.id,
            'name': task.name,
            'state': task.state,
            'type': task.type,
            'description': task.description,
            'create_time': task.create_time
        } for task in task_list]

        return jsonify({"data_list": data_list}), 200
    except Exception as e:
        print(f"Error retrieving tasks: {str(e)}")
        return jsonify({"message": f"Error retrieving tasks: {str(e)}"}), 500


def get_data_list():
    try:
        # 查询数据库获取所有 GWAS 数据
        gwas_data_list = MGwasData.query.all()

        # 将数据转换为字典列表
        data_list = [{
            'gwas_id': data.gwas_id,
            'name': data.name,
            'state': data.state
        } for data in gwas_data_list]
        return jsonify({"data_list": data_list}), 200
    except Exception as e:
        print(f"Error retrieving data: {str(e)}")
        return jsonify({"message": f"Error retrieving data: {str(e)}"}), 500


def upload_and_save_file(form, filename):
    try:
        file_info = {
            'gwas_id': form.get('gwas_id'),
            'name': form.get('name'),
            'state': DataStatus.RAW.value
        }
        gwas_data = save_gwas_data(file_info)
        logging_info = f"File {filename} is saved and information stored in database."
        print(logging_info)
        return jsonify({
            "message": logging_info,
            "data": {
                'gwas_id': gwas_data.gwas_id,
                'name': gwas_data.name,
                'state': DataStatus.RAW.value
            }
        }), 201
    except GWASDataExistsError as e:
        return jsonify({"message": str(e)}), 409
    except Exception as e:
        print(f"File info save error: {str(e)}")
        return jsonify({"message": f"File info save error: {str(e)}"}), 500


def save_gwas_data(file_info):
    try:
        gwas_data = MGwasData.query.filter_by(gwas_id=file_info.get('gwas_id')).first()

        if gwas_data:
            gwas_data.state = file_info.get('state', 'PURIFYING')
        else:
            gwas_data = MGwasData(
                gwas_id=file_info.get('gwas_id'),
                name=file_info.get('name'),
                state=file_info.get('state', 'PURIFYING')
            )
            mr_db.session.add(gwas_data)

        mr_db.session.commit()
        print("File information saved/updated successfully.")
        return gwas_data
    except SQLAlchemyError as e:
        mr_db.session.rollback()
        print(f"Database error occurred: {str(e)}")
        raise


def init_preset_data():
    file_info = {
        'gwas_id': 'ieu-a-2',
        'name': 'BMI',
        'state': 'MISSING'
    }

    save_gwas_data(file_info)

    file_info = {
        'gwas_id': 'ieu-a-7',
        'name': '冠心病',
        'state': 'MISSING'
    }

    save_gwas_data(file_info)


def setup():
    # 创建所有表
    mr_db.create_all()

    # with mr_db.engine.connect() as conn:
    #     conn.execute(text('ALTER TABLE m_gwas_data DROP CONSTRAINT IF EXISTS m_gwas_data_pkey'))
    #     conn.execute(text('ALTER TABLE m_gwas_data DROP COLUMN IF EXISTS id'))
    #     conn.execute(text('ALTER TABLE m_gwas_data DROP COLUMN IF EXISTS type'))
    #     conn.execute(text('ALTER TABLE m_gwas_data ADD PRIMARY KEY (gwas_id)'))

    # 初始化预设数据
    init_preset_data()

# setup()
