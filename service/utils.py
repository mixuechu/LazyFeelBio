import os
from flask import jsonify
from sqlalchemy.exc import SQLAlchemyError
from mr_db import *
from sqlalchemy import text


class GWASDataExistsError(Exception):
    pass


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


def get_data_list():
    try:
        # 查询数据库获取所有 GWAS 数据
        gwas_data_list = MGwasData.query.all()

        # 将数据转换为字典列表
        data_list = [{
            # 'id': data.id,
            # 'type': data.type,
            'gwas_id': data.gwas_id,
            'name': data.name,
            'is_downloaded': data.is_downloaded
        } for data in gwas_data_list]
        return jsonify({"data_list": data_list}), 200
    except Exception as e:
        print(f"Error retrieving data: {str(e)}")
        return jsonify({"message": f"Error retrieving data: {str(e)}"}), 500


def upload_and_save_file(form, file):
    try:
        filename = file.filename
        file_path = os.path.join('./gwas_data', filename)

        # 检查文件是否已经存在
        if os.path.exists(file_path):
            print(f"Warning: File {filename} already exists and will be overwritten")

        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        file.save(file_path)
        print(f"File {filename} saved to path {file_path}")
    except Exception as e:
        print(f"File save error: {str(e)}")
        return jsonify({"message": f"File save error: {str(e)}"}), 400

    try:
        file_info = {
            # 'type': form.get('type'),
            'gwas_id': form.get('gwas_id'),
            'name': form.get('name'),
            'is_downloaded': True
        }
        gwas_data = save_gwas_data(file_info)
        logging_info = f"File {filename} is saved and information stored in database."
        print(logging_info)
        return jsonify({
            "message": logging_info,
            "data": {
                # 'id': gwas_data.id,
                # 'type': gwas_data.type,
                'gwas_id': gwas_data.gwas_id,
                'name': gwas_data.name,
                'is_downloaded': gwas_data.is_downloaded,
            }
        }), 201
    except GWASDataExistsError as e:
        return jsonify({"message": str(e)}), 409
    except Exception as e:
        print(f"File info save error: {str(e)}")
        return jsonify({"message": f"File info save error: {str(e)}"}), 500


def save_gwas_data(file_info):
    try:
        # 检查 gwas_id 是否已经存在
        existing_gwas = MGwasData.query.filter_by(gwas_id=file_info.get('gwas_id')).first()
        if existing_gwas:
            raise GWASDataExistsError(f"GWAS data with gwas_id {file_info.get('gwas_id')} already exists.")

        # 如果不存在，则创建新的记录
        gwas_data = MGwasData(
            # type=file_info.get('type'),
            gwas_id=file_info.get('gwas_id'),
            name=file_info.get('name'),
            is_downloaded=file_info.get('is_downloaded')
        )
        mr_db.session.add(gwas_data)
        mr_db.session.commit()
        print("File information saved successfully.")
        return gwas_data
    except GWASDataExistsError as e:
        print(str(e))
        raise
    except SQLAlchemyError as e:
        mr_db.session.rollback()
        print(f"Database error occurred: {str(e)}")
        raise


def init_preset_data():
    file_info = {
        # 'type': 'GWAS',
        'gwas_id': 'ieu-a-2',
        'name': 'BMI',
        'is_downloaded': True
    }

    save_gwas_data(file_info)

    file_info = {
        # 'type': 'GWAS',
        'gwas_id': 'ieu-a-7',
        'name': '冠心病',
        'is_downloaded': True
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


setup()
