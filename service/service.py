# service.py

import os
from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import logging
from enum import Enum

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
app = Flask(__name__)
app.config['JSON_AS_ASCII'] = False
CORS(app, supports_credentials=True)

from task_handler import create_task as demo_create_task
from utils import *

# In-memory storage for tasks and results
tasks = []
results = []

UPLOAD_FOLDER = 'gwas_data/'
CHUNKS_FOLDER = 'chunks/'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['CHUNKS_FOLDER'] = CHUNKS_FOLDER

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

if not os.path.exists(CHUNKS_FOLDER):
    os.makedirs(CHUNKS_FOLDER)


# 浏览数据
@app.route('/data', methods=['GET'])
def get_data():
    return get_data_list()


@app.route('/data/<string:data_id>', methods=['POST'])
def purify_data(data_id):
    try:
        entity_name = request.json.get('entity_name')
        if not entity_name:
            return jsonify({"error": "Missing entity_name"}), 400

        purification_req = {"data_id": data_id, "task_type": "DATA_PURIFICATION", "entity_name": entity_name}
        return demo_create_task(purification_req)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/data/<string:data_id>', methods=['DELETE'])
def delete_data(data_id):
    return delete_data_from_db(data_id)


@app.route('/upload-chunk', methods=['POST'])
def upload_chunk():
    if 'file' not in request.files:
        return jsonify({"message": "No file part in the request"}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({"message": "No selected file"}), 400

    chunk_number = int(request.form['chunkNumber'])
    total_chunks = int(request.form['totalChunks'])
    gwas_id = request.form['gwas_id']
    name = request.form['name']

    # 临时保存分块
    # chunk_filename = f"{file.filename}_part{chunk_number}"
    chunk_filepath = os.path.join(app.config['CHUNKS_FOLDER'], file.filename)
    file.save(chunk_filepath)

    original_filename = file.filename.split("_part")[0]

    # 检查是否所有分块都上传完毕
    # uploaded_chunks = len([f for f in os.listdir(app.config['CHUNKS_FOLDER']) if f.startswith(original_filename)])
    # print(uploaded_chunks)
    if chunk_number == total_chunks:
        # 合并所有分块
        complete_filename = original_filename
        complete_filepath = os.path.join(app.config['UPLOAD_FOLDER'], complete_filename)
        with open(complete_filepath, 'wb') as complete_file:
            for i in range(1, total_chunks + 1):
                chunk_filepath = os.path.join(app.config['CHUNKS_FOLDER'], file.filename).replace(
                    f"_part{total_chunks}", f"_part{i}")
                with open(chunk_filepath, 'rb') as chunk_file:
                    complete_file.write(chunk_file.read())
                os.remove(chunk_filepath.replace(f"_part{total_chunks}", f"_part{i}"))

        # 保存文件信息到数据库
        file_info = {
            'gwas_id': gwas_id,
            'name': name,
            'state': 'RAW'
        }
        return upload_and_save_file(file_info, complete_filename)

    return jsonify({"message": "Chunk uploaded successfully"}), 200


# 浏览任务
@app.route('/tasks', methods=['GET'])
def get_tasks():
    return get_task_list()


@app.route('/task/<string:task_id>', methods=['DELETE'])
def delete_task(task_id):
    return delete_task_from_db(task_id)


# 创建并任务
@app.route('/task', methods=['POST'])
def create_task():
    task_create_req = request.json
    task_id = demo_create_task(task_create_req)
    logging_info = f"Task {task_id} created and executed..."
    logging.info(logging_info)
    return jsonify({"message": logging_info}), 201


# 任务明细
@app.route('/task', methods=['GET'])
def get_results():
    return jsonify(results)


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=False, port=7001)