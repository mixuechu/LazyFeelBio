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


# 浏览数据
@app.route('/data', methods=['GET'])
def get_data():
    return get_data_list()


@app.route('/data/<string:data_id>', methods=['DELETE'])
def delete_data(data_id):
    return delete_data_from_db(data_id)

# 上传文件并保存到数据库的方法
@app.route('/data', methods=['POST'])
def upload_data():
    if 'file' not in request.files:
        return jsonify({"message": "No file part in the request"}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({"message": "No selected file"}), 400

    return upload_and_save_file(request.form, file)


# 浏览任务
@app.route('/tasks', methods=['GET'])
def get_tasks():
    return jsonify(tasks)


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