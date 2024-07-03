// BrowseTasks.js

import React, { useEffect, useState } from 'react';
import { deleteTask, fetchTasks, createTask } from '../services/api';

const TASK_STATUES_MAP = {
    "READY": "就绪",
    "PROCESSING": "在冲了，别急...",
    "FINISHED": "已完成",
    "FAILED": "失败",
};

const TASK_TYPE_MAP = {
    "BASIC_MR": "基础MR分析",
    "ZIP_DEMO": "DEMO 打包",
    "DATA_PURIFICATION": "数据净化"
};

function BrowseTasks() {
    const [tasks, setTasks] = useState([]);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState(null);
    const [uploadMessage, setUploadMessage] = useState('');
    const [newTask, setNewTask] = useState({
        name: '',
        type: 'BASIC_MR',
        description: ''
    });

    useEffect(() => {
        const getTasks = async () => {
            try {
                const response = await fetchTasks();
                setTasks(response.data.data_list);
                setLoading(false);
            } catch (error) {
                setError('Error fetching tasks');
                setLoading(false);
            }
        };

        getTasks();
    }, []);

    if (loading) {
        return <div>Loading...</div>;
    }

    if (error) {
        return <div>{error}</div>;
    }

    const handleDelete = async (task_id) => {
        try {
            await deleteTask(task_id);
            setTasks(prevData => prevData.filter(item => item.id !== task_id));
            setUploadMessage('任务删除成功');
        } catch (error) {
            console.error('删除任务时出错：', error);
            setError('删除任务时出错，请稍后重试。');
        }
    };

    const handleCreateTask = async (e) => {
        e.preventDefault();
        try {
            const response = await createTask(newTask);
            setTasks(prevTasks => [...prevTasks, response.data]);
            setNewTask({ name: '', type: 'BASIC_MR', description: '' });
            setUploadMessage('任务创建成功');
        } catch (error) {
            console.error('创建任务时出错：', error);
            setError('创建任务时出错，请稍后重试。');
        }
    };

    const handleInputChange = (e) => {
        const { name, value } = e.target;
        setNewTask(prevState => ({
            ...prevState,
            [name]: value
        }));
    };

    const renderFormFields = () => {
        switch (newTask.type) {
            case 'DATA_PURIFICATION':
                return (
                    <>
                        <div>
                            <label>GWAS ID:</label>
                            <input
                                type="text"
                                name="gwas_id"
                                value={newTask.gwas_id || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                        <div>
                            <label>Task ID:</label>
                            <input
                                type="text"
                                name="task_id"
                                value={newTask.task_id || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                        <div>
                            <label>Entity Name:</label>
                            <input
                                type="text"
                                name="entity_name"
                                value={newTask.entity_name || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                    </>
                );
            case 'BASIC_MR':
                return (
                    <>
                        <div>
                            <label>Outcome Name:</label>
                            <input
                                type="text"
                                name="outcome_name"
                                value={newTask.outcome_name || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                        <div>
                            <label>Outcome File:</label>
                            <input
                                type="text"
                                name="outcome_file"
                                value={newTask.outcome_file || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                        <div>
                            <label>Exposure Name:</label>
                            <input
                                type="text"
                                name="exposure_name"
                                value={newTask.exposure_name || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                        <div>
                            <label>Exposure File:</label>
                            <input
                                type="text"
                                name="exposure_file"
                                value={newTask.exposure_file || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                        <div>
                            <label>ID:</label>
                            <input
                                type="text"
                                name="id"
                                value={newTask.id || ''}
                                onChange={handleInputChange}
                                required
                            />
                        </div>
                    </>
                );
            default:
                return null;
        }
    };

    return (
        <div>
            <h1>任务列表</h1>
            <form onSubmit={handleCreateTask}>
                <div>
                    <label>任务名称:</label>
                    <input
                        type="text"
                        name="name"
                        value={newTask.name}
                        onChange={handleInputChange}
                        required
                    />
                </div>
                <div>
                    <label>任务类型:</label>
                    <select
                        name="type"
                        value={newTask.type}
                        onChange={handleInputChange}
                    >
                        {Object.keys(TASK_TYPE_MAP).map(key => (
                            <option key={key} value={key}>
                                {TASK_TYPE_MAP[key]}
                            </option>
                        ))}
                    </select>
                </div>
                {renderFormFields()}
                <div>
                    <label>任务描述:</label>
                    <input
                        type="text"
                        name="description"
                        value={newTask.description}
                        onChange={handleInputChange}
                    />
                </div>
                <button type="submit">创建任务</button>
            </form>
            {uploadMessage && <p>{uploadMessage}</p>}
            <table>
                <thead>
                <tr>
                    <th>Name</th>
                    <th>State</th>
                    <th>Type</th>
                    <th>Description</th>
                    <th>Create Time</th>
                    <th>操作</th>
                </tr>
                </thead>
                <tbody>
                {tasks.map(task => (
                    <tr key={task.id}>
                        <td>{task.name}</td>
                        <td>{TASK_STATUES_MAP[task.state]}</td>
                        <td>{TASK_TYPE_MAP[task.type]}</td>
                        <td>{task.description}</td>
                        <td>{new Date(task.create_time).toLocaleString()}</td>
                        <td>
                            <button onClick={() => handleDelete(task.id)}>删除</button>
                        </td>
                    </tr>
                ))}
                </tbody>
            </table>
        </div>
    );
}

export default BrowseTasks;