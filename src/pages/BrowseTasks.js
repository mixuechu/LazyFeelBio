import React, { useEffect, useState } from 'react';
import { deleteTask, fetchTasks, createTask, fetchData } from '../services/api';

const TASK_STATUSES_MAP = {
    "READY": "就绪",
    "PROCESSING": "在冲了，别急...",
    "FINISHED": "已完成",
    "FAILED": "失败",
};

const TASK_TYPE_MAP = {
    "BASIC_MR": "基础MR分析",
    "ZIP_DEMO": "DEMO 打包"
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
    const [dataList, setDataList] = useState([]);
    const [showDataList, setShowDataList] = useState(false);
    const [currentInput, setCurrentInput] = useState('');

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

    const handleInputChange = (e) => {
        const { name, value } = e.target;
        setNewTask(prevState => ({
            ...prevState,
            [name]: value
        }));
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

    const fetchDataList = async () => {
        try {
            const response = await fetchData();
            setDataList(response.data.data_list.filter(data => data.state === 'PURIFIED'));
        } catch (error) {
            console.error('获取数据列表时出错：', error);
            setError('获取数据列表时出错，请稍后重试。');
        }
    };

    const handleShowDataList = async (inputName) => {
        await fetchDataList();
        setCurrentInput(inputName);
        setShowDataList(true);
    };

    const handleSelectData = (id) => {
        setNewTask(prevState => ({
            ...prevState,
            [currentInput]: id
        }));
        setShowDataList(false);
    };

    const renderFormFields = () => {
        switch (newTask.type) {
            case 'BASIC_MR':
                return (
                    <>
                        <div>
                            <label>Exposure ID:</label>
                            <input
                                type="text"
                                name="exposure_id"
                                value={newTask.exposure_id || ''}
                                onClick={() => handleShowDataList('exposure_id')}
                                readOnly
                                required
                            />
                        </div>
                        <div>
                            <label>Outcome ID:</label>
                            <input
                                type="text"
                                name="outcome_id"
                                value={newTask.outcome_id || ''}
                                onClick={() => handleShowDataList('outcome_id')}
                                readOnly
                                required
                            />
                        </div>
                    </>
                );
            case 'ZIP_DEMO':
                return (
                    <>
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
                        <td>{TASK_STATUSES_MAP[task.state]}</td>
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

            {showDataList && (
                <div className="data-list-modal">
                    <h2>选择 {currentInput === 'exposure_id' ? 'Exposure ID' : 'Outcome ID'}</h2>
                    <ul>
                        {dataList.map(data => (
                            <li key={data.gwas_id} onClick={() => handleSelectData(data.gwas_id)}>
                                {data.gwas_id} - {data.name}
                            </li>
                        ))}
                    </ul>
                    <button onClick={() => setShowDataList(false)}>关闭</button>
                </div>
            )}
        </div>
    );
}

export default BrowseTasks;
