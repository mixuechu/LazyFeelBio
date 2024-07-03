// BrowseData.js

import React, {useState, useEffect} from 'react';
import {fetchData, deleteData, purifyData, uploadFileChunked} from '../services/api';
import {CircularProgressbar, buildStyles} from 'react-circular-progressbar';
import 'react-circular-progressbar/dist/styles.css';
import './BrowseData.css'; // 导入样式文件

const DATA_STATUS_MAP = {
    "RAW": "原始",
    "PURIFIED": "已清洗",
    "MISSING": "缺失",
    "PURIFYING": "正在清洗",
    "FAILED": "清洗失败"
};

function BrowseData() {
    const [data, setData] = useState([]);
    const [file, setFile] = useState(null);
    const [formData, setFormData] = useState({
        gwas_id: '',
        name: ''
    });
    const [uploadMessage, setUploadMessage] = useState('');
    const [error, setError] = useState('');
    const [isUploading, setIsUploading] = useState(false);
    const [uploadProgress, setUploadProgress] = useState(0);

    useEffect(() => {
        fetchData()
            .then(response => {
                setData(response.data.data_list);
            })
            .catch(error => {
                console.error('获取数据列表时出错：', error);
            });
    }, []);

    const handleFileChange = event => {
        setFile(event.target.files[0]);
    };

    const handleChange = event => {
        const {name, value} = event.target;
        setFormData(prevState => ({
            ...prevState,
            [name]: value
        }));
    };

    const handleSubmit = async event => {
        event.preventDefault();

        if (!file) {
            setUploadMessage('请选择要上传的文件！');
            return;
        }

        setIsUploading(true);

        try {
            await uploadFileChunked(file, formData, (percentCompleted) => {
                setUploadProgress(percentCompleted);
            });
            setUploadMessage('文件上传成功');
            setIsUploading(false);
            setUploadProgress(0);
            // 重新获取数据列表
            fetchData()
                .then(response => {
                    setData(response.data.data_list);
                })
                .catch(error => {
                    console.error('获取数据列表时出错：', error);
                });
        } catch (error) {
            setIsUploading(false);
            setUploadProgress(0);
            if (error.response && error.response.data) {
                setUploadMessage(error.response.data.message);
            } else {
                console.error('上传文件时出错：', error);
                setUploadMessage('上传文件时出错，请稍后重试。');
            }
        }
    };

    const handleDelete = async (gwas_id) => {
        try {
            await deleteData(gwas_id);
            setData(prevData => prevData.filter(item => item.gwas_id !== gwas_id));
            setUploadMessage('数据删除成功');
        } catch (error) {
            console.error('删除数据时出错：', error);
            setError('删除数据时出错，请稍后重试。');
        }
    };

const handlePurify = async (gwas_id, entity_data) => {
    setData(prevData =>
        prevData.map(item =>
            item.gwas_id === gwas_id ? { ...item, state: 'PURIFYING' } : item
        )
    );
    try {
        await purifyData(gwas_id, entity_data);
        // 重新获取数据列表
        fetchData()
            .then(response => {
                setData(response.data.data_list);
            })
            .catch(error => {
                console.error('获取数据列表时出错：', error);
            });
    } catch (error) {
        console.error('净化数据时出错：', error);
        setError('净化数据时出错，请稍后重试。');
    }
};


    return (
        <div>
            <h1>查看数据</h1>
            {isUploading && (
                <div style={{width: '50px', height: '50px', margin: 'auto'}}>
                    <CircularProgressbar
                        value={uploadProgress}
                        text={`${uploadProgress}%`}
                        styles={buildStyles({
                            textSize: '20px',
                            pathColor: `rgba(62, 152, 199, ${uploadProgress / 100})`,
                            textColor: '#f88',
                            trailColor: '#d6d6d6',
                            backgroundColor: '#3e98c7',
                        })}
                    />
                </div>
            )}
            <table>
                <thead>
                <tr>
                    <th>GWAS ID</th>
                    <th>名称</th>
                    <th>状态</th>
                    <th>操作</th>
                </tr>
                </thead>
                <tbody>
                {data.map(item => (
                    <tr key={item.gwas_id}>
                        <td>{item.gwas_id}</td>
                        <td>{item.name}</td>
                        <td>{DATA_STATUS_MAP[item.state]}</td>
                        <td>
                            <button onClick={() => handleDelete(item.gwas_id)}>删除</button>
                            <button
                                onClick={() => handlePurify(item.gwas_id, item.name)}
                                disabled={item.state !== 'RAW' && item.state !== 'FAILED' && item.state !== 'PROCESSING'}
                            >
                                净化
                            </button>
                        </td>
                    </tr>
                ))}
                </tbody>
            </table>

            <h2>上传文件</h2>
            <form onSubmit={handleSubmit}>
                <div>
                    <label>GWAS ID:</label>
                    <input type="text" name="gwas_id" value={formData.gwas_id} onChange={handleChange}/>
                </div>
                <div>
                    <label>选择文件:</label>
                    <input type="file" onChange={handleFileChange}/>
                </div>
                <div>
                    <label>名称:</label>
                    <input type="text" name="name" value={formData.name} onChange={handleChange}/>
                </div>
                <button type="submit">上传</button>
            </form>
            {uploadMessage && <p>{uploadMessage}</p>}
            {error && <p>{error}</p>}
        </div>
    );
}

export default BrowseData;