import React, {useState, useEffect} from 'react';
import {fetchData, uploadFile, deleteData} from '../services/api';
import './BrowseData.css'; // 导入样式文件

function BrowseData() {
    const [data, setData] = useState([]);
    const [file, setFile] = useState(null);
    const [formData, setFormData] = useState({
        gwas_id: '', name: ''
    });
    const [uploadMessage, setUploadMessage] = useState('');
    const [error, setError] = useState('');

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
            ...prevState, [name]: value
        }));
    };

    const handleSubmit = async event => {
        event.preventDefault();

        if (!file) {
            setUploadMessage('请选择要上传的文件！');
            return;
        }

        const formDataObject = new FormData();
        formDataObject.append('file', file);
        formDataObject.append('gwas_id', formData.gwas_id);
        formDataObject.append('name', formData.name);

        try {
            const response = await uploadFile(formDataObject);
            setUploadMessage(response.data.message);
            // 重新获取数据列表
            fetchData()
                .then(response => {
                    setData(response.data.data_list);
                })
                .catch(error => {
                    console.error('获取数据列表时出错：', error);
                });
        } catch (error) {
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

    return (<div>
        <h1>查看数据</h1>
        <table>
            <thead>
            <tr>
                <th>GWAS ID</th>
                <th>名称</th>
                <th>是否已下载</th>
                <th>操作</th>
            </tr>
            </thead>
            <tbody>
            {data.map(item => (<tr key={item.gwas_id}>
                <td>{item.gwas_id}</td>
                <td>{item.name}</td>
                <td>{item.is_downloaded ? '是' : '否'}</td>
                <td>
                    <button onClick={() => handleDelete(item.gwas_id)}>删除</button>
                </td>
            </tr>))}
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
    </div>);
}

export default BrowseData;
