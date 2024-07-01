import axios from 'axios';

const api = axios.create({
    baseURL: 'http://172.16.211.159:7001', // 替换为你的后端URL
});

export const fetchData = () => api.get('/data');
export const uploadFile = (formData) => api.post('/data', formData, {
    headers: {
        'Content-Type': 'multipart/form-data',
    },
});
export const deleteData = (id) => api.delete(`/data/${id}`);
export const fetchTasks = () => api.get('/tasks');
export const createTask = (task) => api.post('/tasks', task);
export const fetchResults = () => api.get('/results');
