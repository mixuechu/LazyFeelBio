// api.js

import axios from 'axios';

const api = axios.create({
    baseURL: 'http://192.168.31.158:7001', // 替换为你的后端URL
});

export const fetchData = () => api.get('/data');
export const uploadFile = (formData) => api.post('/data', formData, {
    headers: {
        'Content-Type': 'multipart/form-data',
    },
});
export const deleteData = (id) => api.delete(`/data/${id}`);
export const deleteTask = (taskId) => api.delete(`/task/${taskId}`);
export const purifyData = (data_id, entity_name) => api.post(`/data/${data_id}`, {entity_name});
export const fetchTasks = () => api.get('/tasks');
export const createTask = (task) => api.post('/tasks', task);
export const fetchResults = () => api.get('/results');


export const uploadFileChunked = async (file, formData, onUploadProgress) => {
    const CHUNK_SIZE = 5 * 1024 * 1024; // 5MB
    const totalChunks = Math.ceil(file.size / CHUNK_SIZE);
    let uploadedChunks = 0;

    for (let start = 0; start < file.size; start += CHUNK_SIZE) {
        const chunk = file.slice(start, start + CHUNK_SIZE);
        const chunkFormData = new FormData();
        const chunkIndex = Math.floor(start / CHUNK_SIZE) + 1;
        chunkFormData.append('file', chunk, `${file.name}_part${chunkIndex}`);
        chunkFormData.append('gwas_id', formData.gwas_id);
        chunkFormData.append('name', formData.name);
        chunkFormData.append('chunkNumber', chunkIndex);
        chunkFormData.append('totalChunks', totalChunks);

        await api.post('/upload-chunk', chunkFormData, {
            onUploadProgress: (progressEvent) => {
                uploadedChunks++;
                const percentCompleted = Math.round((uploadedChunks / totalChunks) * 100);
                onUploadProgress(percentCompleted);
            }
        });
    }
};