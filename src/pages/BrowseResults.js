import React, {useEffect, useState} from 'react';
import {fetchResults} from '../services/api';

function BrowseResults() {
    const [results, setResults] = useState([]); // 将初始状态设置为一个空数组


    useEffect(() => {
        fetchResults()
            .then(response => {
                setResults(response.data);
            })
            .catch(error => {
                console.error('获取结果列表时出错：', error);
            });
    }, []);


    return (
        <div>
            <h1>Browse Results</h1>
            {results && results.length > 0 ? ( // 检查 results 是否为 undefined 或 null
                <table>
                    <thead>
                    <tr>
                        <th>ID</th>
                        <th>Exposure ID</th>
                        <th>Outcome ID</th>
                        <th>Task ID</th>
                        <th>Finished</th>
                    </tr>
                    </thead>
                    <tbody>
                    {results.map(result => (
                        <tr key={result.id}>
                            <td>{result.id}</td>
                            <td>{result.exposure_id}</td>
                            <td>{result.outcome_id}</td>
                            <td>{result.task_id}</td>
                            <td>{result.finished ? 'Yes' : 'No'}</td>
                        </tr>
                    ))}
                    </tbody>
                </table>
            ) : (
                <p>No results available</p>
            )}
        </div>
    );
}

export default BrowseResults;