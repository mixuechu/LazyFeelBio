/* App.js */
import React, { useState } from 'react';
import BrowseData from './pages/BrowseData';
import BrowseTasks from './pages/BrowseTasks';
import CreateTask from './pages/CreateTask';
import BrowseResults from './pages/BrowseResults';
import './App.css'; // 导入主样式文件

function App() {
  const [activeTab, setActiveTab] = useState('browse-data');

  const renderTabContent = () => {
    switch (activeTab) {
      case 'browse-data':
        return <BrowseData />;
      case 'browse-tasks':
        return <BrowseTasks />;
      case 'create-task':
        return <CreateTask />;
      case 'browse-results':
        return <BrowseResults />;
      default:
        return null;
    }
  };

  return (
    <div className="App">
      <h1>功能管理</h1>
      <div className="tabs">
        <button className={activeTab === 'browse-data' ? 'active' : ''} onClick={() => setActiveTab('browse-data')}>
          查看数据
        </button>
        <button className={activeTab === 'browse-tasks' ? 'active' : ''} onClick={() => setActiveTab('browse-tasks')}>
          查看任务列表
        </button>
        <button className={activeTab === 'create-task' ? 'active' : ''} onClick={() => setActiveTab('create-task')}>
          创建任务
        </button>
        <button className={activeTab === 'browse-results' ? 'active' : ''} onClick={() => setActiveTab('browse-results')}>
          查看任务结果
        </button>
      </div>

      <div className="tab-content">
        {renderTabContent()}
      </div>
    </div>
  );
}

export default App;