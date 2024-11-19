import React, { useState } from 'react';
import { DataTable } from './DataTable';
import { GrotrianDiagram } from './GrotrianDiagram';
import { Table, LineChart } from 'lucide-react';

export function TabPanel() {
  const [activeTab, setActiveTab] = useState<'levels' | 'transitions' | 'diagram'>('diagram');

  return (
    <div className="bg-white rounded-lg shadow-md">
      <div className="border-b border-gray-200">
        <nav className="flex -mb-px" aria-label="Tabs">
          <button
            onClick={() => setActiveTab('levels')}
            className={`${
              activeTab === 'levels'
                ? 'border-blue-500 text-blue-600'
                : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
            } flex items-center px-4 py-2 border-b-2 font-medium text-sm`}
          >
            <Table className="w-5 h-5 mr-2" />
            Levels
          </button>
          <button
            onClick={() => setActiveTab('transitions')}
            className={`${
              activeTab === 'transitions'
                ? 'border-blue-500 text-blue-600'
                : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
            } flex items-center px-4 py-2 border-b-2 font-medium text-sm`}
          >
            <Table className="w-5 h-5 mr-2" />
            Transitions
          </button>
          <button
            onClick={() => setActiveTab('diagram')}
            className={`${
              activeTab === 'diagram'
                ? 'border-blue-500 text-blue-600'
                : 'border-transparent text-gray-500 hover:text-gray-700 hover:border-gray-300'
            } flex items-center px-4 py-2 border-b-2 font-medium text-sm`}
          >
            <LineChart className="w-5 h-5 mr-2" />
            Grotrian Diagram
          </button>
        </nav>
      </div>
      <div className="p-4">
        {activeTab === 'levels' && <DataTable type="levels" />}
        {activeTab === 'transitions' && <DataTable type="transitions" />}
        {activeTab === 'diagram' && <GrotrianDiagram />}
      </div>
    </div>
  );
}