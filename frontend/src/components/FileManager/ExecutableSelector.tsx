import React from 'react';
import { Upload, Play } from 'lucide-react';

interface ExecutableSelectorProps {
  executablePath: string;
  onExecutableSelect: (event: React.ChangeEvent<HTMLInputElement>) => void;
  onRunSimulation: () => void;
}

export function ExecutableSelector({ executablePath, onExecutableSelect, onRunSimulation }: ExecutableSelectorProps) {
  return (
    <div className="bg-white p-4 rounded-lg shadow-sm mb-6">
      <h3 className="text-md font-medium mb-4">Simulation Executable</h3>
      <div className="flex items-center gap-4">
        <label className="flex-1">
          <div className="flex items-center px-4 py-2 bg-gray-100 text-gray-700 rounded-md hover:bg-gray-200 cursor-pointer">
            <Upload className="w-4 h-4 mr-2" />
            {executablePath ? 'Selected: ' + executablePath : 'Select Executable'}
            <input
              type="file"
              className="hidden"
              accept=".exe"
              onChange={onExecutableSelect}
            />
          </div>
        </label>
        <button
          onClick={onRunSimulation}
          disabled={!executablePath}
          className={`flex items-center px-4 py-2 rounded-md ${
            executablePath 
              ? 'bg-green-600 text-white hover:bg-green-700' 
              : 'bg-gray-300 text-gray-500 cursor-not-allowed'
          }`}
        >
          <Play className="w-4 h-4 mr-2" />
          Run Simulation
        </button>
      </div>
    </div>
  );
}