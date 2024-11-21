import React from 'react';
import { Upload } from 'lucide-react';

interface ExecutableSelectorProps {
  executablePath: string;
  onExecutableSelect: (event: React.ChangeEvent<HTMLInputElement>) => void;
}

export function ExecutableSelector({ executablePath, onExecutableSelect }: ExecutableSelectorProps) {
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
      </div>
      <p className="mt-2 text-sm text-gray-500">
        Select the simulation executable from your local computer. The simulation will be run locally when you click the "Run Simulation" button in the header.
      </p>
    </div>
  );
}