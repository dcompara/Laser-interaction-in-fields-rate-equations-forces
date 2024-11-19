import React from 'react';
import { Download, Upload } from 'lucide-react';
import { FileConfig } from './types';

interface ParameterFileProps {
  onUpload: (event: React.ChangeEvent<HTMLInputElement>, file: FileConfig) => void;
  onDownload: (file: FileConfig) => void;
}

export function ParameterFile({ onUpload, onDownload }: ParameterFileProps) {
  const paramFile: FileConfig = {
    name: 'Liste_Param.h',
    description: 'Parameter file',
    path: 'Liste_Param.h',
    type: 'input'
  };

  return (
    <div className="space-y-4">
      <h3 className="text-md font-medium">Parameter File</h3>
      <div className="flex gap-4">
        <label className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 cursor-pointer">
          <Upload className="w-4 h-4 mr-2" />
          Upload Liste_Param.h
          <input
            type="file"
            className="hidden"
            accept=".h"
            onChange={(e) => onUpload(e, paramFile)}
          />
        </label>
        <button
          onClick={() => onDownload(paramFile)}
          className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
        >
          <Download className="w-4 h-4 mr-2" />
          Download Liste_Param.h
        </button>
      </div>
    </div>
  );
}