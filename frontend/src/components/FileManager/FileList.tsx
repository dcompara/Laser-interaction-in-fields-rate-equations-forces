import React from 'react';
import { Download, Upload, FileText } from 'lucide-react';
import { FileConfig } from './types';

interface FileListProps {
  title: string;
  files: FileConfig[];
  onUpload?: (event: React.ChangeEvent<HTMLInputElement>, file: FileConfig) => void;
  onDownload?: (file: FileConfig) => void;
}

export function FileList({ title, files, onUpload, onDownload }: FileListProps) {
  return (
    <div className="space-y-4">
      <h3 className="text-md font-medium">{title}</h3>
      <div className="space-y-2">
        {files.map((fileConfig) => (
          <div key={fileConfig.path} className="flex items-center justify-between bg-white p-3 rounded-lg shadow-sm">
            <div className="flex items-center">
              <FileText className="w-5 h-5 text-gray-400 mr-2" />
              <div>
                <p className="text-sm font-medium">{fileConfig.name}</p>
                <p className="text-xs text-gray-500">{fileConfig.path}</p>
              </div>
            </div>
            {onUpload && (
              <label className="flex items-center px-3 py-1 bg-blue-100 text-blue-700 rounded hover:bg-blue-200 cursor-pointer">
                <Upload className="w-4 h-4 mr-1" />
                Upload
                <input
                  type="file"
                  className="hidden"
                  accept=".dat,.txt"
                  onChange={(e) => onUpload(e, fileConfig)}
                />
              </label>
            )}
            {onDownload && (
              <button
                onClick={() => onDownload(fileConfig)}
                className="flex items-center px-3 py-1 bg-blue-100 text-blue-700 rounded hover:bg-blue-200"
              >
                <Download className="w-4 h-4 mr-1" />
                Download
              </button>
            )}
          </div>
        ))}
      </div>
    </div>
  );
}