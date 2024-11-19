import React, { useState } from 'react';
import { Download, Upload, FileText, Play } from 'lucide-react';
import { saveAs } from 'file-saver';
import { useParticleStore } from '../store/useParticleStore';
import { useLaserStore } from '../store/useLaserStore';
import { useFieldStore } from '../store/useFieldStore';
import { useAlgorithmStore } from '../store/useAlgorithmStore';
import { useOutputStore } from '../store/useOutputStore';
import { useDataStore } from '../store/useDataStore';
import { parseLevelsFile, parseTransitionsFile } from '../utils/fileHandlers';
import {
  generateHeader,
  generateParticleSection,
  generateLaserSection,
  generateFieldSection,
  generateOutputSection,
  generateAlgorithmSection,
  generateScanSection,
} from '../utils/parameterFileGenerator';

interface FileConfig {
  name: string;
  description: string;
  path: string;
  type: 'input' | 'output';
}

const fileConfigs: FileConfig[] = [
  {
    name: 'Energy Levels',
    description: 'Contains energy level data',
    path: 'Data/n/n_levels.dat',
    type: 'input'
  },
  {
    name: 'Transition Lines',
    description: 'Contains transition line data',
    path: 'Data/n/n_lines.dat',
    type: 'input'
  },
  {
    name: 'Laser Spectrum',
    description: 'Laser spectrum data',
    path: 'Data/n/Laser_Spectrum.dat',
    type: 'input'
  },
  {
    name: 'Laser Intensity',
    description: 'Laser timing data',
    path: 'Data/n/Laser_Intensity.dat',
    type: 'input'
  },
  {
    name: 'Molecule Data',
    description: 'Output file for molecule data',
    path: 'Data/donnee_Mol.dat',
    type: 'output'
  },
  {
    name: 'Rate Data',
    description: 'Output file for rate data',
    path: 'Data/sortie_rate.dat',
    type: 'output'
  },
  {
    name: 'Data Card',
    description: 'Output file for data card',
    path: 'Data/data_card.dat',
    type: 'output'
  },
  {
    name: 'Random Generator State',
    description: 'Random number generator state',
    path: 'Data/random_gen.txt',
    type: 'output'
  }
];

export function FileManager() {
  const { particles, updateParticle } = useParticleStore();
  const { lasers, globalSettings: laserSettings, updateLaser, updateGlobalSettings } = useLaserStore();
  const { fields, updateFields } = useFieldStore();
  const { settings: algorithmSettings, updateSettings: updateAlgorithmSettings } = useAlgorithmStore();
  const { settings: outputSettings, updateSettings: updateOutputSettings } = useOutputStore();
  const { setLevels, setTransitions } = useDataStore();
  const [executablePath, setExecutablePath] = useState<string>('');

  const generateListeParamH = () => {
    let content = generateHeader();
    content += generateParticleSection(particles);
    content += generateLaserSection(lasers, laserSettings);
    content += generateFieldSection(fields);
    content += generateOutputSection(outputSettings);
    content += generateAlgorithmSection(algorithmSettings);
    content += generateScanSection({
      tauModif: '1e-3',
      isRandomScan: false,
      scanParams: []
    });
    return content;
  };

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>, fileConfig: FileConfig) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      
      try {
        switch (fileConfig.name) {
          case 'Energy Levels':
            const levels = parseLevelsFile(content);
            setLevels(levels);
            break;
          case 'Transition Lines':
            const transitions = parseTransitionsFile(content);
            setTransitions(transitions);
            break;
          case 'Liste_Param.h':
            // Parse parameter file and update stores
            break;
          default:
            console.log(`Handler for ${fileConfig.name} not implemented yet`);
        }
      } catch (error) {
        console.error(`Error parsing ${fileConfig.name}:`, error);
      }
    };
    reader.readAsText(file);
  };

  const handleExecutableSelect = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      setExecutablePath(file.path);
    }
  };

  const handleRunSimulation = () => {
    if (!executablePath) {
      alert('Please select the simulation executable first.');
      return;
    }
    // In a real application, you would use a backend API to:
    // 1. Save the parameter file
    // 2. Execute the simulation with the selected executable
    // 3. Monitor the simulation progress
    // 4. Load results when complete
    console.log('Running simulation with executable:', executablePath);
  };

  const handleDownload = (fileConfig: FileConfig) => {
    let content = '';
    switch (fileConfig.name) {
      case 'Liste_Param.h':
        content = generateListeParamH();
        break;
      default:
        content = '';
    }
    const blob = new Blob([content], { type: 'text/plain;charset=utf-8' });
    saveAs(blob, fileConfig.path.split('/').pop() || fileConfig.name);
  };

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center mb-6">
        <h2 className="text-lg font-semibold">File Management</h2>
      </div>

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
                onChange={handleExecutableSelect}
              />
            </div>
          </label>
          <button
            onClick={handleRunSimulation}
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

      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
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
                onChange={(e) => handleFileUpload(e, { name: 'Liste_Param.h', description: 'Parameter file', path: 'Liste_Param.h', type: 'input' })}
              />
            </label>
            <button
              onClick={() => handleDownload({ name: 'Liste_Param.h', description: 'Parameter file', path: 'Liste_Param.h', type: 'input' })}
              className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
            >
              <Download className="w-4 h-4 mr-2" />
              Download Liste_Param.h
            </button>
          </div>
        </div>

        <div className="space-y-4">
          <h3 className="text-md font-medium">Input Files</h3>
          <div className="space-y-2">
            {fileConfigs.filter(f => f.type === 'input').map((fileConfig) => (
              <div key={fileConfig.path} className="flex items-center justify-between bg-white p-3 rounded-lg shadow-sm">
                <div className="flex items-center">
                  <FileText className="w-5 h-5 text-gray-400 mr-2" />
                  <div>
                    <p className="text-sm font-medium">{fileConfig.name}</p>
                    <p className="text-xs text-gray-500">{fileConfig.path}</p>
                  </div>
                </div>
                <label className="flex items-center px-3 py-1 bg-blue-100 text-blue-700 rounded hover:bg-blue-200 cursor-pointer">
                  <Upload className="w-4 h-4 mr-1" />
                  Upload
                  <input
                    type="file"
                    className="hidden"
                    accept=".dat,.txt"
                    onChange={(e) => handleFileUpload(e, fileConfig)}
                  />
                </label>
              </div>
            ))}
          </div>
        </div>

        <div className="space-y-4">
          <h3 className="text-md font-medium">Output Files</h3>
          <div className="space-y-2">
            {fileConfigs.filter(f => f.type === 'output').map((fileConfig) => (
              <div key={fileConfig.path} className="flex items-center justify-between bg-white p-3 rounded-lg shadow-sm">
                <div className="flex items-center">
                  <FileText className="w-5 h-5 text-gray-400 mr-2" />
                  <div>
                    <p className="text-sm font-medium">{fileConfig.name}</p>
                    <p className="text-xs text-gray-500">{fileConfig.path}</p>
                  </div>
                </div>
                <button
                  onClick={() => handleDownload(fileConfig)}
                  className="flex items-center px-3 py-1 bg-blue-100 text-blue-700 rounded hover:bg-blue-200"
                >
                  <Download className="w-4 h-4 mr-1" />
                  Download
                </button>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}