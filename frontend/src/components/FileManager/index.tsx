import React, { useState } from 'react';
import { saveAs } from 'file-saver';
import { useParticleStore } from '../../store/useParticleStore';
import { useLaserStore } from '../../store/useLaserStore';
import { useFieldStore } from '../../store/useFieldStore';
import { useAlgorithmStore } from '../../store/useAlgorithmStore';
import { useOutputStore } from '../../store/useOutputStore';
import { useDataStore } from '../../store/useDataStore';
import { parseLevelsFile, parseTransitionsFile } from '../../utils/fileHandlers';
import { generateHeader, generateParticleSection, generateLaserSection, generateFieldSection, generateOutputSection, generateAlgorithmSection, generateScanSection } from '../../utils/parameterFileGenerator';
import { ExecutableSelector } from './ExecutableSelector';
import { ParameterFile } from './ParameterFile';
import { FileList } from './FileList';
import { FileConfig, fileConfigs } from './types';

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

  const inputFiles = fileConfigs.filter(f => f.type === 'input');
  const outputFiles = fileConfigs.filter(f => f.type === 'output');

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center mb-6">
        <h2 className="text-lg font-semibold">File Management</h2>
      </div>

      <ExecutableSelector
        executablePath={executablePath}
        onExecutableSelect={handleExecutableSelect}
        onRunSimulation={handleRunSimulation}
      />

      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <ParameterFile
          onUpload={handleFileUpload}
          onDownload={handleDownload}
        />

        <FileList
          title="Input Files"
          files={inputFiles}
          onUpload={handleFileUpload}
        />

        <FileList
          title="Output Files"
          files={outputFiles}
          onDownload={handleDownload}
        />
      </div>
    </div>
  );
}