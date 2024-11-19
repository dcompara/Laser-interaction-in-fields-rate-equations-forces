import React from 'react';
import { Download, Upload } from 'lucide-react';
import { saveAs } from 'file-saver';
import { useParticleStore } from '../store/useParticleStore';
import { useLaserStore } from '../store/useLaserStore';
import { useFieldStore } from '../store/useFieldStore';
import { useAlgorithmStore } from '../store/useAlgorithmStore';
import { useOutputStore } from '../store/useOutputStore';
import {
  generateHeader,
  generateParticleSection,
  generateLaserSection,
  generateFieldSection,
  generateOutputSection,
  generateAlgorithmSection,
  generateScanSection,
} from '../utils/parameterFileGenerator';
import { parseParameters } from '../utils/parameterFileParser';

export function FileManager() {
  const { particles, updateParticle } = useParticleStore();
  const { lasers, globalSettings: laserSettings, updateLaser, updateGlobalSettings } = useLaserStore();
  const { fields, updateFields } = useFieldStore();
  const { settings: algorithmSettings, updateSettings: updateAlgorithmSettings } = useAlgorithmStore();
  const { settings: outputSettings, updateSettings: updateOutputSettings } = useOutputStore();

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

  const handleUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      const {
        particles,
        lasers,
        fields,
        algorithm,
        output,
        laserGlobalSettings,
      } = parseParameters(content);

      // Update all stores with the parsed data
      particles.forEach((particle, index) => {
        updateParticle(index, particle);
      });

      lasers.forEach((laser, index) => {
        updateLaser(index, laser);
      });

      updateGlobalSettings(laserGlobalSettings);
      updateFields(fields);
      updateAlgorithmSettings(algorithm);
      updateOutputSettings(output);
    };
    reader.readAsText(file);
  };

  const handleDownload = () => {
    const content = generateListeParamH();
    const blob = new Blob([content], { type: 'text/plain;charset=utf-8' });
    saveAs(blob, 'Liste_Param.h');
  };

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Parameter File Management</h2>
      </div>

      <div className="flex gap-4">
        <label className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 cursor-pointer">
          <Upload className="w-4 h-4 mr-2" />
          Upload Liste_Param.h
          <input
            type="file"
            className="hidden"
            accept=".h"
            onChange={handleUpload}
          />
        </label>

        <button
          onClick={handleDownload}
          className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
        >
          <Download className="w-4 h-4 mr-2" />
          Download Liste_Param.h
        </button>
      </div>
    </div>
  );
}