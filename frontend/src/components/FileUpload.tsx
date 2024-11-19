import React from 'react';
import { Upload } from 'lucide-react';
import { useDataStore } from '../store/useDataStore';
import { EnergyLevel, Transition } from '../types';

export function FileUpload() {
  const { setLevels, setTransitions } = useDataStore();

  const parseLevelsFile = (content: string) => {
    const lines = content.trim().split('\n');
    const levels: EnergyLevel[] = lines.map(line => {
      const [
        manifold, magneticNumber, boundLevel, stateNumber,
        population, vibrational, angularMomentum, energy,
        electricShift, magneticShift
      ] = line.trim().split(/\s+/).map(Number);

      return {
        manifold,
        magneticNumber,
        boundLevel,
        stateNumber,
        population,
        vibrational,
        angularMomentum,
        energy,
        electricShift,
        magneticShift
      };
    });
    setLevels(levels);
  };

  const parseTransitionsFile = (content: string) => {
    const lines = content.trim().split('\n');
    const transitions: Transition[] = lines.map(line => {
      const [
        upperManifold, upperM, upperBound, upperState,
        lowerManifold, lowerM, lowerBound, lowerState,
        upperEnergy, lowerEnergy, dipoleStrength
      ] = line.trim().split(/\s+/).map(Number);

      return {
        upperState: {
          manifold: upperManifold,
          magneticNumber: upperM,
          boundLevel: upperBound,
          stateNumber: upperState
        },
        lowerState: {
          manifold: lowerManifold,
          magneticNumber: lowerM,
          boundLevel: lowerBound,
          stateNumber: lowerState
        },
        upperEnergy,
        lowerEnergy,
        dipoleStrength
      };
    });
    setTransitions(transitions);
  };

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>, type: 'levels' | 'transitions') => {
    const file = event.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target?.result as string;
      if (type === 'levels') {
        parseLevelsFile(content);
      } else {
        parseTransitionsFile(content);
      }
    };
    reader.readAsText(file);
  };

  return (
    <div className="flex gap-4 p-4 bg-white rounded-lg shadow-md">
      <div className="flex-1">
        <label className="flex flex-col items-center p-4 border-2 border-dashed border-gray-300 rounded-lg hover:border-blue-500 cursor-pointer">
          <Upload className="w-8 h-8 text-gray-400" />
          <span className="mt-2 text-sm text-gray-500">Upload Levels File</span>
          <input
            type="file"
            className="hidden"
            accept=".txt,.dat"
            onChange={(e) => handleFileUpload(e, 'levels')}
          />
        </label>
      </div>
      <div className="flex-1">
        <label className="flex flex-col items-center p-4 border-2 border-dashed border-gray-300 rounded-lg hover:border-blue-500 cursor-pointer">
          <Upload className="w-8 h-8 text-gray-400" />
          <span className="mt-2 text-sm text-gray-500">Upload Transitions File</span>
          <input
            type="file"
            className="hidden"
            accept=".txt,.dat"
            onChange={(e) => handleFileUpload(e, 'transitions')}
          />
        </label>
      </div>
    </div>
  );
}