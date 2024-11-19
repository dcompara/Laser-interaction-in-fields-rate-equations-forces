import React from 'react';
import { useDataStore } from '../store/useDataStore';
import { Download, Filter } from 'lucide-react';

export function Controls() {
  const { levels, transitions } = useDataStore();

  const handleSave = (type: 'levels' | 'transitions') => {
    const data = type === 'levels' ? levels : transitions;
    const content = type === 'levels'
      ? levels.map(level => 
          `${level.manifold} ${level.magneticNumber} ${level.boundLevel} ${level.stateNumber} ` +
          `${level.population} ${level.vibrational || 0} ${level.angularMomentum || 0} ${level.energy} ` +
          `${level.electricShift} ${level.magneticShift}`
        ).join('\n')
      : transitions.map(transition =>
          `${transition.upperState.manifold} ${transition.upperState.magneticNumber} ` +
          `${transition.upperState.boundLevel} ${transition.upperState.stateNumber} ` +
          `${transition.lowerState.manifold} ${transition.lowerState.magneticNumber} ` +
          `${transition.lowerState.boundLevel} ${transition.lowerState.stateNumber} ` +
          `${transition.upperEnergy} ${transition.lowerEnergy} ${transition.dipoleStrength}`
        ).join('\n');

    const blob = new Blob([content], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${type}.txt`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  return (
    <div className="flex gap-4 p-4 bg-white rounded-lg shadow-md">
      <button
        onClick={() => handleSave('levels')}
        className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
      >
        <Download className="w-4 h-4 mr-2" />
        Save Levels
      </button>
      <button
        onClick={() => handleSave('transitions')}
        className="flex items-center px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
      >
        <Download className="w-4 h-4 mr-2" />
        Save Transitions
      </button>
    </div>
  );
}