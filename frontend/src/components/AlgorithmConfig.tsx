import React from 'react';
import { useAlgorithmStore } from '../store/useAlgorithmStore';
import { ParameterInput } from './common/ParameterInput';

export function AlgorithmConfig() {
  const { settings, updateSettings } = useAlgorithmStore();

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Algorithm Configuration</h2>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Monte Carlo Settings</h3>
        <ParameterInput
          label="Monte Carlo Algorithm"
          name="Choix_algorithme_Monte_Carlo"
          value={settings.monteCarlo.algorithm}
          type="select"
          options={[
            'None',
            'Kinetic Monte Carlo',
            'Random Selection Method',
            'First Reaction Method',
            'Fast Rough Method'
          ]}
          onChange={(value) => {
            const algorithmMap: { [key: string]: number } = {
              'None': -1,
              'Kinetic Monte Carlo': 0,
              'Random Selection Method': 1,
              'First Reaction Method': 2,
              'Fast Rough Method': 3
            };
            updateSettings({
              ...settings,
              monteCarlo: {
                ...settings.monteCarlo,
                algorithm: algorithmMap[value]
              }
            });
          }}
        />
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">N-Body Algorithm</h3>
        <ParameterInput
          label="N-Body Algorithm"
          name="Choix_algorithme_N_corps"
          value={settings.nBody.algorithm}
          type="select"
          options={[
            'None',
            'Verlet (Acceleration)',
            'Verlet (Potential)',
            'Yoshida6 (Acceleration)',
            'Yoshida6 (Potential)',
            'Verlet (High Order)',
            'Boris-Buneman'
          ]}
          onChange={(value) => {
            const algorithmMap: { [key: string]: number } = {
              'None': -1,
              'Verlet (Acceleration)': 1,
              'Verlet (Potential)': 2,
              'Yoshida6 (Acceleration)': 3,
              'Yoshida6 (Potential)': 4,
              'Verlet (High Order)': 6,
              'Boris-Buneman': 7
            };
            updateSettings({
              ...settings,
              nBody: {
                ...settings.nBody,
                algorithm: algorithmMap[value]
              }
            });
          }}
        />

        <ParameterInput
          label="Position Epsilon (m)"
          name="choix_epsilon"
          value={settings.nBody.epsilon}
          type="number"
          step="1e-9"
          onChange={(value) => updateSettings({
            ...settings,
            nBody: {
              ...settings.nBody,
              epsilon: Number(value)
            }
          })}
        />
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Time Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Dynamic Time Step"
            name="dt_dyn_epsilon_param"
            value={settings.time.dynamicTimeStep}
            type="number"
            step="1e-9"
            onChange={(value) => updateSettings({
              ...settings,
              time: {
                ...settings.time,
                dynamicTimeStep: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Final Time (s)"
            name="t_fin"
            value={settings.time.finalTime}
            type="number"
            step="1e-6"
            onChange={(value) => updateSettings({
              ...settings,
              time: {
                ...settings.time,
                finalTime: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Diagnostic Interval (s)"
            name="dt_dia"
            value={settings.time.diagnosticInterval}
            type="number"
            step="1e-6"
            onChange={(value) => updateSettings({
              ...settings,
              time: {
                ...settings.time,
                diagnosticInterval: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Output Interval (s)"
            name="dt_out"
            value={settings.time.outputInterval}
            type="number"
            step="1e-6"
            onChange={(value) => updateSettings({
              ...settings,
              time: {
                ...settings.time,
                outputInterval: Number(value)
              }
            })}
          />
        </div>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Random Number Generator</h3>
        <ParameterInput
          label="Random Seed"
          name="Seed_Init_Random_Number_Generator"
          value={settings.randomSeed}
          type="number"
          onChange={(value) => updateSettings({
            ...settings,
            randomSeed: Number(value)
          })}
        />
      </div>
    </div>
  );
}