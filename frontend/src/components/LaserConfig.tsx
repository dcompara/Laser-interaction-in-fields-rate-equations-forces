import React from 'react';
import { useLaserStore } from '../store/useLaserStore';
import { ParameterInput } from './common/ParameterInput';
import { Plus, Trash2 } from 'lucide-react';
import { parameterDescriptions } from '../utils/parameterDescriptions';

export function LaserConfig() {
  const { lasers, addLaser, removeLaser, updateLaser, globalSettings, updateGlobalSettings } = useLaserStore();

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Laser Configuration</h2>
        <button
          onClick={addLaser}
          className="px-3 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 flex items-center"
        >
          <Plus className="w-4 h-4 mr-2" />
          Add Laser
        </button>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Global Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Power Scale"
            name="scale_Power"
            value={globalSettings.powerScale}
            type="number"
            min={0}
            step="0.1"
            onChange={(value) => updateGlobalSettings({ powerScale: Number(value) })}
            tooltip={parameterDescriptions['scale_Power']}
          />
          <ParameterInput
            label="Detuning Offset"
            name="Offset_Detuning_cm"
            value={globalSettings.detuningOffset}
            type="number"
            step="0.01"
            unit="cm⁻¹"
            onChange={(value) => updateGlobalSettings({ detuningOffset: Number(value) })}
            tooltip={parameterDescriptions['Offset_Detuning_cm']}
          />
          <ParameterInput
            label="Linewidth Scale"
            name="scale_Gamma"
            value={globalSettings.linewidthScale}
            type="number"
            min={0}
            step="0.1"
            onChange={(value) => updateGlobalSettings({ linewidthScale: Number(value) })}
            tooltip={parameterDescriptions['scale_Gamma']}
          />
          <ParameterInput
            label="Force Rates"
            name="is_forced_rates"
            value={globalSettings.forcedRates}
            type="boolean"
            onChange={(value) => updateGlobalSettings({ forcedRates: Boolean(value) })}
            tooltip={parameterDescriptions['is_forced_rates']}
          />
        </div>
      </div>

      {lasers.map((laser, index) => (
        <div key={index} className="bg-white p-4 rounded-lg shadow-sm space-y-4">
          <div className="flex justify-between items-center">
            <h3 className="font-medium">Laser {index + 1}</h3>
            <button
              onClick={() => removeLaser(index)}
              className="text-red-600 hover:text-red-700"
            >
              <Trash2 className="w-4 h-4" />
            </button>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Position X"
              name={`waist_pos_x[${index}]`}
              value={laser.position.x}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateLaser(index, {
                position: { ...laser.position, x: Number(value) }
              })}
              tooltip={parameterDescriptions['waist_pos']}
            />
            <ParameterInput
              label="Position Y"
              name={`waist_pos_y[${index}]`}
              value={laser.position.y}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateLaser(index, {
                position: { ...laser.position, y: Number(value) }
              })}
              tooltip={parameterDescriptions['waist_pos']}
            />
            <ParameterInput
              label="Position Z"
              name={`waist_pos_z[${index}]`}
              value={laser.position.z}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateLaser(index, {
                position: { ...laser.position, z: Number(value) }
              })}
              tooltip={parameterDescriptions['waist_pos']}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Direction X"
              name={`direction_x[${index}]`}
              value={laser.direction.x}
              type="number"
              step="0.1"
              min={-1}
              max={1}
              onChange={(value) => updateLaser(index, {
                direction: { ...laser.direction, x: Number(value) }
              })}
              tooltip={parameterDescriptions['direction']}
            />
            <ParameterInput
              label="Direction Y"
              name={`direction_y[${index}]`}
              value={laser.direction.y}
              type="number"
              step="0.1"
              min={-1}
              max={1}
              onChange={(value) => updateLaser(index, {
                direction: { ...laser.direction, y: Number(value) }
              })}
              tooltip={parameterDescriptions['direction']}
            />
            <ParameterInput
              label="Direction Z"
              name={`direction_z[${index}]`}
              value={laser.direction.z}
              type="number"
              step="0.1"
              min={-1}
              max={1}
              onChange={(value) => updateLaser(index, {
                direction: { ...laser.direction, z: Number(value) }
              })}
              tooltip={parameterDescriptions['direction']}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <ParameterInput
              label="Waist"
              name={`waist[${index}]`}
              value={laser.waist}
              type="number"
              min={0}
              step="0.000001"
              unit="m"
              onChange={(value) => updateLaser(index, { waist: Number(value) })}
              tooltip={parameterDescriptions['waist']}
            />
            <ParameterInput
              label="Energy"
              name={`Energie_cm[${index}]`}
              value={laser.energy}
              type="number"
              step="0.000001"
              unit="cm⁻¹"
              onChange={(value) => updateLaser(index, { energy: Number(value) })}
              tooltip={parameterDescriptions['Energie_cm']}
            />
            <ParameterInput
              label="Linewidth"
              name={`Gamma_L_MHz[${index}]`}
              value={laser.linewidth}
              type="number"
              min={0}
              step="0.1"
              unit="MHz"
              onChange={(value) => updateLaser(index, { linewidth: Number(value) })}
              tooltip={parameterDescriptions['Gamma_L_MHz']}
            />
            <ParameterInput
              label="Power"
              name={`Power[${index}]`}
              value={laser.power}
              type="number"
              min={0}
              step="0.000001"
              unit="W"
              onChange={(value) => updateLaser(index, { power: Number(value) })}
              tooltip={parameterDescriptions['Power']}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Left Circular Polarization"
              name={`Pol_circulaire_left_sp[${index}]`}
              value={laser.polarization.leftCircular}
              type="number"
              step="0.000001"
              onChange={(value) => updateLaser(index, {
                polarization: { ...laser.polarization, leftCircular: Number(value) }
              })}
            />
            <ParameterInput
              label="Right Circular Polarization"
              name={`Pol_circulaire_right_sm[${index}]`}
              value={laser.polarization.rightCircular}
              type="number"
              step="0.000001"
              onChange={(value) => updateLaser(index, {
                polarization: { ...laser.polarization, rightCircular: Number(value) }
              })}
            />
            <ParameterInput
              label="Polarization Angle"
              name={`polar_angle_degree[${index}]`}
              value={laser.polarAngle}
              type="number"
              min={0}
              max={360}
              step="0.1"
              unit="°"
              onChange={(value) => updateLaser(index, { polarAngle: Number(value) })}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <ParameterInput
              label="Laser Type"
              name={`type_laser[${index}]`}
              value={laser.type}
              type="select"
              options={[
                'CW',
                'Femtosecond',
                'Multi-mode',
                'Pulse',
                'Gaussian',
                'Lorentzian',
                'Frequency Comb'
              ]}
              onChange={(value) => {
                const typeMap: { [key: string]: string } = {
                  'CW': '0',
                  'Femtosecond': '1',
                  'Multi-mode': '2',
                  'Pulse': '3',
                  'Gaussian': '5',
                  'Lorentzian': '6',
                  'Frequency Comb': '7'
                };
                updateLaser(index, { type: typeMap[value] });
              }}
              tooltip={parameterDescriptions['type_laser']}
            />
            <ParameterInput
              label="Coherent with Laser"
              name={`coherent_avec_laser_num[${index}]`}
              value={laser.coherentLaser}
              type="number"
              min={-1}
              step={1}
              onChange={(value) => updateLaser(index, { coherentLaser: Number(value) })}
              tooltip={parameterDescriptions['coherent_avec_laser_num']}
            />
          </div>

          {laser.type === '1' && (
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
              <ParameterInput
                label="Offset Frequency"
                name={`nu_offset_MHz[${index}]`}
                value={laser.femtosecond?.offsetFreq || 0}
                type="number"
                step="0.1"
                unit="MHz"
                onChange={(value) => updateLaser(index, {
                  femtosecond: { ...laser.femtosecond, offsetFreq: Number(value) }
                })}
              />
              <ParameterInput
                label="Repetition Rate"
                name={`nu_repetition_MHz[${index}]`}
                value={laser.femtosecond?.repRate || 80}
                type="number"
                step="0.1"
                unit="MHz"
                onChange={(value) => updateLaser(index, {
                  femtosecond: { ...laser.femtosecond, repRate: Number(value) }
                })}
              />
              <ParameterInput
                label="Comb Line Spacing"
                name={`nu_individual_comb_line_MHz[${index}]`}
                value={laser.femtosecond?.combLineSpacing || 80}
                type="number"
                step="0.1"
                unit="MHz"
                onChange={(value) => updateLaser(index, {
                  femtosecond: { ...laser.femtosecond, combLineSpacing: Number(value) }
                })}
              />
            </div>
          )}
        </div>
      ))}
    </div>
  );
}