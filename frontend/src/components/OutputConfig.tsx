import React from 'react';
import { useOutputStore } from '../store/useOutputStore';
import { ParameterInput } from './common/ParameterInput';

export function OutputConfig() {
  const { settings, updateSettings } = useOutputStore();

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Output Configuration</h2>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Graphics Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Display Size"
            name="SIZE_affichage"
            value={settings.graphics.displaySize}
            type="number"
            step="0.01"
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                displaySize: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Display Wait Time (s)"
            name="t_wait_affichage"
            value={settings.graphics.waitTime}
            type="number"
            step="0.1"
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                waitTime: Number(value)
              }
            })}
          />
        </div>

        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          <ParameterInput
            label="Rotation Axis X"
            name="rot_axe_x"
            value={settings.graphics.rotation.x}
            type="number"
            min={0}
            max={1}
            step={1}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                rotation: { ...settings.graphics.rotation, x: Number(value) }
              }
            })}
          />
          <ParameterInput
            label="Rotation Axis Y"
            name="rot_axe_y"
            value={settings.graphics.rotation.y}
            type="number"
            min={0}
            max={1}
            step={1}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                rotation: { ...settings.graphics.rotation, y: Number(value) }
              }
            })}
          />
          <ParameterInput
            label="Rotation Axis Z"
            name="rot_axe_z"
            value={settings.graphics.rotation.z}
            type="number"
            min={0}
            max={1}
            step={1}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                rotation: { ...settings.graphics.rotation, z: Number(value) }
              }
            })}
          />
        </div>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Output Files</h3>
        <div className="grid grid-cols-1 gap-4">
          <ParameterInput
            label="Data Card Output"
            name="is_DataCard_out"
            value={settings.files.dataCardOutput}
            type="select"
            options={['No Output', 'Single File', 'Split Files']}
            onChange={(value) => {
              const outputMap: { [key: string]: number } = {
                'No Output': 0,
                'Single File': 1,
                'Split Files': 2
              };
              updateSettings({
                ...settings,
                files: {
                  ...settings.files,
                  dataCardOutput: outputMap[value]
                }
              });
            }}
          />
        </div>
      </div>
    </div>
  );
}