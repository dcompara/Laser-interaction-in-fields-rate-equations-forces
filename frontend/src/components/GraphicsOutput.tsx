import React from 'react';
import { useOutputStore } from '../store/useOutputStore';
import { ParameterInput } from './common/ParameterInput';
import { parameterDescriptions } from '../utils/parameterDescriptions';

export function GraphicsOutput() {
  const { settings, updateSettings } = useOutputStore();

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Graphics & Output Configuration</h2>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Display Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Display Size"
            name="SIZE_affichage"
            value={settings.graphics.displaySize}
            type="number"
            step="0.01"
            unit="m"
            tooltip={parameterDescriptions['SIZE_affichage']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                displaySize: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Display Wait Time"
            name="t_wait_affichage"
            value={settings.graphics.waitTime}
            type="number"
            step="0.1"
            unit="s"
            tooltip={parameterDescriptions['t_wait_affichage']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                waitTime: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Enable Graphics"
            name="Graphics"
            value={settings.graphics.enabled}
            type="boolean"
            tooltip={parameterDescriptions['Graphics']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                enabled: Boolean(value)
              }
            })}
          />
        </div>

        <h3 className="font-medium mt-6">Rotation Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
          <ParameterInput
            label="Rotation X"
            name="rot_axe_x"
            value={settings.graphics.rotation.x}
            type="number"
            min={0}
            max={1}
            step={1}
            tooltip={parameterDescriptions['rot_axe']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                rotation: { ...settings.graphics.rotation, x: Number(value) }
              }
            })}
          />
          <ParameterInput
            label="Rotation Y"
            name="rot_axe_y"
            value={settings.graphics.rotation.y}
            type="number"
            min={0}
            max={1}
            step={1}
            tooltip={parameterDescriptions['rot_axe']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                rotation: { ...settings.graphics.rotation, y: Number(value) }
              }
            })}
          />
          <ParameterInput
            label="Rotation Z"
            name="rot_axe_z"
            value={settings.graphics.rotation.z}
            type="number"
            min={0}
            max={1}
            step={1}
            tooltip={parameterDescriptions['rot_axe']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                rotation: { ...settings.graphics.rotation, z: Number(value) }
              }
            })}
          />
        </div>

        <h3 className="font-medium mt-6">Molecule Drawing Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Max Vibration Level"
            name="v_max"
            value={settings.graphics.vMax}
            type="number"
            min={0}
            step={1}
            tooltip={parameterDescriptions['v_max']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                vMax: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Max Rotation Level (2J)"
            name="two_J_max"
            value={settings.graphics.twoJMax}
            type="number"
            min={0}
            step={2}
            tooltip={parameterDescriptions['two_J_max']}
            onChange={(value) => updateSettings({
              ...settings,
              graphics: {
                ...settings.graphics,
                twoJMax: Number(value)
              }
            })}
          />
        </div>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Output Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Manifold Exclusion"
            name="num_manifold_not_studied"
            value={settings.output.excludedManifold}
            type="number"
            tooltip={parameterDescriptions['num_manifold_not_studied']}
            onChange={(value) => updateSettings({
              ...settings,
              output: {
                ...settings.output,
                excludedManifold: Number(value)
              }
            })}
          />
          <ParameterInput
            label="Study Level"
            name="num_niveau_etudie"
            value={settings.output.studyLevel}
            type="number"
            tooltip={parameterDescriptions['num_niveau_etudie']}
            onChange={(value) => updateSettings({
              ...settings,
              output: {
                ...settings.output,
                studyLevel: Number(value)
              }
            })}
          />
        </div>
      </div>
    </div>
  );
}