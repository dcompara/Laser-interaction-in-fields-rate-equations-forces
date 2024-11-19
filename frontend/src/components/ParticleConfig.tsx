import React from 'react';
import { useParticleStore } from '../store/useParticleStore';
import { ParameterInput } from './common/ParameterInput';
import { Plus, Trash2 } from 'lucide-react';
import { parameterDescriptions } from '../utils/parameterDescriptions';

export function ParticleConfig() {
  const { particles, addParticle, removeParticle, updateParticle } = useParticleStore();

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Particle Configuration</h2>
        <button
          onClick={() => addParticle()}
          className="px-3 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 flex items-center"
        >
          <Plus className="w-4 h-4 mr-2" />
          Add Particle Type
        </button>
      </div>

      {particles.map((particle, index) => (
        <div key={index} className="bg-white p-4 rounded-lg shadow-sm space-y-4">
          <div className="flex justify-between items-center">
            <h3 className="font-medium">Particle Type {index + 1}</h3>
            {index > 0 && (
              <button
                onClick={() => removeParticle(index)}
                className="text-red-600 hover:text-red-700"
              >
                <Trash2 className="w-4 h-4" />
              </button>
            )}
          </div>

          <ParameterInput
            label="Particle Name"
            name={`Nom_Mol[${index}]`}
            value={particle.name}
            type="select"
            options={['BaF', 'Cs2', 'NH', 'Cs', 'CO', 'Li6Cs', 'Laminus', 'Li7Cs', 'Rb85Cs', 'Rb87Cs', 'Ps', 'C2minus', 'Ps_minus', 'P_bar', 'Osminus']}
            onChange={(value) => updateParticle(index, { name: value as string })}
            tooltip={parameterDescriptions['Nom_Mol']}
          />

          <ParameterInput
            label="Number of Particles"
            name={`N_Mol[${index}]`}
            value={particle.count}
            type="number"
            min={1}
            onChange={(value) => updateParticle(index, { count: Number(value) })}
            tooltip={parameterDescriptions['N_Mol']}
          />

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Initial Temperature X"
              name={`Temp_ini_x[${index}]`}
              value={particle.temperature.x}
              type="number"
              step="0.000001"
              unit="K"
              onChange={(value) => updateParticle(index, {
                temperature: { ...particle.temperature, x: Number(value) }
              })}
              tooltip={parameterDescriptions['Temp_ini']}
            />
            <ParameterInput
              label="Initial Temperature Y"
              name={`Temp_ini_y[${index}]`}
              value={particle.temperature.y}
              type="number"
              step="0.000001"
              unit="K"
              onChange={(value) => updateParticle(index, {
                temperature: { ...particle.temperature, y: Number(value) }
              })}
              tooltip={parameterDescriptions['Temp_ini']}
            />
            <ParameterInput
              label="Initial Temperature Z"
              name={`Temp_ini_z[${index}]`}
              value={particle.temperature.z}
              type="number"
              step="0.000001"
              unit="K"
              onChange={(value) => updateParticle(index, {
                temperature: { ...particle.temperature, z: Number(value) }
              })}
              tooltip={parameterDescriptions['Temp_ini']}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Initialization X"
              name={`Procedure_init_x[${index}]`}
              value={particle.initialization.x}
              type="select"
options={[
      '-1: Fixed size and ordered positions initially (one axis fixed, others randomized)', 
      '0: Fixed size determined by @size parameters (Gaussian distribution)', 
      '1: Magnetic linear potential (Laplace distribution, coefficient determined by F1)', 
      '2: Magnetic quadratic potential for neutral particles (Gaussian distribution)', 
      '3: Magnetic quadratic potential for charged particles (Gaussian distribution)', 
      '4: Quadratic electric potential (linear electric field for charged particles, e.g., in a Paul trap --> Gaussian distribution)', 
      '5: Perfectly ordered Gaussian in velocity (randomized positions)', 
      '6: Effusive beam (one component x, y, or z set to 5 to orient the beam)', 
      '7: Collimated effusive beam (one component x, y, or z set to 6 for orientation, others set to 0 for transverse temperature)'
    ]}
              onChange={(value) => updateParticle(index, {
                initialization: { ...particle.initialization, x: Number(value.split(':')[0]) }
              })}
              tooltip={parameterDescriptions['Procedure_init']}
            />
            <ParameterInput
              label="Initialization Y"
              name={`Procedure_init_y[${index}]`}
              value={particle.initialization.y}
              type="select"
options={[
      '-1: Fixed size and ordered positions initially (one axis fixed, others randomized)', 
      '0: Fixed size determined by @size parameters (Gaussian distribution)', 
      '1: Magnetic linear potential (Laplace distribution, coefficient determined by F1)', 
      '2: Magnetic quadratic potential for neutral particles (Gaussian distribution)', 
      '3: Magnetic quadratic potential for charged particles (Gaussian distribution)', 
      '4: Quadratic electric potential (linear electric field for charged particles, e.g., in a Paul trap --> Gaussian distribution)', 
      '5: Perfectly ordered Gaussian in velocity (randomized positions)', 
      '6: Effusive beam (one component x, y, or z set to 5 to orient the beam)', 
      '7: Collimated effusive beam (one component x, y, or z set to 6 for orientation, others set to 0 for transverse temperature)'
    ]}
              onChange={(value) => updateParticle(index, {
                initialization: { ...particle.initialization, y: Number(value.split(':')[0]) }
              })}
              tooltip={parameterDescriptions['Procedure_init']}
            />
            <ParameterInput
              label="Initialization Z"
              name={`Procedure_init_z[${index}]`}
              value={particle.initialization.z}
              type="select"
options={[
      '-1: Fixed size and ordered positions initially (one axis fixed, others randomized)', 
      '0: Fixed size determined by @size parameters (Gaussian distribution)', 
      '1: Magnetic linear potential (Laplace distribution, coefficient determined by F1)', 
      '2: Magnetic quadratic potential for neutral particles (Gaussian distribution)', 
      '3: Magnetic quadratic potential for charged particles (Gaussian distribution)', 
      '4: Quadratic electric potential (linear electric field for charged particles, e.g., in a Paul trap --> Gaussian distribution)', 
      '5: Perfectly ordered Gaussian in velocity (randomized positions)', 
      '6: Effusive beam (one component x, y, or z set to 5 to orient the beam)', 
      '7: Collimated effusive beam (one component x, y, or z set to 6 for orientation, others set to 0 for transverse temperature)'
    ]}
              onChange={(value) => updateParticle(index, {
                initialization: { ...particle.initialization, z: Number(value.split(':')[0]) }
              })}
              tooltip={parameterDescriptions['Procedure_init']}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Size X"
              name={`size_x[${index}]`}
              value={particle.size.x}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateParticle(index, {
                size: { ...particle.size, x: Number(value) }
              })}
              tooltip={parameterDescriptions['size']}
            />
            <ParameterInput
              label="Size Y"
              name={`size_y[${index}]`}
              value={particle.size.y}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateParticle(index, {
                size: { ...particle.size, y: Number(value) }
              })}
              tooltip={parameterDescriptions['size']}
            />
            <ParameterInput
              label="Size Z"
              name={`size_z[${index}]`}
              value={particle.size.z}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateParticle(index, {
                size: { ...particle.size, z: Number(value) }
              })}
              tooltip={parameterDescriptions['size']}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Position Offset X"
              name={`offset_x[${index}]`}
              value={particle.position.x}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateParticle(index, {
                position: { ...particle.position, x: Number(value) }
              })}
              tooltip={parameterDescriptions['offset']}
            />
            <ParameterInput
              label="Position Offset Y"
              name={`offset_y[${index}]`}
              value={particle.position.y}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateParticle(index, {
                position: { ...particle.position, y: Number(value) }
              })}
              tooltip={parameterDescriptions['offset']}
            />
            <ParameterInput
              label="Position Offset Z"
              name={`offset_z[${index}]`}
              value={particle.position.z}
              type="number"
              step="0.001"
              unit="m"
              onChange={(value) => updateParticle(index, {
                position: { ...particle.position, z: Number(value) }
              })}
              tooltip={parameterDescriptions['offset']}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Initial Velocity X"
              name={`v0_x[${index}]`}
              value={particle.velocity.x}
              type="number"
              step="0.1"
              unit="m/s"
              onChange={(value) => updateParticle(index, {
                velocity: { ...particle.velocity, x: Number(value) }
              })}
              tooltip={parameterDescriptions['v0']}
            />
            <ParameterInput
              label="Initial Velocity Y"
              name={`v0_y[${index}]`}
              value={particle.velocity.y}
              type="number"
              step="0.1"
              unit="m/s"
              onChange={(value) => updateParticle(index, {
                velocity: { ...particle.velocity, y: Number(value) }
              })}
              tooltip={parameterDescriptions['v0']}
            />
            <ParameterInput
              label="Initial Velocity Z"
              name={`v0_z[${index}]`}
              value={particle.velocity.z}
              type="number"
              step="0.1"
              unit="m/s"
              onChange={(value) => updateParticle(index, {
                velocity: { ...particle.velocity, z: Number(value) }
              })}
              tooltip={parameterDescriptions['v0']}
            />
          </div>
        </div>
      ))}
    </div>
  );
}