import React from 'react';
import { useFieldStore } from '../store/useFieldStore';
import { ParameterInput } from './common/ParameterInput';
import { parameterDescriptions } from '../utils/parameterDescriptions';

export function FieldConfig() {
  const { fields, updateFields, updateMagneticField, updateElectricField } = useFieldStore();

  const fieldTypes = [
    '0: Analytical (2nd order + nth order)',
    '1: Helmholtz/anti-Helmholtz coils',
    '2: Cylindrical 3D map (r,z,Fr,Fz)',
    '3: Cylindrical 3D map with derivatives',
    '4: 3D field map (x,y,z,Fx,Fy,Fz)'
  ];

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Field Configuration</h2>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Internal State Shift</h3>
        <ParameterInput
          label="Field Type for Internal State Shift"
          name="type_of_field_for_internal_state_shift"
          value={fields.typeForStateShift}
          type="select"
          options={[
            '0: Magnetic (Zeeman shift)',
            '1: Electric (Stark shift)'
          ]}
          onChange={(value) => updateFields({
            typeForStateShift: value === '0: Magnetic (Zeeman shift)' ? 0 : 1
          })}
          tooltip="Choose between magnetic (Zeeman) or electric (Stark) field for internal state shifts"
        />
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Field Definition Methods</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Electric Field Method"
            name="type_field_read_E"
            value={fields.electric.type}
            type="select"
            options={fieldTypes}
            onChange={(value) => updateElectricField({
              type: fieldTypes.indexOf(value)
            })}
            tooltip="Method for defining the electric field"
          />
          <ParameterInput
            label="Magnetic Field Method"
            name="type_field_read_B"
            value={fields.magnetic.type}
            type="select"
            options={fieldTypes}
            onChange={(value) => updateMagneticField({
              type: fieldTypes.indexOf(value)
            })}
            tooltip="Method for defining the magnetic field"
          />
        </div>
      </div>

      {fields.magnetic.type === 1 && (
        <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
          <h3 className="font-medium">Helmholtz Coil Configuration</h3>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <ParameterInput
              label="Number of Coils"
              name="Nb_bobines"
              value={fields.magnetic.helmholtz?.coilCount || 5}
              type="number"
              min={1}
              tooltip="Number of coils along the z-axis"
              onChange={(value) => updateMagneticField({
                helmholtz: {
                  ...fields.magnetic.helmholtz,
                  coilCount: Number(value)
                }
              })}
            />
            <ParameterInput
              label="Coil Configuration"
              name="is_Helmholtz"
              value={fields.magnetic.helmholtz?.isHelmholtz || 0}
              type="select"
              options={[
                '0: Single set of coils',
                '1: Helmholtz configuration',
                '-1: Anti-Helmholtz configuration'
              ]}
              onChange={(value) => updateMagneticField({
                helmholtz: {
                  ...fields.magnetic.helmholtz,
                  isHelmholtz: Number(value.split(':')[0])
                }
              })}
              tooltip="Configuration type for the coils"
            />
            <ParameterInput
              label="Coil Spacing"
              name="gap_bobines"
              value={fields.magnetic.helmholtz?.coilSpacing || 4e-2}
              type="number"
              step="0.001"
              unit="m"
              tooltip="Spacing between coils"
              onChange={(value) => updateMagneticField({
                helmholtz: {
                  ...fields.magnetic.helmholtz,
                  coilSpacing: Number(value)
                }
              })}
            />
            <ParameterInput
              label="Coil Current"
              name="courant_bobines"
              value={fields.magnetic.helmholtz?.current || 16000}
              type="number"
              unit="A"
              tooltip="Current in the coils"
              onChange={(value) => updateMagneticField({
                helmholtz: {
                  ...fields.magnetic.helmholtz,
                  current: Number(value)
                }
              })}
            />
            <ParameterInput
              label="Coil Radius"
              name="rayon_bobines"
              value={fields.magnetic.helmholtz?.radius || 1e-2}
              type="number"
              step="0.001"
              unit="m"
              tooltip="Radius of the coils"
              onChange={(value) => updateMagneticField({
                helmholtz: {
                  ...fields.magnetic.helmholtz,
                  radius: Number(value)
                }
              })}
            />
          </div>
        </div>
      )}

      {fields.magnetic.type === 0 && (
        <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
          <h3 className="font-medium">Magnetic Field Components</h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            {/* Field components */}
            <ParameterInput
              label="B_x"
              name="B_x"
              value={fields.magnetic.field.x}
              type="number"
              step="1e-10"
              unit="T"
              onChange={(value) => updateMagneticField({
                field: { ...fields.magnetic.field, x: Number(value) }
              })}
              tooltip="Magnetic field x component"
            />
            {/* Add similar inputs for y and z components */}
            {/* Add gradient and second gradient inputs */}
          </div>
        </div>
      )}

      {fields.electric.type === 0 && (
        <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
          <h3 className="font-medium">Electric Field Components</h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            {/* Similar structure to magnetic field components */}
          </div>
        </div>
      )}

      {(fields.magnetic.type >= 2 || fields.electric.type >= 2) && (
        <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
          <h3 className="font-medium">Grid Configuration</h3>
          <div className="grid grid-cols-1 gap-4">
            <ParameterInput
              label="Field Map File"
              name="field_map_file"
              value={fields.magnetic.gridConfig?.filePath || ''}
              type="text"
              tooltip="Path to the field map file"
              onChange={(value) => updateMagneticField({
                gridConfig: {
                  ...fields.magnetic.gridConfig,
                  filePath: String(value)
                }
              })}
            />
          </div>
        </div>
      )}
    </div>
  );
}