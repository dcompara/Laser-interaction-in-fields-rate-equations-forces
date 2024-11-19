import React from 'react';
import { Plus, Trash2 } from 'lucide-react';
import { useScanStore, ScanParameter } from '../store/useScanStore';
import { ParameterInput } from './common/ParameterInput';

export function ScanConfig() {
  const {
    tauModif,
    isRandomScan,
    parameters,
    addParameter,
    removeParameter,
    updateParameter,
    updateTauModif,
    updateIsRandomScan
  } = useScanStore();

  const handleAddParameter = () => {
    addParameter({
      name: '',
      minValue: 0,
      maxValue: 1,
      steps: 10,
      isScanned: false,
      isTimeDependent: false
    });
  };

  return (
    <div className="space-y-6 p-4">
      <div className="flex justify-between items-center">
        <h2 className="text-lg font-semibold">Scan & Time-Varying Parameters</h2>
        <button
          onClick={handleAddParameter}
          className="px-3 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 flex items-center"
        >
          <Plus className="w-4 h-4 mr-2" />
          Add Parameter
        </button>
      </div>

      <div className="bg-white p-4 rounded-lg shadow-sm space-y-4">
        <h3 className="font-medium">Global Settings</h3>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <ParameterInput
            label="Default Time Constant (τ)"
            name="Tau_Modif"
            value={tauModif}
            type="number"
            step="0.0001"
            unit="s"
            tooltip="Default time constant for exponential modifications if not specified for individual parameters"
            onChange={(value) => updateTauModif(Number(value))}
          />
          <ParameterInput
            label="Random Scanning"
            name="is_Scan_Random"
            value={isRandomScan}
            type="boolean"
            tooltip="Enable random scanning of parameters instead of ordered scanning"
            onChange={(value) => updateIsRandomScan(Boolean(value))}
          />
        </div>
      </div>

      {parameters.map((param, index) => (
        <div key={index} className="bg-white p-4 rounded-lg shadow-sm space-y-4">
          <div className="flex justify-between items-center">
            <h3 className="font-medium">Parameter {index + 1}</h3>
            <button
              onClick={() => removeParameter(index)}
              className="text-red-600 hover:text-red-700"
            >
              <Trash2 className="w-4 h-4" />
            </button>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <ParameterInput
              label="Parameter Name"
              name={`param_name_${index}`}
              value={param.name}
              type="text"
              tooltip="Name of the parameter to scan or vary with time"
              onChange={(value) => updateParameter(index, { name: String(value) })}
            />
            <ParameterInput
              label="Number of Steps"
              name={`param_steps_${index}`}
              value={param.steps}
              type="number"
              min={1}
              tooltip="Number of intervals between min and max values"
              onChange={(value) => updateParameter(index, { steps: Number(value) })}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <ParameterInput
              label="Minimum Value"
              name={`param_min_${index}`}
              value={param.minValue}
              type="number"
              step="any"
              tooltip="Minimum value of the parameter range"
              onChange={(value) => updateParameter(index, { minValue: Number(value) })}
            />
            <ParameterInput
              label="Maximum Value"
              name={`param_max_${index}`}
              value={param.maxValue}
              type="number"
              step="any"
              tooltip="Maximum value of the parameter range"
              onChange={(value) => updateParameter(index, { maxValue: Number(value) })}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <ParameterInput
              label="Enable Scanning"
              name={`param_scan_${index}`}
              value={param.isScanned}
              type="boolean"
              tooltip="Enable parameter scanning within the specified range"
              onChange={(value) => updateParameter(index, { isScanned: Boolean(value) })}
            />
            <ParameterInput
              label="Time Dependent"
              name={`param_time_${index}`}
              value={param.isTimeDependent}
              type="boolean"
              tooltip="Enable time-dependent variation of the parameter"
              onChange={(value) => updateParameter(index, { isTimeDependent: Boolean(value) })}
            />
            {param.isTimeDependent && (
              <ParameterInput
                label="Time Constant (τ)"
                name={`param_tau_${index}`}
                value={param.tau || tauModif}
                type="number"
                step="0.0001"
                unit="s"
                tooltip="Time constant for exponential changes (if not specified, global τ is used)"
                onChange={(value) => updateParameter(index, { tau: Number(value) })}
              />
            )}
          </div>
        </div>
      ))}
    </div>
  );
}