import React from 'react';
import Tippy from '@tippy.js/react';
import { HelpCircle } from 'lucide-react';
import 'tippy.js/dist/tippy.css';

interface ParameterInputProps {
  label: string;
  name: string;
  value: string | number;
  type: 'text' | 'number' | 'select' | 'boolean';
  onChange: (value: string | number | boolean) => void;
  options?: string[];
  min?: number;
  max?: number;
  step?: string;
  tooltip?: string;
  unit?: string;
}

export function ParameterInput({
  label,
  name,
  value,
  type,
  onChange,
  options = [],
  min,
  max,
  step,
  tooltip,
  unit
}: ParameterInputProps) {
  const id = `param-${name}`;

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
    if (type === 'boolean') {
      onChange((e.target as HTMLInputElement).checked);
    } else if (type === 'number') {
      onChange(Number(e.target.value));
    } else {
      onChange(e.target.value);
    }
  };

  return (
    <div className="space-y-1">
      <div className="flex items-center gap-2">
        <label htmlFor={id} className="block text-sm font-medium text-gray-700">
          {label}
          {unit && <span className="text-gray-500 text-xs ml-1">({unit})</span>}
        </label>
        {tooltip && (
          <Tippy 
            content={tooltip}
            placement="right"
            interactive={true}
            className="bg-white text-gray-900 shadow-lg rounded-lg p-2 text-sm max-w-md"
          >
            <span className="text-gray-400 hover:text-gray-600 cursor-help">
              <HelpCircle size={16} />
            </span>
          </Tippy>
        )}
      </div>
      {type === 'select' ? (
        <select
          id={id}
          value={value}
          onChange={handleChange}
          className="mt-1 block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-500 focus:ring-blue-500 sm:text-sm"
        >
          {options.map((option) => (
            <option key={option} value={option}>
              {option}
            </option>
          ))}
        </select>
      ) : type === 'boolean' ? (
        <input
          type="checkbox"
          id={id}
          checked={value as boolean}
          onChange={handleChange}
          className="mt-1 h-4 w-4 rounded border-gray-300 text-blue-600 focus:ring-blue-500"
        />
      ) : (
        <div className="relative mt-1 rounded-md shadow-sm">
          <input
            type={type}
            id={id}
            name={name}
            value={value}
            onChange={handleChange}
            min={min}
            max={max}
            step={step}
            className="block w-full rounded-md border-gray-300 shadow-sm focus:border-blue-500 focus:ring-blue-500 sm:text-sm"
          />
          {unit && (
            <div className="pointer-events-none absolute inset-y-0 right-0 flex items-center pr-3">
              <span className="text-gray-500 sm:text-sm">{unit}</span>
            </div>
          )}
        </div>
      )}
    </div>
  );
}