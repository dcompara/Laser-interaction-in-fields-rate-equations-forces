export interface Parameter {
  name: string;
  value: string | number | boolean;
  description: string;
  type: 'number' | 'text' | 'boolean' | 'select';
  options?: string[];
  min?: number;
  max?: number;
  step?: number;
  unit?: string;
}

export interface ParameterCategory {
  name: string;
  description: string;
  parameters: Parameter[];
}

export interface EnergyLevel {
  manifold: number;
  magneticNumber: number;
  boundLevel: number;
  stateNumber: number;
  population: number;
  vibrational?: number;
  angularMomentum?: number;
  energy: number;
  electricShift: number;
  magneticShift: number;
}

export interface Transition {
  upperState: {
    manifold: number;
    magneticNumber: number;
    boundLevel: number;
    stateNumber: number;
  };
  lowerState: {
    manifold: number;
    magneticNumber: number;
    boundLevel: number;
    stateNumber: number;
  };
  upperEnergy: number;
  lowerEnergy: number;
  dipoleStrength: number;
}