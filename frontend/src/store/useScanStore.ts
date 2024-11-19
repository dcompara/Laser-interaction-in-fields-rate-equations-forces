import { create } from 'zustand';

export interface ScanParameter {
  name: string;
  minValue: number;
  maxValue: number;
  steps: number;
  isScanned: boolean;
  isTimeDependent: boolean;
  tau?: number;
}

interface ScanState {
  tauModif: number;
  isRandomScan: boolean;
  parameters: ScanParameter[];
  addParameter: (param: ScanParameter) => void;
  removeParameter: (index: number) => void;
  updateParameter: (index: number, param: Partial<ScanParameter>) => void;
  updateTauModif: (value: number) => void;
  updateIsRandomScan: (value: boolean) => void;
}

const defaultParameters: ScanParameter[] = [
  {
    name: 'Offset_Detuning_cm',
    minValue: -0.01,
    maxValue: 0.01,
    steps: 20,
    isScanned: false,
    isTimeDependent: false
  },
  {
    name: 'scale_Gamma',
    minValue: 0.3,
    maxValue: 0.6,
    steps: 3,
    isScanned: false,
    isTimeDependent: false
  },
  {
    name: 'Tau_Modif',
    minValue: 0.5e-3,
    maxValue: 2e-3,
    steps: 2,
    isScanned: false,
    isTimeDependent: false
  },
  {
    name: 'E_z',
    minValue: 0,
    maxValue: 1e6,
    steps: 1000,
    isScanned: false,
    isTimeDependent: false
  },
  {
    name: 'scale_Power',
    minValue: 0.1,
    maxValue: 5.1,
    steps: 10,
    isScanned: false,
    isTimeDependent: false
  }
];

export const useScanStore = create<ScanState>((set) => ({
  tauModif: 1e-3,
  isRandomScan: false,
  parameters: defaultParameters,
  addParameter: (param) => set((state) => ({
    parameters: [...state.parameters, param]
  })),
  removeParameter: (index) => set((state) => ({
    parameters: state.parameters.filter((_, i) => i !== index)
  })),
  updateParameter: (index, param) => set((state) => ({
    parameters: state.parameters.map((p, i) =>
      i === index ? { ...p, ...param } : p
    )
  })),
  updateTauModif: (value) => set({ tauModif: value }),
  updateIsRandomScan: (value) => set({ isRandomScan: value })
}));