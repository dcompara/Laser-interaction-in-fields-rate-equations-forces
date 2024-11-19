import { create } from 'zustand';

interface AlgorithmSettings {
  monteCarlo: {
    algorithm: number;
  };
  nBody: {
    algorithm: number;
    epsilon: number;
  };
  time: {
    dynamicTimeStep: number;
    finalTime: number;
    diagnosticInterval: number;
    outputInterval: number;
  };
  randomSeed: number;
}

interface AlgorithmStore {
  settings: AlgorithmSettings;
  updateSettings: (settings: AlgorithmSettings) => void;
}

const defaultSettings: AlgorithmSettings = {
  monteCarlo: {
    algorithm: 0, // Kinetic Monte Carlo
  },
  nBody: {
    algorithm: 1, // Verlet (Acceleration)
    epsilon: 1e-6,
  },
  time: {
    dynamicTimeStep: 1e-6,
    finalTime: 50000e-6,
    diagnosticInterval: 100e-6,
    outputInterval: 100e-6,
  },
  randomSeed: 2,
};

export const useAlgorithmStore = create<AlgorithmStore>((set) => ({
  settings: defaultSettings,
  updateSettings: (settings) => set({ settings }),
}));