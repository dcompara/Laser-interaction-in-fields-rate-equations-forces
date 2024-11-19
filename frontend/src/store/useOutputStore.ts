import { create } from 'zustand';

interface OutputSettings {
  graphics: {
    displaySize: number;
    waitTime: number;
    enabled: boolean;
    rotation: {
      x: number;
      y: number;
      z: number;
    };
    vMax: number;
    twoJMax: number;
  };
  output: {
    excludedManifold: number;
    studyLevel: number;
  };
  files: {
    dataCardOutput: number;
    paramScannedOutput: boolean;
  };
}

interface OutputStore {
  settings: OutputSettings;
  updateSettings: (settings: OutputSettings) => void;
}

const defaultSettings: OutputSettings = {
  graphics: {
    displaySize: 10e-2,
    waitTime: 1e-1,
    enabled: true,
    rotation: {
      x: 1,
      y: 0,
      z: 0,
    },
    vMax: 2,
    twoJMax: 4,
  },
  output: {
    excludedManifold: -10,
    studyLevel: -1,
  },
  files: {
    dataCardOutput: 2,
    paramScannedOutput: true,
  },
};

export const useOutputStore = create<OutputStore>((set) => ({
  settings: defaultSettings,
  updateSettings: (settings) => set({ settings }),
}));