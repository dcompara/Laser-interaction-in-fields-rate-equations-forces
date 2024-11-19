import { create } from 'zustand';

interface Vector3D {
  x: number;
  y: number;
  z: number;
}

interface FemtosecondParams {
  offsetFreq: number;
  repRate: number;
  combLineSpacing: number;
}

export interface LaserState {
  position: Vector3D;
  direction: Vector3D;
  waist: number;
  energy: number;
  linewidth: number;
  power: number;
  polarization: {
    leftCircular: number;
    rightCircular: number;
  };
  polarAngle: number;
  type: string;
  coherentLaser: number;
  femtosecond?: FemtosecondParams;
}

interface LaserGlobalSettings {
  powerScale: number;
  detuningOffset: number;
  linewidthScale: number;
  forcedRates: boolean;
}

interface LaserStore {
  lasers: LaserState[];
  globalSettings: LaserGlobalSettings;
  addLaser: () => void;
  removeLaser: (index: number) => void;
  updateLaser: (index: number, laser: Partial<LaserState>) => void;
  updateGlobalSettings: (settings: Partial<LaserGlobalSettings>) => void;
}

const defaultLaser: LaserState = {
  position: { x: 0, y: 0, z: 0 },
  direction: { x: 1, y: 0, z: 0 },
  waist: 3e-3,
  energy: 0.92,
  linewidth: 10,
  power: 1e-10,
  polarization: {
    leftCircular: 0.707107,
    rightCircular: 0.707107,
  },
  polarAngle: 0,
  type: '5',
  coherentLaser: -1,
};

const defaultGlobalSettings: LaserGlobalSettings = {
  powerScale: 1.0,
  detuningOffset: 0.0,
  linewidthScale: 1.0,
  forcedRates: false,
};

export const useLaserStore = create<LaserStore>((set) => ({
  lasers: [defaultLaser],
  globalSettings: defaultGlobalSettings,
  addLaser: () =>
    set((state) => ({
      lasers: [...state.lasers, { ...defaultLaser }],
    })),
  removeLaser: (index) =>
    set((state) => ({
      lasers: state.lasers.filter((_, i) => i !== index),
    })),
  updateLaser: (index, laser) =>
    set((state) => ({
      lasers: state.lasers.map((l, i) =>
        i === index ? { ...l, ...laser } : l
      ),
    })),
  updateGlobalSettings: (settings) =>
    set((state) => ({
      globalSettings: { ...state.globalSettings, ...settings },
    })),
}));