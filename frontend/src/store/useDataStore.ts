import { create } from 'zustand';
import { EnergyLevel, Transition } from '../types';

interface DataState {
  levels: EnergyLevel[];
  transitions: Transition[];
  setLevels: (levels: EnergyLevel[]) => void;
  setTransitions: (transitions: Transition[]) => void;
  selectedLevel: EnergyLevel | null;
  setSelectedLevel: (level: EnergyLevel | null) => void;
}

export const useDataStore = create<DataState>((set) => ({
  levels: [],
  transitions: [],
  selectedLevel: null,
  setLevels: (levels) => set({ levels }),
  setTransitions: (transitions) => set({ transitions }),
  setSelectedLevel: (level) => set({ selectedLevel: level }),
}));