import { create } from 'zustand';

interface Vector3D {
  x: number;
  y: number;
  z: number;
}

export interface ParticleState {
  name: string;
  count: number;
  temperature: Vector3D;
  initialization: Vector3D;
  size: Vector3D;
  position: Vector3D;
  velocity: Vector3D;
}

interface ParticleStore {
  particles: ParticleState[];
  addParticle: () => void;
  removeParticle: (index: number) => void;
  updateParticle: (index: number, particle: Partial<ParticleState>) => void;
}

const defaultParticle: ParticleState = {
  name: 'Cs',
  count: 100,
  temperature: {
    x: 0.00001,
    y: 1e-10,
    z: 1e-10,
  },
  initialization: {
    x: -1,
    y: 0,
    z: 0,
  },
  size: {
    x: 1e-3,
    y: 1e-3,
    z: 1e-3,
  },
  position: {
    x: -4e-2,
    y: 0,
    z: 0,
  },
  velocity: {
    x: 10.0,
    y: 0.0,
    z: 0.0,
  },
};

export const useParticleStore = create<ParticleStore>((set) => ({
  particles: [defaultParticle],
  addParticle: () =>
    set((state) => ({
      particles: [...state.particles, { ...defaultParticle }],
    })),
  removeParticle: (index) =>
    set((state) => ({
      particles: state.particles.filter((_, i) => i !== index),
    })),
  updateParticle: (index, particle) =>
    set((state) => ({
      particles: state.particles.map((p, i) =>
        i === index ? { ...p, ...particle } : p
      ),
    })),
}));