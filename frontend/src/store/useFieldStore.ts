import { create } from 'zustand';

interface Vector3D {
  x: number;
  y: number;
  z: number;
}

interface HelmholtzConfig {
  coilCount: number;
  isHelmholtz: number;
  coilSpacing: number;
  current: number;
  radius: number;
}

interface FieldDefinition {
  type: number;
  field: Vector3D;
  gradient: Vector3D;
  secondGradient: Vector3D;
  helmholtz?: HelmholtzConfig;
  gridConfig?: {
    filePath: string;
    type: '2D-cylindrical' | '2D-cylindrical-derivatives' | '3D';
  };
}

export interface FieldState {
  typeForStateShift: number;
  magnetic: FieldDefinition;
  electric: FieldDefinition;
}

interface FieldStore {
  fields: FieldState;
  updateFields: (fields: Partial<FieldState>) => void;
  updateMagneticField: (field: Partial<FieldDefinition>) => void;
  updateElectricField: (field: Partial<FieldDefinition>) => void;
}

const defaultFields: FieldState = {
  typeForStateShift: 0,
  magnetic: {
    type: 0,
    field: { x: 1e-10, y: 2e-10, z: 0 },
    gradient: { x: 0, y: 0, z: 0 },
    secondGradient: { x: 0, y: 0, z: 0 },
  },
  electric: {
    type: 0,
    field: { x: 0, y: 0, z: 0 },
    gradient: { x: 0, y: 0, z: 0 },
    secondGradient: { x: 0, y: 0, z: 0 },
  },
};

export const useFieldStore = create<FieldStore>((set) => ({
  fields: defaultFields,
  updateFields: (newFields) => set((state) => ({
    fields: { ...state.fields, ...newFields }
  })),
  updateMagneticField: (field) => set((state) => ({
    fields: {
      ...state.fields,
      magnetic: { ...state.fields.magnetic, ...field }
    }
  })),
  updateElectricField: (field) => set((state) => ({
    fields: {
      ...state.fields,
      electric: { ...state.fields.electric, ...field }
    }
  })),
}));