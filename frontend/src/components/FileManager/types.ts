export interface FileConfig {
  name: string;
  description: string;
  path: string;
  type: 'input' | 'output';
}

export const fileConfigs: FileConfig[] = [
  {
    name: 'Energy Levels',
    description: 'Contains energy level data',
    path: 'Data/n/n_levels.dat',
    type: 'input'
  },
  {
    name: 'Transition Lines',
    description: 'Contains transition line data',
    path: 'Data/n/n_lines.dat',
    type: 'input'
  },
  {
    name: 'Laser Spectrum',
    description: 'Laser spectrum data',
    path: 'Data/n/Laser_Spectrum.dat',
    type: 'input'
  },
  {
    name: 'Laser Intensity',
    description: 'Laser timing data',
    path: 'Data/n/Laser_Intensity.dat',
    type: 'input'
  },
  {
    name: 'Molecule Data',
    description: 'Output file for molecule data',
    path: 'Data/donnee_Mol.dat',
    type: 'output'
  },
  {
    name: 'Rate Data',
    description: 'Output file for rate data',
    path: 'Data/sortie_rate.dat',
    type: 'output'
  },
  {
    name: 'Data Card',
    description: 'Output file for data card',
    path: 'Data/data_card.dat',
    type: 'output'
  },
  {
    name: 'Random Generator State',
    description: 'Random number generator state',
    path: 'Data/random_gen.txt',
    type: 'output'
  }
];