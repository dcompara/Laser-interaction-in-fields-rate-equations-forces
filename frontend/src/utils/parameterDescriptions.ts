export const parameterDescriptions = {
  // Particle parameters
  'Nb_type_of_Mol': 'Number of different particle types. Currently only the first type is laser cooled.',
  'Nom_Mol': 'Particle name selection. The choice only affects mass, not the file names which need to be changed separately.',
  'N_Mol': 'Number of particles of this type that are laser cooled.',
  'Temp_ini': 'Initial temperature in each direction (x, y, z)',
  'Procedure_init': `Position initialization method:
    -1: Fixed size, ordered positions at start
    0: Fixed Gaussian size
    1: Linear magnetic potential (Laplace)
    2: Quadratic magnetic potential (neutral)
    3: Quadratic potential (charged)
    4: Perfect ordered Gaussian in velocity
    5: Effusive beam (orientation defined by component)
    6: Collimated effusive beam (orientation defined by component)`,
  'size': 'Fixed size in each direction if fixed size is chosen',
  'offset': 'Initial position offset added to random position',
  'v0': 'Initial velocity offset added to random velocity',

  // Time parameters
  'dt_dyn_epsilon_param': 'Dynamic time step size for convergence control. Typical values: 0.001*waist/velocity or 0.001*cyclotron period',
  't_fin': 'Final simulation time',
  'dt_dia': 'Time interval between diagnostic outputs',
  'dt_out': 'Time interval between snapshot outputs',

  // Graphics parameters
  'SIZE_affichage': 'Display size for visualization',
  't_wait_affichage': 'Wait time between displays to control display speed',
  'Graphics': 'Enable/disable graphics output',
  'rot_axe': 'Rotation axis for 3D visualization. Default (1,0,0) = x right, y in screen, z up',
  'v_max': 'Maximum vibration level for molecule drawing',
  'two_J_max': 'Maximum rotation level (2J) for molecule drawing',

  // Output parameters
  'num_manifold_not_studied': 'Manifold number to exclude from statistics. Use -10 to include all.',
  'num_niveau_etudie': 'Level number for statistics. Use -1 for all molecules.',
  'is_DataCard_out': 'Output datacard: 0=none, 1=single file, 2=split files',
  'is_param_scanned_out': 'Enable output of scanned parameters',

  // File parameters
  'nom_file_Levels': 'File containing energy levels data',
  'nom_file_Lines': 'File containing transition lines data',
  'nom_file_Laser_Spectrum': 'Base name for laser spectrum files',
  'nom_file_Laser_Intensity': 'Base name for laser timing files',
  'nom_sortie_donnees': 'Output file for molecule data',
  'nom_sortie_rate': 'Output file for rate data',
  'nom_sortie_donnees_Data': 'Output file for data card',
  'nom_fichier_random_gen': 'File for random number generator state',

  // Field parameters
  'type_of_field_for_internal_state_shift': 'Choose between magnetic (0) or electric (1) field for internal state shifts',
  'type_field_read_E': `Electric field definition method:
    0: 2nd order plus nth order (analytical)
    1: Helmholtz coils
    2: 2D cylindrical symmetry (r,z,Fr,Fz)
    3: 2D cylindrical with derivatives
    4: 3D field map`,
  'type_field_read_B': 'Magnetic field definition method (same options as electric)',
  'B': 'Magnetic field components (x, y, z). Never use 0, use small value like 1e-10 to maintain quantization axis',
  'grad_B': 'Magnetic field gradient',
  'grad_grad_B': 'Second-order magnetic field gradient',
  'E': 'Electric field components (x, y, z)',
  'grad_E': 'Electric field gradient',
  'grad_grad_E': 'Second-order electric field gradient',

  // Laser parameters
  'scale_Power': 'Global scaling factor for all laser powers',
  'Offset_Detuning_cm': 'Global frequency offset for all lasers',
  'scale_Gamma': 'Global scaling factor for laser spectral widths',
  'waist_pos': 'Laser waist position',
  'direction': 'Laser propagation direction',
  'waist': 'Laser beam waist (1/e² radius)',
  'Energie_cm': 'Laser energy in cm⁻¹',
  'Gamma_L_MHz': 'Laser linewidth in MHz',
  'Power': 'Laser power in Watts',
  'type_laser': `Laser type:
    -1: Spontaneous emission
    0: CW
    1: Femtosecond
    2: Multi-mode
    3: Pulse
    4: Shaped (obsolete)
    5: Gaussian (most common)
    6: Lorentzian
    7: Frequency comb
    8: Pseudo-BBR
    9: Field ionization
    10: Collimated rectangular super-Gaussian`,
  'coherent_avec_laser_num': 'Laser interference settings. -1 for no interference, or number of laser to interfere with',

  // Algorithm parameters
  'Choix_algorithme_Monte_Carlo': `Monte Carlo algorithm selection:
    -1: None (no rate calculation)
    0: Kinetic Monte Carlo
    1: Random Selection Method
    2: First Reaction Method
    3: Fast Rough Method`,
  'Choix_algorithme_N_corps': `N-body algorithm selection:
    -1: None (with photon recoil)
    1: Verlet acceleration (no dipole force)
    2: Verlet potential (with dipole potential)
    6: High-order Verlet potential
    3: Yoshida6 acceleration
    4: Yoshida6 potential
    7: Boris-Buneman (charged particles in magnetic field)`,
  'choix_epsilon': 'Position epsilon for potential derivative calculation. Aim for 1/100 potential variation.',
  'Seed_Init_Random_Number_Generator': 'Random seed for reproducible simulations. Same value gives same sequence.',

  // Scan parameters
  'Tau_Modif': 'Default time constant for parameter modifications if not specified',
  'is_Scan_Random': 'Enable random scanning of parameters instead of ordered',
  'SCAN_parameters': 'Parameters that can be scanned or time-varied. Format: name, min, max, steps, is_scanned, is_time_dependent, tau'
};