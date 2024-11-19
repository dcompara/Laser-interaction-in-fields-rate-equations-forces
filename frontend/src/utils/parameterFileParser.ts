import { ParticleState } from '../store/useParticleStore';
import { LaserState } from '../store/useLaserStore';
import { FieldState } from '../store/useFieldStore';
import { AlgorithmSettings } from '../store/useAlgorithmStore';
import { OutputSettings } from '../store/useOutputStore';

export function parseParameters(content: string): {
  particles: Partial<ParticleState>[];
  lasers: Partial<LaserState>[];
  fields: Partial<FieldState>;
  algorithm: Partial<AlgorithmSettings>;
  output: Partial<OutputSettings>;
  laserGlobalSettings: any;
  scanSettings: any;
} {
  const lines = content.split('\n');
  const params = new Map<string, string>();

  // Parse parameters into a map
  lines.forEach(line => {
    if (line.startsWith('@')) {
      const [param, ...valueParts] = line.split('\t');
      const value = valueParts.join('\t').trim();
      params.set(param.trim(), value);
    }
  });

  return {
    particles: parseParticles(params),
    lasers: parseLasers(params),
    fields: parseFields(params),
    algorithm: parseAlgorithm(params),
    output: parseOutput(params),
    laserGlobalSettings: parseLaserGlobalSettings(params),
    scanSettings: parseScanSettings(content),
  };
}

function parseParticles(params: Map<string, string>): Partial<ParticleState>[] {
  const particles: Partial<ParticleState>[] = [];
  const particleCount = Number(params.get('@Nb_type_of_Mol') || '1');

  for (let i = 0; i < particleCount; i++) {
    particles.push({
      name: params.get(`@Nom_Mol[${i}]`) || 'Cs',
      count: Number(params.get(`@N_Mol[${i}]`) || '100'),
      temperature: {
        x: Number(params.get(`@Temp_ini_x[${i}]`) || '0.00001'),
        y: Number(params.get(`@Temp_ini_y[${i}]`) || '1e-10'),
        z: Number(params.get(`@Temp_ini_z[${i}]`) || '1e-10'),
      },
      initialization: {
        x: Number(params.get(`@Procedure_init_x[${i}]`) || '-1'),
        y: Number(params.get(`@Procedure_init_y[${i}]`) || '0'),
        z: Number(params.get(`@Procedure_init_z[${i}]`) || '0'),
      },
      size: {
        x: Number(params.get(`@size_x[${i}]`) || '1e-3'),
        y: Number(params.get(`@size_y[${i}]`) || '1e-3'),
        z: Number(params.get(`@size_z[${i}]`) || '1e-3'),
      },
      position: {
        x: Number(params.get(`@offset_x[${i}]`) || '0'),
        y: Number(params.get(`@offset_y[${i}]`) || '0'),
        z: Number(params.get(`@offset_z[${i}]`) || '0'),
      },
      velocity: {
        x: Number(params.get(`@v0_x[${i}]`) || '0'),
        y: Number(params.get(`@v0_y[${i}]`) || '0'),
        z: Number(params.get(`@v0_z[${i}]`) || '0'),
      },
    });
  }

  return particles;
}

function parseLaserGlobalSettings(params: Map<string, string>) {
  return {
    powerScale: Number(params.get('@scale_Power') || '1.0'),
    detuningOffset: Number(params.get('@Offset_Detuning_cm') || '0.0'),
    linewidthScale: Number(params.get('@scale_Gamma') || '1.0'),
    forcedRates: params.get('@is_forced_rates') === '1',
  };
}

function parseLasers(params: Map<string, string>): Partial<LaserState>[] {
  const lasers: Partial<LaserState>[] = [];
  const laserCount = Number(params.get('@Nb_laser') || '1');

  for (let i = 0; i < laserCount; i++) {
    const laser: Partial<LaserState> = {
      position: {
        x: Number(params.get(`@waist_pos_x[${i}]`) || '0'),
        y: Number(params.get(`@waist_pos_y[${i}]`) || '0'),
        z: Number(params.get(`@waist_pos_z[${i}]`) || '0'),
      },
      direction: {
        x: Number(params.get(`@direction_x[${i}]`) || '1'),
        y: Number(params.get(`@direction_y[${i}]`) || '0'),
        z: Number(params.get(`@direction_z[${i}]`) || '0'),
      },
      waist: Number(params.get(`@waist[${i}]`) || '3e-3'),
      energy: Number(params.get(`@Energie_cm[${i}]`) || '0.92'),
      linewidth: Number(params.get(`@Gamma_L_MHz[${i}]`) || '10'),
      power: Number(params.get(`@Power[${i}]`) || '1e-10'),
      type: params.get(`@type_laser[${i}]`) || '5',
      polarAngle: Number(params.get(`@polar_angle_degree[${i}]`) || '0'),
      coherentLaser: Number(params.get(`@coherent_avec_laser_num[${i}]`) || '-1'),
    };

    if (laser.type === '1') {
      laser.femtosecond = {
        offsetFreq: Number(params.get(`@nu_offset_MHz[${i}]`) || '0.0'),
        repRate: Number(params.get(`@nu_repetition_MHz[${i}]`) || '80.0'),
        combLineSpacing: Number(params.get(`@nu_individual_comb_line_MHz[${i}]`) || '80.0'),
      };
    }

    lasers.push(laser);
  }

  return lasers;
}

function parseFields(params: Map<string, string>): Partial<FieldState> {
  const fields: Partial<FieldState> = {
    typeForStateShift: Number(params.get('@type_of_field_for_internal_state_shift') || '0'),
    magnetic: {
      type: Number(params.get('@type_field_read_B') || '0'),
      B: {
        x: Number(params.get('@B_x') || '1e-10'),
        y: Number(params.get('@B_y') || '2e-10'),
        z: Number(params.get('@B_z') || '0'),
      },
      gradient: {
        x: Number(params.get('@grad_B_x') || '0'),
        y: Number(params.get('@grad_B_y') || '0'),
        z: Number(params.get('@grad_B_z') || '0'),
      },
      secondGradient: {
        x: Number(params.get('@grad_grad_B_x') || '0'),
        y: Number(params.get('@grad_grad_B_y') || '0'),
        z: Number(params.get('@grad_grad_B_z') || '0'),
      },
    },
    electric: {
      type: Number(params.get('@type_field_read_E') || '0'),
      E: {
        x: Number(params.get('@E_x') || '0'),
        y: Number(params.get('@E_y') || '0'),
        z: Number(params.get('@E_z') || '0'),
      },
      gradient: {
        x: Number(params.get('@grad_E_x') || '0'),
        y: Number(params.get('@grad_E_y') || '0'),
        z: Number(params.get('@grad_E_z') || '0'),
      },
      secondGradient: {
        x: Number(params.get('@grad_grad_E_x') || '0'),
        y: Number(params.get('@grad_grad_E_y') || '0'),
        z: Number(params.get('@grad_grad_E_z') || '0'),
      },
    },
  };

  if (fields.magnetic?.type === 1) {
    fields.magnetic.helmholtz = {
      coilCount: Number(params.get('@Nb_bobines') || '5'),
      isHelmholtz: Number(params.get('@is_Helmholtz') || '0'),
      coilSpacing: Number(params.get('@gap_bobines') || '4e-2'),
      current: Number(params.get('@courant_bobines') || '16000'),
      radius: Number(params.get('@rayon_bobines') || '1e-2'),
    };
  }

  return fields;
}

function parseAlgorithm(params: Map<string, string>): Partial<AlgorithmSettings> {
  return {
    monteCarlo: {
      algorithm: Number(params.get('@Choix_algorithme_Monte_Carlo') || '0'),
    },
    nBody: {
      algorithm: Number(params.get('@Choix_algorithme_N_corps') || '1'),
      epsilon: Number(params.get('@choix_epsilon') || '1e-6'),
    },
    randomSeed: Number(params.get('@Seed_Init_Random_Number_Generator') || '2'),
    diagonalization: Number(params.get('@is_Levels_Lines_Diagonalized') || '0'),
  };
}

function parseOutput(params: Map<string, string>): Partial<OutputSettings> {
  return {
    graphics: {
      displaySize: Number(params.get('@SIZE_affichage') || '10e-2'),
      waitTime: Number(params.get('@t_wait_affichage') || '1e-1'),
      enabled: params.get('@Graphics') === '1',
      rotation: {
        x: Number(params.get('@rot_axe_x') || '1'),
        y: Number(params.get('@rot_axe_y') || '0'),
        z: Number(params.get('@rot_axe_z') || '0'),
      },
      vMax: Number(params.get('@v_max') || '2'),
      twoJMax: Number(params.get('@two_J_max') || '4'),
    },
    files: {
      dataCardOutput: Number(params.get('@is_DataCard_out') || '2'),
      paramScannedOutput: params.get('@is_param_scanned_out') === '1',
      levelFile: params.get('@nom_file_Levels') || 'Data/n/n_levels.dat',
      lineFile: params.get('@nom_file_Lines') || 'Data/n/n_lines.dat',
      laserSpectrumFile: params.get('@nom_file_Laser_Spectrum') || 'Data/n/Laser_Spectrum.dat',
      laserIntensityFile: params.get('@nom_file_Laser_Intensity') || 'Data/n/Laser_Intensity.dat',
      outputDataFile: params.get('@nom_sortie_donnees') || 'Data/donnee_Mol.dat',
      outputRateFile: params.get('@nom_sortie_rate') || 'Data/sortie_rate.dat',
      outputDataCardFile: params.get('@nom_sortie_donnees_Data') || 'Data/data_card.dat',
      randomGenFile: params.get('@nom_fichier_random_gen') || 'Data/random_gen.txt',
    },
  };
}

function parseScanSettings(content: string): any {
  const scanParams: any[] = [];
  let inFitParams = false;
  const lines = content.split('\n');

  const tauModif = lines.find(line => line.startsWith('@Tau_Modif'))?.split('\t')[1] || '1e-3';
  const isRandomScan = lines.find(line => line.startsWith('@is_Scan_Random'))?.split('\t')[1] === 'true';

  for (const line of lines) {
    if (line.includes('BEGIN_OF_FITPARAMS')) {
      inFitParams = true;
      continue;
    }
    if (line.includes('END_OF_FITPARAMS')) {
      inFitParams = false;
      continue;
    }
    if (inFitParams && line.startsWith('@SCAN_')) {
      const [name, minValue, maxValue, steps, isScanned, isTimeDependent, tau] = line
        .substring(6)
        .split('\t')
        .map(s => s.trim());
      
      scanParams.push({
        name,
        minValue: Number(minValue),
        maxValue: Number(maxValue),
        steps: Number(steps),
        isScanned: isScanned === 'true',
        isTimeDependent: isTimeDependent === 'true',
        tau: tau ? Number(tau) : undefined,
      });
    }
  }

  return {
    tauModif,
    isRandomScan,
    scanParams,
  };
}