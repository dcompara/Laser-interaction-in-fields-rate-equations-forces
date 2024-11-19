import { ParticleState } from '../store/useParticleStore';
import { LaserState } from '../store/useLaserStore';
import { FieldState } from '../store/useFieldStore';
import { AlgorithmSettings } from '../store/useAlgorithmStore';
import { OutputSettings } from '../store/useOutputStore';

export function generateHeader() {
  return `/************************* Parameter File ***************************

This file is used to define the parameters for the program.
Example of usage:
double Name_Parameter = params.LocateParam("Name_Parameter")->val

# 0 = false ; 1 = true

***********************************************************************/\n\n`;
}

export function generateParticleSection(particles: ParticleState[]) {
  let content = `########## \tParticles\t###########################################################################\n`;
  content += `########## \tAll values are in SI units (unless explicitly stated otherwise) \t########################\n#\n\n`;
  
  content += `// Number of particle types. Currently, only the first type can be laser cooled (Levels and Lines are specific to this type).\n`;
  content += `@Nb_type_of_Mol\t${particles.length}\n\n`;
  
  content += `// Choice of particle name (e.g., BaF, Cs2, NH, Cs, CO, Li6Cs, Laminus, Li7Cs, Rb85Cs, Rb87Cs, Ps, C2minus, Ps_minus, P_bar, Osminus).\n`;
  content += `// Note: This affects only the mass; the corresponding file names must still be updated manually.\n\n`;

  particles.forEach((particle, index) => {
    content += `// ${index === 0 ? 'First' : 'Second'} type of particle\n`;
    content += `// Particles range from Mol[${index}] to Mol[Nom_Mol[${index}]-1].\n`;
    content += `@Nom_Mol[${index}]\t${particle.name} // Name of the ${index === 0 ? 'first' : 'second'} particle type\n`;
    content += `@N_Mol[${index}]\t${particle.count} // Number of particles to be laser cooled\n`;
    content += `@Temp_ini_x[${index}]\t${particle.temperature.x} // Initial temperature along x-axis\n`;
    content += `@Temp_ini_y[${index}]\t${particle.temperature.y} // Initial temperature along y-axis\n`;
    content += `@Temp_ini_z[${index}]\t${particle.temperature.z} // Initial temperature along z-axis\n\n`;

    content += `@Procedure_init_x[${index}]\t${particle.initialization.x} // Positioning procedure along x-axis\n`;
    content += `@Procedure_init_y[${index}]\t${particle.initialization.y} // Positioning procedure along y-axis\n`;
    content += `@Procedure_init_z[${index}]\t${particle.initialization.z} // Positioning procedure along z-axis\n\n`;

    content += `// Size (x, y, z) if fixed size positioning is chosen\n`;
    content += `@size_x[${index}]\t${particle.size.x} // Size along x-axis\n`;
    content += `@size_y[${index}]\t${particle.size.y} // Size along y-axis\n`;
    content += `@size_z[${index}]\t${particle.size.z} // Size along z-axis\n\n`;

    content += `// Offset added to the initial randomized position\n`;
    content += `@offset_x[${index}]\t${particle.position.x} // Offset along x-axis\n`;
    content += `@offset_y[${index}]\t${particle.position.y} // Offset along y-axis\n`;
    content += `@offset_z[${index}]\t${particle.position.z} // Offset along z-axis\n\n`;

    content += `// Initial velocity added to the randomized velocity\n`;
    content += `@v0_x[${index}]\t${particle.velocity.x} // Initial velocity along x-axis\n`;
    content += `@v0_y[${index}]\t${particle.velocity.y} // Initial velocity along y-axis\n`;
    content += `@v0_z[${index}]\t${particle.velocity.z} // Initial velocity along z-axis\n\n`;
  });

  return content;
}

export function generateLaserSection(lasers: LaserState[], globalSettings: any) {
  let content = `########## \tLASERS\t########################################\n#\n`;
  
  content += `// Scaling factor for all laser powers\n`;
  content += `@scale_Power\t${globalSettings.powerScale}\n\n`;
  
  content += `// Additive offset to the frequency of all lasers\n`;
  content += `@Offset_Detuning_cm\t${globalSettings.detuningOffset}\n\n`;
  
  content += `// Scaling factor for the spectral width of all lasers\n`;
  content += `@scale_Gamma\t${globalSettings.linewidthScale}\n\n`;
  
  content += `// Number of lasers used\n`;
  content += `@Nb_laser\t${lasers.length}\n\n`;
  
  content += `// Force rates not linked to lasers\n`;
  content += `@is_forced_rates\t${globalSettings.forcedRates ? '1' : '0'}\n\n`;

  lasers.forEach((laser, index) => {
    content += `# Laser ${index + 1}\n`;
    content += `@waist_pos_x[${index}]\t${laser.position.x}\n`;
    content += `@waist_pos_y[${index}]\t${laser.position.y}\n`;
    content += `@waist_pos_z[${index}]\t${laser.position.z}\n`;
    content += `@direction_x[${index}]\t${laser.direction.x}\n`;
    content += `@direction_y[${index}]\t${laser.direction.y}\n`;
    content += `@direction_z[${index}]\t${laser.direction.z}\n\n`;
    
    content += `@waist[${index}]\t${laser.waist}\n`;
    content += `@Energie_cm[${index}]\t${laser.energy}\n`;
    content += `@Gamma_L_MHz[${index}]\t${laser.linewidth}\n`;
    content += `@Power[${index}]\t${laser.power}\n`;
    content += `@polar_angle_degree[${index}]\t${laser.polarAngle}\n`;
    content += `@type_laser[${index}]\t${laser.type}\n`;
    content += `@coherent_avec_laser_num[${index}]\t${laser.coherentLaser}\n\n`;

    if (laser.type === '1') {
      content += `@nu_offset_MHz[${index}]\t${laser.femtosecond?.offsetFreq}\n`;
      content += `@nu_repetition_MHz[${index}]\t${laser.femtosecond?.repRate}\n`;
      content += `@nu_individual_comb_line_MHz[${index}]\t${laser.femtosecond?.combLineSpacing}\n\n`;
    }
  });

  return content;
}

export function generateFieldSection(fields: FieldState) {
  let content = `########## \tEXTERNAL FIELDS\t################################\n#\n`;
  
  content += `// Type of field for internal state shift (0: magnetic, 1: electric)\n`;
  content += `@type_of_field_for_internal_state_shift\t${fields.typeForStateShift}\n\n`;
  
  content += `// Field definition methods\n`;
  content += `@type_field_read_E\t${fields.electric.type}\n`;
  content += `@type_field_read_B\t${fields.magnetic.type}\n\n`;

  if (fields.magnetic.type === 1) {
    content += `## MAGNETIC FIELD (Helmholtz coils) ##\n`;
    content += `@Nb_bobines\t${fields.magnetic.helmholtz?.coilCount}\n`;
    content += `@is_Helmholtz\t${fields.magnetic.helmholtz?.isHelmholtz}\n`;
    content += `@gap_bobines\t${fields.magnetic.helmholtz?.coilSpacing}\n`;
    content += `@courant_bobines\t${fields.magnetic.helmholtz?.current}\n`;
    content += `@rayon_bobines\t${fields.magnetic.helmholtz?.radius}\n\n`;
  }

  content += `## MAGNETIC FIELD ##\n`;
  content += `@B_x\t${fields.magnetic.B.x}\n`;
  content += `@B_y\t${fields.magnetic.B.y}\n`;
  content += `@B_z\t${fields.magnetic.B.z}\n\n`;
  
  content += `@grad_B_x\t${fields.magnetic.gradient.x}\n`;
  content += `@grad_B_y\t${fields.magnetic.gradient.y}\n`;
  content += `@grad_B_z\t${fields.magnetic.gradient.z}\n\n`;
  
  content += `@grad_grad_B_x\t${fields.magnetic.secondGradient.x}\n`;
  content += `@grad_grad_B_y\t${fields.magnetic.secondGradient.y}\n`;
  content += `@grad_grad_B_z\t${fields.magnetic.secondGradient.z}\n\n`;

  content += `## ELECTRIC FIELD ##\n`;
  content += `@E_x\t${fields.electric.E.x}\n`;
  content += `@E_y\t${fields.electric.E.y}\n`;
  content += `@E_z\t${fields.electric.E.z}\n\n`;
  
  content += `@grad_E_x\t${fields.electric.gradient.x}\n`;
  content += `@grad_E_y\t${fields.electric.gradient.y}\n`;
  content += `@grad_E_z\t${fields.electric.gradient.z}\n\n`;
  
  content += `@grad_grad_E_x\t${fields.electric.secondGradient.x}\n`;
  content += `@grad_grad_E_y\t${fields.electric.secondGradient.y}\n`;
  content += `@grad_grad_E_z\t${fields.electric.secondGradient.z}\n\n`;

  return content;
}

export function generateOutputSection(settings: OutputSettings) {
  let content = `########## \tOUTPUT PARAMETERS\t##############################\n#\n`;
  
  content += `// Graphics settings\n`;
  content += `@SIZE_affichage\t${settings.graphics.displaySize}\n`;
  content += `@t_wait_affichage\t${settings.graphics.waitTime}\n`;
  content += `@Graphics\t${settings.graphics.enabled ? '1' : '0'}\n\n`;
  
  content += `// Rotation settings\n`;
  content += `@rot_axe_x\t${settings.graphics.rotation.x}\n`;
  content += `@rot_axe_y\t${settings.graphics.rotation.y}\n`;
  content += `@rot_axe_z\t${settings.graphics.rotation.z}\n\n`;
  
  content += `// Molecule drawing settings\n`;
  content += `@v_max\t${settings.graphics.vMax}\n`;
  content += `@two_J_max\t${settings.graphics.twoJMax}\n\n`;
  
  content += `// Output file settings\n`;
  content += `@is_DataCard_out\t${settings.files.dataCardOutput}\n`;
  content += `@is_param_scanned_out\t${settings.files.paramScannedOutput ? '1' : '0'}\n\n`;
  
  content += `// File paths\n`;
  content += `@nom_file_Levels\t${settings.files.levelFile}\n`;
  content += `@nom_file_Lines\t${settings.files.lineFile}\n`;
  content += `@nom_file_Laser_Spectrum\t${settings.files.laserSpectrumFile}\n`;
  content += `@nom_file_Laser_Intensity\t${settings.files.laserIntensityFile}\n`;
  content += `@nom_sortie_donnees\t${settings.files.outputDataFile}\n`;
  content += `@nom_sortie_rate\t${settings.files.outputRateFile}\n`;
  content += `@nom_sortie_donnees_Data\t${settings.files.outputDataCardFile}\n`;
  content += `@nom_fichier_random_gen\t${settings.files.randomGenFile}\n\n`;

  return content;
}

export function generateAlgorithmSection(settings: AlgorithmSettings) {
  let content = `########## \tALGORITHM SETTINGS\t##############################\n#\n`;
  
  content += `// Monte Carlo algorithm\n`;
  content += `@Choix_algorithme_Monte_Carlo\t${settings.monteCarlo.algorithm}\n\n`;
  
  content += `// N-body algorithm\n`;
  content += `@Choix_algorithme_N_corps\t${settings.nBody.algorithm}\n`;
  content += `@choix_epsilon\t${settings.nBody.epsilon}\n\n`;
  
  content += `// Random number generator\n`;
  content += `@Seed_Init_Random_Number_Generator\t${settings.randomSeed}\n\n`;
  
  content += `// Diagonalization\n`;
  content += `@is_Levels_Lines_Diagonalized\t${settings.diagonalization}\n\n`;

  return content;
}

export function generateScanSection(settings: { tauModif: string; isRandomScan: boolean; scanParams: any[] }) {
  let content = `########## \tSCAN PARAMETERS\t##############################\n#\n`;
  
  content += `// Scanning and time-dependent variables are defined in lists.\n`;
  content += `// The list starts with BEGIN_OF... and ends with END_OF....\n`;
  content += `// Each scanned variable must be named as SCAN_<VariableName>, followed by:\n`;
  content += `// - Minimum value, maximum value, number of steps (nbstep).\n`;
  content += `// - Two booleans: is_scanned (true/false) and is_time (true/false).\n`;
  content += `// - tau_var: Time constant for exponential changes if is_time is true.\n`;
  content += `// Nb_steps refers to the number of intervals:\n`;
  content += `// - 1 means two values (min and max), 2 means three values (min, mid, max), etc.\n`;
  content += `//\n`;
  content += `// Value at time t is calculated as:\n`;
  content += `// val_t0 = min_value + steps * (max_value - min_value) / nbstep\n`;
  content += `// tau = time rate of change → val = val_t0 * exp^(-t / tau).\n`;
  content += `//\n`;
  content += `// If tau is unspecified but is_time is true, Tau_Modif will be used as the default.\n`;
  content += `// If multiple tables are needed, specify the table number as [table_index].\n\n`;
  
  content += `# Default time constant for exponential modifications\n`;
  content += `@Tau_Modif ${settings.tauModif}\n\n`;
  
  content += `# Randomized or ordered scans\n`;
  content += `@is_Scan_Random ${settings.isRandomScan}\n\n`;
  
  content += `# Scanned variables\n`;
  content += `BEGIN_OF_FITPARAMS\n`;
  settings.scanParams.forEach(param => {
    content += `@SCAN_${param.name}\t${param.minValue}\t${param.maxValue}\t${param.steps}\t${param.isScanned}\t${param.isTimeDependent}${param.tau ? '\t' + param.tau : ''}\n`;
  });
  content += `END_OF_FITPARAMS\n\n`;

  return content;
}