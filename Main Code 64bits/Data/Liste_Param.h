/************************* Parameter File ***************************

This file is used to define the parameters for the program.
Example of usage:
double Name_Parameter = params.LocateParam("Name_Parameter")->val

# 0 = false ; 1 = true

***********************************************************************/
########## 	Particles  	###########################################################################
########## 	All values are in SI units (unless explicitly stated otherwise) 	########################
#

// Number of particle types. Currently, only the first type can be laser cooled (Levels and Lines are specific to this type).
@Nb_type_of_Mol	1

// Choice of particle name (e.g., BaF, Cs2, NH, Cs, CO, Li6Cs, Laminus, Li7Cs, Rb85Cs, Rb87Cs, Ps, C2minus, Ps_minus, P_bar, Osminus).
// Note: This affects only the mass; the corresponding file names must still be updated manually.

// First type of particle
// Particles range from Mol[0] to Mol[Nom_Mol[0]-1].
@Nom_Mol[0]	n // Name of the first particle type
@N_Mol[0]  100 // Number of particles to be laser cooled
@Temp_ini_x[0] 0.00001 // Initial temperature along x-axis
@Temp_ini_y[0] 1e-10   // Initial temperature along y-axis
@Temp_ini_z[0] 1e-10   // Initial temperature along z-axis


// Choice of particle positioning: fixed size (sigma_pos) or based on density
// -1: Fixed size and ordered positions initially (one axis is fixed, others are randomized).
//  0: Fixed size determined by @size parameters (Gaussian distribution).
//  1: Magnetic linear potential --> Laplace distribution (coefficient determined by F1).
//  2: Magnetic quadratic potential for neutral particles --> Gaussian distribution.
//  3: Magnetic quadratic potential for charged particles --> Gaussian distribution.
//  4: Quadratic electric potential (linear electric field) for CHARGED particles (e.g., in a Paul trap) --> Gaussian distribution.
//  5: Perfectly ordered Gaussian in velocity (randomized positions).
//  6: Effusive beam, where only one component (x, y, or z) is set to 5 to orient the beam.
//  7: Collimated effusive beam, where one component (x, y, or z) is set to 6 for orientation, and others are typically 0 for transverse temperature.

@Procedure_init_x[0]   -1 // Positioning procedure along x-axis
@Procedure_init_y[0]    0 // Positioning procedure along y-axis
@Procedure_init_z[0]    0 // Positioning procedure along z-axis

// Size (x, y, z) if fixed size positioning is chosen
@size_x[0]  1e-3 // Size along x-axis
@size_y[0]  1e-3 // Size along y-axis
@size_z[0]  1e-3 // Size along z-axis

// Offset added to the initial randomized position
@offset_x[0]   -4e-2 // Offset along x-axis
@offset_y[0]    0    // Offset along y-axis
@offset_z[0]    0    // Offset along z-axis

// Initial velocity added to the randomized velocity
@v0_x[0]      10.0  // Initial velocity along x-axis
@v0_y[0]       0.0  // Initial velocity along y-axis
@v0_z[0]       0.0  // Initial velocity along z-axis


#
#2nd type of particle
#
// so Mol[Nom_Mol[0]] to Mol[Nom_Mol[0]+Nom_Mol[1]-1]
@Nom_Mol[1]	P_bar
@N_Mol[1]  10
@Temp_ini_x[1] 4
@Temp_ini_y[1] 4
@Temp_ini_z[1] 4
@Procedure_init_x[1]   0
@Procedure_init_y[1]   0
@Procedure_init_z[1]   0
@size_x[1]	1e-4
@size_y[1]	1e-4
@size_z[1]	6e-4
@offset_x[1]	0.
@offset_y[1]	0.
@offset_z[1]	0.
@v0_x[1]	0.
@v0_y[1]	0.
@v0_z[1]	0.


##########   KMC Time and Output Parameters   ####################
#
// Parameter to control the dynamic time step size (in seconds) for convergence checks.
// Typical values:
// - 0.001 * waist/velocity (or 0.001 * wavelength/velocity for lattices).
// - Near 0.001 * cyclotron period (2π * m / (q * B)) for Leapfrog (in magnetic field cases).
// - 0.1 * m / (q * B) for Boris algorithm (e.g., 10^-8 at 0.0001 T for 3me mass; 2e-8 for C2- at 1 Tesla).
// A good test is to remove lasers and check energy conservation.

// Time step for t < t_scaling_max
@dt_dyn_epsilon_param  1e-6
// Final simulation time
@t_fin  50000e-6
// Time interval for diagnostic output (e.g., console output)
@dt_dia 100e-6
// Time interval for particle snapshots (e.g., visualization updates)
@dt_out 100e-6

###################### GRAPHICS and OUTPUT ###############################

@SIZE_affichage 10e-2 // Size of the display in the visualization window
// Delay time between two consecutive visualizations to avoid overly fast updates
@t_wait_affichage  1e-1
// Enable or disable graphics (0 = false, 1 = true)
@Graphics 1
// Rotation vector (rot_axe_x, rot_axe_y, rot_axe_z) for 90-degree rotations.
// Default axes:
// - (1,0,0): x-axis right, y-axis into the screen, z-axis up (gravity down).
// - (0,1,0): x-axis out of screen, y-axis up, z-axis left.
// - (0,0,1): x-axis up, y-axis left, z-axis out of screen.
@rot_axe_x  1
@rot_axe_y  0
@rot_axe_z  0

// Scaling for molecule visualization based on vibrational (v_max) and rotational (J_max) levels.
@v_maX 2
@two_J_maX 4

# OUTPUT
#
// Excludes specific manifolds from statistical analysis.
// For example, -1 for dead levels or photoionized states.
// Use a value like -10 to include all manifolds.
@num_manifold_not_studied  -10
// Index of the specific level to study for statistics.
// Note: This is the index (starting from 0) in the file, not the manifold number.
// Use -1 to include all molecules.
@num_niveau_etudie  -1

######### 	Diagonalization ################################
// Indicates whether energy levels are diagonalized or calculated using simple analytical formulas (linear, quadratic, or two-level approximations).
// 0: false (use standard formulas).
// 1: true (perform full diagonalization).
// 2 (or other values): true but first-order diagonalization (block diagonalization by manifold for eigenvalues and eigenvectors).
@is_Levels_Lines_Diagonalized  0

######### 	EXTERNAL FIELDS (SI units: T, T/m, T/m^2, etc.) ################################
// For complex cases, refer to is_Levels_Lines_Diagonalized.

// Currently, only one type of shift (Zeeman or Stark) can be applied to internal states due to a single parameter in the Level file.
// Choose the shift type:
// - 0: magnetic (Zeeman shift based on Delta and C parameters in the Level file).
// - 1: electric (Stark shift).
@type_of_field_for_internal_state_shift  0

// Both electric and magnetic fields can be used for charged particle dynamics (Lorentz forces).
// Example: In a Penning trap, magnetic fields are fully treated (Zeeman + Lorentz force),
// but electric fields affect only the Lorentz force, not Stark shifts.

/*** Field Definition Methods ***/
// - 0: Analytical (default) — Second-order plus nth-order field expansions.
// - 1: Helmholtz/anti-Helmholtz coils (usually for magnetic fields).
// Additional grid-based options (2–4) allow for defining fields via files or precomputed data:
@type_field_read_E  0 // Electric field
@type_field_read_B  1 // Magnetic field

// Grid-based field configurations:
// - 2: Cylindrical 3D map with four columns (r, z, F_r(r,z), F_z(r,z)).
// - 3: Cylindrical 3D map with additional derivatives (d/dr, d/dz, d^2/drdz).
// - 4: 3D field map with six columns (x, y, z, Bx, By, Bz) [TO BE IMPLEMENTED].


## MAGNETIC FIELD definition (if using coils) ##
#
// In case @type_field_read_B = 1 (Helmholtz coils):
// Number of coils along the z-axis, spaced by (N-1/2) * @gap_bobines. Add +r offset if anti-Helmholtz.
@Nb_bobines 5 // Number of coils

// If @is_Helmholtz is true (1), double the coils (add offset +r).
// -1: Anti-Helmholtz configuration.
//  0: No additional coils, only the defined coils are created.
@is_Helmholtz 0

// Spacing between coils (in meters)
@gap_bobines 4e-2

// Current in coils (in amperes). Magnetic field at center: µ0 * I / (2 * r).
// Example: 0.63 mT for 1A with r = 1mm.
@courant_bobines 16000

// Radius of the coils (in meters)
@rayon_bobines 1e-2



#
## MAGNETIC FIELD definition (analytical) ##
#
// Magnetic field along x, y, and z decomposed into components.
// Example for Ox: B_x + grad_B_x * x + grad_grad_B_x * x^2 + Bn * x^n.
// NEVER set these values to 0; use a small value like 1e-10 to preserve quantization axis.

@B_x 1e-10      // Base magnetic field along x-axis
@B_y 2e-10      // Base magnetic field along y-axis
@B_z 0.0        // Base magnetic field along z-axis

// Gradients of the magnetic field (first-order derivatives)
@grad_B_x 0.0
@grad_B_y 0.0
@grad_B_z 0.0

// Second-order gradients of the magnetic field
@grad_grad_B_x 0.0
@grad_grad_B_y 0.0
@grad_grad_B_z 0.0

// Higher-order field terms
@n_value_B 3 // Polynomial degree for nth-order terms
@Bn_x 0.0
@Bn_y 0.0
@Bn_z 0.0



## ELECTRIC FIELD ##
#
// Electric field along x, y, and z decomposed into components.
// Example for Ox: E_x + grad_E_x * x + grad_grad_E_x * x^2 + En * x^n.
// NEVER set these values to 0; use a small value like 1e-10 for stability.

@E_x 0.0        // Base electric field along x-axis
@E_y 0.0        // Base electric field along y-axis
@E_z 0.0        // Base electric field along z-axis

// Gradients of the electric field (first-order derivatives)
@grad_E_x 0.0
@grad_E_y 0.0
@grad_E_z 0.0

// Second-order gradients of the electric field
@grad_grad_E_x 0.0
@grad_grad_E_y 0.0
@grad_grad_E_z 0.0

// Higher-order field terms
@n_value_E 3 // Polynomial degree for nth-order terms
@En_x 0.0
@En_y 0.0
@En_z 0.0


########## LASERS ########################################
#
// Scaling factor for all laser powers
@scale_Power 1.0

// Additive offset to the frequency of all lasers
@Offset_Detuning_cm 0.0

// Scaling factor for the spectral width of all lasers
@scale_Gamma 1.0

// Number of lasers used (this may be fewer than the lasers defined below)
@Nb_laser 1

// Force rates not linked to lasers (0 = false, 1 = true)
@is_forced_rates 1

# First laser (Laser 1, referenced as laser 0 in the program)
@waist_pos_x[0] 0.0
@waist_pos_y[0] 0.0
@waist_pos_z[0] 0.0
@direction_x[0] 1.0
@direction_y[0] 0.0
@direction_z[0] 0.0

// Waist size (in meters)
@waist[0] 3e-3

// Energy levels and parameters for the laser
@Energie_cm[0] 0.92
@Gamma_L_MHz[0] 10.0
@Power[0] 1e-10

// Laser polarization vector in the propagation frame
@Pol_circulaire_left_sp[0] 0.707107
@Pol_circulaire_right_sm[0] 0.707107
@polar_angle_degree[0] 0.0

// Laser type (refer to laser.h for details)
// Examples: CW = 0, Gaussian = 5 (most used), Lorentzian = 6, Comb = 7
@type_laser[0] 5

// Femtosecond comb laser parameters (used only if type_laser = 7)
@nu_offset_MHz[0] 0.0
@nu_repetition_MHz[0] 80.0
@nu_individual_comb_line_MHz[0] 80.0

// Laser interference settings
// -1: This laser does not interfere with others.
@coherent_avec_laser_num[0] -1



#Secons laser (if needed). Laser n°2
@waist_pos_x[1]	0.
@waist_pos_y[1]	0.
@waist_pos_z[1]	0
@direction_x[1]	0.
@direction_y[1]	0.
@direction_z[1]	-1.

@waist[1]	2e-3
@Energie_cm[1] 11732.1813

@Gamma_L_MHz[1]	1
@Power[1]	0.01
//Pol_circulaire_left_sp[1]    -0.7071 and Pol_circulaire_right_sm[1]   0.7071 mean pol pi for a laser along z and B field along x
@Pol_circulaire_left_sp[1]    -0.7071
@Pol_circulaire_right_sm[1]   0.7071
@polar_angle_degree[1]  0
@type_laser[1]  5
@coherent_avec_laser_num[1]  -1


########## OUTPUT PARAMETERS ##############################
#
// Save the parameter file (datacard)?
// Options: 0 = no, 1 = yes, 2 = separate files (one for datacard, one for data without datacard).
@is_DataCard_out 2

// Output scanned parameters to a file (0 = no, 1 = yes)
@is_param_scanned_out 1



########## FILE NAMES ####################################
// Files containing energy levels, transitions, and spectra
@nom_file_Levels Data/n/n_levels.dat
@nom_file_Lines Data/n/n_lines.dat
@nom_file_Laser_Spectrum Data/n/Laser_Spectrum.dat
@nom_file_Laser_Intensity Data/n/Laser_Intensity.dat

// Output files
@nom_sortie_donnees Data/donnee_Mol.dat
@nom_sortie_rate Data/sortie_rate.dat
@nom_sortie_donnees_Data Data/data_card.dat
@nom_fichier_random_gen Data/random_gen.txt

############# Monte Carlo and N-body Algorithm Choices ##############
#
// Verify which algorithms are functional. For now, -1 and 0 are working options.
// Monte Carlo algorithms:
// - Aucun_MC = -1 (does not calculate rates).
// - Kinetic_Monte_Carlo = 0 (standard KMC method).
// - Random_Selection_Method = 1 (selects transitions randomly).
// - First_Reaction_Method = 2 (processes the earliest reaction first).
// - Fast_Rough_Method = 3 (new method, evolves ~0.1/rate_max time steps, meaning N_Mol/10 typically evolve per step).
@Choix_algorithme_Monte_Carlo 0 // Default: Kinetic Monte Carlo

// N-body algorithms:
// - Aucun_N_corps = -1 (only photon recoil).
// - Verlet_acc = 1 (without dipolar forces).
// - Verlet_pot = 2 (with dipolar potential).
// - Verlet_pot_gradient_high_order = 6 (dipolar potential with higher-order gradient calculations in one_body).
// - Yoshida6_acc = 3 (6th-order symplectic algorithm based on Verlet_acc).
// - Yoshida6_pot = 4 (6th-order symplectic algorithm based on Verlet_pot).
// - Boris_Buneman = 7 (magnetic field for charged particles).
@Choix_algorithme_N_corps 1 // Default: Verlet_acc (no dipolar forces)

// Epsilon value for potential position derivative (in meters).
// A variation of 1/100 of the potential with epsilon is generally effective.
// Example: 1e-6 for standard lasers, or 1e-9 if interferences are present.
// Test values of epsilon < 0 for further tuning and cross-verify with epsilon_param.
@choix_epsilon 1e-6 // Default epsilon value

// Seed for random number generator. The same value ensures repeatable sequences.
// Use a negative value to load a new sequence from a file for every run (not commonly used).
@Seed_Init_Random_Number_Generator 2 // Default: reproducible sequence



############# Lists of Scanned or Time-varying Variables #############
#
// Scanning and time-dependent variables are defined in lists.
// The list starts with BEGIN_OF... and ends with END_OF....
// Each scanned variable must be named as SCAN_<VariableName>, followed by:
// - Minimum value, maximum value, number of steps (nbstep).
// - Two booleans: is_scanned (true/false) and is_time (true/false).
// - tau_var: Time constant for exponential changes if is_time is true.
// Nb_steps refers to the number of intervals:
// - 1 means two values (min and max), 2 means three values (min, mid, max), etc.
//
// Value at time t is calculated as:
// val_t0 = min_value + steps * (max_value - min_value) / nbstep
// tau = time rate of change → val = val_t0 * exp^(-t / tau).
//
// If tau is unspecified but is_time is true, Tau_Modif will be used as the default.
// If multiple tables are needed, specify the table number as [table_index].

# Default time constant for exponential modifications
@Tau_Modif 1e-3

# Randomized or ordered scans
@is_Scan_Random false // Default: ordered scanning

# Scanned variables
BEGIN_OF_FITPARAMS
@SCAN_Offset_Detuning_cm -0.01 0.01 20 false false // Offset detuning in cm^-1
@SCAN_scale_Gamma 0.3 0.6 3 false false           // Spectral width scaling
@SCAN_Tau_Modif 0.5e-3 2e-3 2 false false         // Time constant for modifications
@SCAN_E_z 0 1e6 1000 false false                  // Electric field along z-axis
@SCAN_scale_Power 0.1 5.1 10 false false          // Power scaling factor
END_OF_FITPARAMS
