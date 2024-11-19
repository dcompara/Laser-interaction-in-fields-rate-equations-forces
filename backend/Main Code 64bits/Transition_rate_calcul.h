<<<<<<< HEAD
/*
  Name: Transition Rate Calculation Algorithm
  Author: Daniel Comparat
  Date: 16/12/08

  Description:
  This file provides functions to calculate transition rates for molecular systems under laser excitation.
  It handles single-photon and multi-photon transitions, spontaneous emissions, ionization, and rate
  calculations for bound and continuum states.

  Key Features:
  - Includes Gaussian, Lorentzian, and Voigt laser profiles.
  - Supports pseudo-blackbody radiation (BBR) spectra and laser frequency combs.
  - Handles interference and interactions between multiple lasers.

  Notes:
  - Units are in cm^-1 for energy and spectral calculations.
  - Lasers are assumed to be continuous wave (CW). Special handling is required for femtosecond lasers.

  Dependencies:
  - Eigen for matrix manipulations.
  - GSL for random number generation and statistical operations.
*/

#ifndef Transition_rate_SEEN
#define Transition_rate_SEEN

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Eigen> // Matrix manipulation
#include "laser.h"
#include "shift_molecule.h"
#include "sortie_donnees.h"
#include "Kinetic_Monte_Carlo.h"
#include "algorithmes.h"
#include "params.h"

using namespace std;
using namespace Eigen;

/************************************************************************/
/******************** Laser-Molecule Interaction Functions **************/
/************************************************************************/

/**
 * Calculates the spontaneous emission rate for a transition.
 * @param dip_Debye Dipole moment in Debye.
 * @param energie_cm Transition energy in cm^-1.
 * @return Spontaneous emission rate.
 */
double Gamma_spon(const double dip_Debye, const double energie_cm);

/**
 * Gaussian laser profile.
 * @param I_tot Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma FWHM of the Gaussian profile.
 * @return Intensity at the specified detuning.
 */
inline double Gauss(const double I_tot, const double delta, const double Gamma);

/**
 * Flat laser profile.
 * @param I_tot Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma FWHM of the laser.
 * @return Intensity at the specified detuning.
 */
inline double intensity_flat(const double I_tot, const double delta, const double Gamma);

/**
 * Lorentzian laser profile.
 * @param I_tot Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma FWHM of the Lorentzian profile.
 * @return Intensity at the specified detuning.
 */
inline double Lorentz(const double I_tot, const double delta, const double Gamma);

/**
 * Pseudo-Voigt profile: Combines Lorentzian and Gaussian profiles.
 * @param delta Detuning from the transition frequency.
 * @param G_L Lorentzian FWHM.
 * @param G_G Gaussian FWHM.
 * @return Intensity at the specified detuning.
 */
double Pseudo_Voigt(double delta, const double G_L, const double G_G);

/**
 * Comb spectrum: Gaussian laser with individual Lorentzian comb lines.
 * @param I Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param G_L Lorentzian FWHM.
 * @param G_G Gaussian FWHM.
 * @param n_las Laser index.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param params Simulation parameters.
 * @return Intensity at the specified detuning.
 */
double comb_shape(double I, double delta, const double G_L, const double G_G, const int n_las,
                  const double Energy_transition_cm, FitParams &params);

/**
 * Pseudo-blackbody radiation (BBR) spectrum.
 * @param I Total intensity.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param T Temperature in Kelvin.
 * @return Intensity at the specified transition energy.
 */
inline double pseudo_BBR_intensity(const double I, const double Energy_transition_cm, const double T);

/**
 * Local intensity with linewidth convolution.
 * @param I Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma_atomic Atomic linewidth (natural linewidth).
 * @param Gamma_Laser Laser linewidth.
 * @param laser_type Type of laser (CW, pulse, etc.).
 * @param n_las Laser index.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param Energy_transition_laser_cm Laser energy in cm^-1.
 * @param params Simulation parameters.
 * @return Intensity at the specified detuning.
 */
double intensity_Convolution_linewidth(const double I, const double delta, const double Gamma_atomic,
                                       const double Gamma_Laser, const int laser_type, const int n_las,
                                       const double Energy_transition_cm, const double Energy_transition_laser_cm,
                                       FitParams &params);

/**
 * Two-level system rate calculation.
 * @param I_loc Local intensity.
 * @param dip_Debye Dipole moment in Debye.
 * @return Transition rate.
 */
double rate_two_level(const double I_loc, const double dip_Debye);

/**
 * Rate for photo-ionization.
 * @param I_tot Total intensity.
 * @param dip_Debye Dipole moment in Debye.
 * @return Ionization rate.
 */
double rate_ionization(const double I_tot, const double dip_Debye);

/**
 * Field ionization rate for a given molecule under a laser field.
 * @param my_laser Laser configuration.
 * @param Mol Vector of molecules.
 * @param n_mol Molecule index.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @return Field ionization rate.
 */
double rate_field_ionization(const Laser &my_laser, const vector<Molecule> &Mol, const int n_mol,
                             const Field &fieldB, const Field &fieldE);


/**
 * Calculates the excitation rate for a given transition (bound or ionization) caused by a laser (or a class of coherent lasers).
 * Adds the calculated rate to the reaction list only if `is_rate_recorded` is true.
 * @param reaction_list List of reactions (output).
 * @param rate Vector to store calculated rates (output).
 * @param n_las Index of the laser.
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param k Wave vector of the laser.
 * @param Itot_loc Local intensity of the laser.
 * @param Itot Total intensity of the laser.
 * @param dipole_debye Dipole moment in Debye.
 * @param delta Detuning from the transition energy.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param is_bound_transition Indicates if the transition is bound.
 * @param Internal_state_out Final internal state of the transition.
 * @param is_rate_recorded Determines whether the rate is added to the reaction list.
 * @return The calculated excitation rate.
 */
double rate_excitation(vector<type_codage_react>& reaction_list, vector<double>& rate, const int n_las,
                       const vector<Molecule>& Mol, const int n_mol, const Vecteur3D k, const double Itot_loc,
                       const double Itot, const double dipole_debye, const double delta, const double Energy_transition_cm,
                       const Field& fieldB, const Field& fieldE, const int is_bound_transition,
                       const Internal_state Internal_state_out, const bool is_rate_recorded);

/************************************************************************/
/************************ Molecule Transition Functions *****************/
/************************************************************************/

/**
 * Calculates all spontaneous emission rates for a molecule.
 * @return The number of calculated spontaneous emission rates.
 */
int rates_molecule_spon(vector<Internal_state>& Level, vector<type_codage_react>& reaction_list,
                        vector<double>& rate, const Molecule& my_mol, const Field& fieldB, const Field& fieldE,
                        const int num_mol, FitParams& params, SelfAdjointEigenSolver<MatrixXcd>& es, MatrixXcd& H,
                        MatrixXcd& E0_cm, MatrixXcd& Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[]);

/**
 * Calculates rates for a single molecule transitioning between levels under the influence of a specific laser.
 * Includes light shift effects in the dipolar potential.
 */
int rates_single_molecule_laser_level(const int n_las, double dipole, double& delta, double& Gamma_spon_tot,
                                      double sqrt_intensity_loc[], vector<Internal_state>& Level,
                                      Internal_state& Internal_state_in, Internal_state& Internal_state_out,
                                      vector<type_codage_react>& reaction_list, vector<double>& rate,
                                      const vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                                      const Field& fieldE, const Laser& my_laser, const double t,
                                      double& delta_pot_dipolaire, FitParams& params, MatrixXcd d[],
                                      bool is_rate_recorded, int is_bound_transition = 1, const int n_level_in = 0,
                                      const int n_level_out = 1, const double Gamma_in = 0.);

/**
 * Calculates forced rates for a single molecule for purposes such as repumping or Sisyphus cooling.
 */
int forced_rates_single_molecule_laser_level(const int n_las, double dipole, double& delta, double& Gamma_spon_tot,
                                             double sqrt_intensity_loc[], vector<Internal_state>& Level,
                                             Internal_state& Internal_state_in, Internal_state& Internal_state_out,
                                             vector<type_codage_react>& reaction_list, vector<double>& rate,
                                             const vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                                             const Field& fieldE, const Laser& my_laser, const double t,
                                             double& delta_pot_dipolaire, FitParams& params, MatrixXcd d[],
                                             bool is_rate_recorded, int is_bound_transition = 1,
                                             const int n_level_in = 0, const int n_level_out = 1,
                                             const double Gamma_in = 0.);

/**
 * Calculates all transition rates for a molecule under all lasers and adds them to the existing rates.
 */
int rates_molecule(vector<Internal_state>& Level, vector<type_codage_react>& reaction_list, vector<double>& rate,
                   const vector<Molecule>& Mol, const int n_mol, const Field& fieldB, const Field& fieldE,
                   const vector<Laser>& laser, const double t, double& delta_pot_dipolaire, FitParams& params,
                   MatrixXcd& H, MatrixXcd& E0_cm, MatrixXcd& Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[],
                   bool is_rate_recorded = true);

/**
 * Calculates rates considering laser interference effects for a given molecule.
 */
void rates_single_molecule_lasers_interference_level(vector<Internal_state>& Level, vector<type_codage_react>& reaction_list,
                                                     vector<double>& rate, const vector<Molecule>& Mol,
                                                     const int n_mol, const Field& fieldB, const Field& fieldE,
                                                     const vector<Laser>& laser, double sqrt_intensity_loc[],
                                                     Internal_state& Internal_state_in, Internal_state& Internal_state_out,
                                                     double delta[], double dipole_debye, double& delta_pot_dipolaire,
                                                     double Gamma_spon_tot, double t, bool is_bound_transition,
                                                     bool is_rate_recorded);

/**
 * Calculates the decay rate from a given energy level.
 */
double Gamma_Level_from_diagonalized_dipole(vector<Internal_state>& Level, MatrixXcd d[], const int i);

/**
 * Copies all rates except for a specific molecule.
 */
int copie_rates_molecules(vector<type_codage_react>& reaction_list, vector<double>& rate, const int numero_mol = -1,
                          const int nombre_old_rate = 0);

/**
 * Calculates transition rates for all molecules or for a specific molecule if `numero_mol` is provided.
 */
int calcul_rates_molecules(vector<Internal_state>& Level, MC_algorithmes Algorithme_MC,
                           vector<type_codage_react>& reaction_list, vector<double>& rate, const vector<Molecule>& Mol,
                           const Field& fieldB, const Field& fieldE, const vector<Laser>& laser, const double t,
                           const int numero_mol, const int N_Mol, FitParams& params, MatrixXcd& H, MatrixXcd& E0_cm,
                           MatrixXcd& Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[]);

/************************************************************************/
/************************ Reaction Handling Functions *******************/
/************************************************************************/

/**
 * Scans all rates and executes a reaction if the probability matches.
 */
int do_reaction_FastRoughMethod(const MC_algorithmes Algorithme_MC, const gsl_rng* r, const vector<double>& rate,
                                const vector<type_codage_react>& reaction_list, vector<Molecule>& Mol, const int n_reac,
                                const vector<Laser>& laser, const double dt_KMC, ofstream& file_rate,
                                int& number_photons, FitParams& params);

/**
 * Executes a reaction but does not update potentials. Potentials must be updated separately.
 * @return The index of the affected molecule.
 */
int do_reaction(const MC_algorithmes Algorithme_MC, const gsl_rng* r, const vector<double>& rate,
                const vector<type_codage_react>& reaction_list, vector<Molecule>& Mol, const int n_reac,
                const vector<Laser>& my_laser, const double dt_KMC, ofstream& file_rate, bool first_call,
                int& number_photons, FitParams& params);

/**
 * Generates a random unit vector for spontaneous emission based on quantization axis or polarization vector.
 */
Vecteur3D get_unit_vector_spontaneous_emission(const gsl_rng* r, complex<double> e_pol_dipole_transition[3],
                                               Vecteur3D quantization_axis, int delta_M, FitParams& params);

#endif
=======
/*
  Name: Transition Rate Calculation Algorithm
  Author: Daniel Comparat
  Date: 16/12/08

  Description:
  This file provides functions to calculate transition rates for molecular systems under laser excitation.
  It handles single-photon and multi-photon transitions, spontaneous emissions, ionization, and rate
  calculations for bound and continuum states.

  Key Features:
  - Includes Gaussian, Lorentzian, and Voigt laser profiles.
  - Supports pseudo-blackbody radiation (BBR) spectra and laser frequency combs.
  - Handles interference and interactions between multiple lasers.

  Notes:
  - Units are in cm^-1 for energy and spectral calculations.
  - Lasers are assumed to be continuous wave (CW). Special handling is required for femtosecond lasers.

  Dependencies:
  - Eigen for matrix manipulations.
  - GSL for random number generation and statistical operations.
*/

#ifndef Transition_rate_SEEN
#define Transition_rate_SEEN

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Eigen> // Matrix manipulation
#include "laser.h"
#include "shift_molecule.h"
#include "sortie_donnees.h"
#include "Kinetic_Monte_Carlo.h"
#include "algorithmes.h"
#include "params.h"

using namespace std;
using namespace Eigen;

/************************************************************************/
/******************** Laser-Molecule Interaction Functions **************/
/************************************************************************/

/**
 * Calculates the spontaneous emission rate for a transition.
 * @param dip_Debye Dipole moment in Debye.
 * @param energie_cm Transition energy in cm^-1.
 * @return Spontaneous emission rate.
 */
double Gamma_spon(const double dip_Debye, const double energie_cm);

/**
 * Gaussian laser profile.
 * @param I_tot Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma FWHM of the Gaussian profile.
 * @return Intensity at the specified detuning.
 */
inline double Gauss(const double I_tot, const double delta, const double Gamma);

/**
 * Flat laser profile.
 * @param I_tot Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma FWHM of the laser.
 * @return Intensity at the specified detuning.
 */
inline double intensity_flat(const double I_tot, const double delta, const double Gamma);

/**
 * Lorentzian laser profile.
 * @param I_tot Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma FWHM of the Lorentzian profile.
 * @return Intensity at the specified detuning.
 */
inline double Lorentz(const double I_tot, const double delta, const double Gamma);

/**
 * Pseudo-Voigt profile: Combines Lorentzian and Gaussian profiles.
 * @param delta Detuning from the transition frequency.
 * @param G_L Lorentzian FWHM.
 * @param G_G Gaussian FWHM.
 * @return Intensity at the specified detuning.
 */
double Pseudo_Voigt(double delta, const double G_L, const double G_G);

/**
 * Comb spectrum: Gaussian laser with individual Lorentzian comb lines.
 * @param I Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param G_L Lorentzian FWHM.
 * @param G_G Gaussian FWHM.
 * @param n_las Laser index.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param params Simulation parameters.
 * @return Intensity at the specified detuning.
 */
double comb_shape(double I, double delta, const double G_L, const double G_G, const int n_las,
                  const double Energy_transition_cm, FitParams &params);

/**
 * Pseudo-blackbody radiation (BBR) spectrum.
 * @param I Total intensity.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param T Temperature in Kelvin.
 * @return Intensity at the specified transition energy.
 */
inline double pseudo_BBR_intensity(const double I, const double Energy_transition_cm, const double T);

/**
 * Local intensity with linewidth convolution.
 * @param I Total intensity.
 * @param delta Detuning from the transition frequency.
 * @param Gamma_atomic Atomic linewidth (natural linewidth).
 * @param Gamma_Laser Laser linewidth.
 * @param laser_type Type of laser (CW, pulse, etc.).
 * @param n_las Laser index.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param Energy_transition_laser_cm Laser energy in cm^-1.
 * @param params Simulation parameters.
 * @return Intensity at the specified detuning.
 */
double intensity_Convolution_linewidth(const double I, const double delta, const double Gamma_atomic,
                                       const double Gamma_Laser, const int laser_type, const int n_las,
                                       const double Energy_transition_cm, const double Energy_transition_laser_cm,
                                       FitParams &params);

/**
 * Two-level system rate calculation.
 * @param I_loc Local intensity.
 * @param dip_Debye Dipole moment in Debye.
 * @return Transition rate.
 */
double rate_two_level(const double I_loc, const double dip_Debye);

/**
 * Rate for photo-ionization.
 * @param I_tot Total intensity.
 * @param dip_Debye Dipole moment in Debye.
 * @return Ionization rate.
 */
double rate_ionization(const double I_tot, const double dip_Debye);

/**
 * Field ionization rate for a given molecule under a laser field.
 * @param my_laser Laser configuration.
 * @param Mol Vector of molecules.
 * @param n_mol Molecule index.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @return Field ionization rate.
 */
double rate_field_ionization(const Laser &my_laser, const vector<Molecule> &Mol, const int n_mol,
                             const Field &fieldB, const Field &fieldE);


/**
 * Calculates the excitation rate for a given transition (bound or ionization) caused by a laser (or a class of coherent lasers).
 * Adds the calculated rate to the reaction list only if `is_rate_recorded` is true.
 * @param reaction_list List of reactions (output).
 * @param rate Vector to store calculated rates (output).
 * @param n_las Index of the laser.
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param k Wave vector of the laser.
 * @param Itot_loc Local intensity of the laser.
 * @param Itot Total intensity of the laser.
 * @param dipole_debye Dipole moment in Debye.
 * @param delta Detuning from the transition energy.
 * @param Energy_transition_cm Transition energy in cm^-1.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param is_bound_transition Indicates if the transition is bound.
 * @param Internal_state_out Final internal state of the transition.
 * @param is_rate_recorded Determines whether the rate is added to the reaction list.
 * @return The calculated excitation rate.
 */
double rate_excitation(vector<type_codage_react>& reaction_list, vector<double>& rate, const int n_las,
                       const vector<Molecule>& Mol, const int n_mol, const Vecteur3D k, const double Itot_loc,
                       const double Itot, const double dipole_debye, const double delta, const double Energy_transition_cm,
                       const Field& fieldB, const Field& fieldE, const int is_bound_transition,
                       const Internal_state Internal_state_out, const bool is_rate_recorded);

/************************************************************************/
/************************ Molecule Transition Functions *****************/
/************************************************************************/

/**
 * Calculates all spontaneous emission rates for a molecule.
 * @return The number of calculated spontaneous emission rates.
 */
int rates_molecule_spon(vector<Internal_state>& Level, vector<type_codage_react>& reaction_list,
                        vector<double>& rate, const Molecule& my_mol, const Field& fieldB, const Field& fieldE,
                        const int num_mol, FitParams& params, SelfAdjointEigenSolver<MatrixXcd>& es, MatrixXcd& H,
                        MatrixXcd& E0_cm, MatrixXcd& Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[]);

/**
 * Calculates rates for a single molecule transitioning between levels under the influence of a specific laser.
 * Includes light shift effects in the dipolar potential.
 */
int rates_single_molecule_laser_level(const int n_las, double dipole, double& delta, double& Gamma_spon_tot,
                                      double sqrt_intensity_loc[], vector<Internal_state>& Level,
                                      Internal_state& Internal_state_in, Internal_state& Internal_state_out,
                                      vector<type_codage_react>& reaction_list, vector<double>& rate,
                                      const vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                                      const Field& fieldE, const Laser& my_laser, const double t,
                                      double& delta_pot_dipolaire, FitParams& params, MatrixXcd d[],
                                      bool is_rate_recorded, int is_bound_transition = 1, const int n_level_in = 0,
                                      const int n_level_out = 1, const double Gamma_in = 0.);

/**
 * Calculates forced rates for a single molecule for purposes such as repumping or Sisyphus cooling.
 */
int forced_rates_single_molecule_laser_level(const int n_las, double dipole, double& delta, double& Gamma_spon_tot,
                                             double sqrt_intensity_loc[], vector<Internal_state>& Level,
                                             Internal_state& Internal_state_in, Internal_state& Internal_state_out,
                                             vector<type_codage_react>& reaction_list, vector<double>& rate,
                                             const vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                                             const Field& fieldE, const Laser& my_laser, const double t,
                                             double& delta_pot_dipolaire, FitParams& params, MatrixXcd d[],
                                             bool is_rate_recorded, int is_bound_transition = 1,
                                             const int n_level_in = 0, const int n_level_out = 1,
                                             const double Gamma_in = 0.);

/**
 * Calculates all transition rates for a molecule under all lasers and adds them to the existing rates.
 */
int rates_molecule(vector<Internal_state>& Level, vector<type_codage_react>& reaction_list, vector<double>& rate,
                   const vector<Molecule>& Mol, const int n_mol, const Field& fieldB, const Field& fieldE,
                   const vector<Laser>& laser, const double t, double& delta_pot_dipolaire, FitParams& params,
                   MatrixXcd& H, MatrixXcd& E0_cm, MatrixXcd& Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[],
                   bool is_rate_recorded = true);

/**
 * Calculates rates considering laser interference effects for a given molecule.
 */
void rates_single_molecule_lasers_interference_level(vector<Internal_state>& Level, vector<type_codage_react>& reaction_list,
                                                     vector<double>& rate, const vector<Molecule>& Mol,
                                                     const int n_mol, const Field& fieldB, const Field& fieldE,
                                                     const vector<Laser>& laser, double sqrt_intensity_loc[],
                                                     Internal_state& Internal_state_in, Internal_state& Internal_state_out,
                                                     double delta[], double dipole_debye, double& delta_pot_dipolaire,
                                                     double Gamma_spon_tot, double t, bool is_bound_transition,
                                                     bool is_rate_recorded);

/**
 * Calculates the decay rate from a given energy level.
 */
double Gamma_Level_from_diagonalized_dipole(vector<Internal_state>& Level, MatrixXcd d[], const int i);

/**
 * Copies all rates except for a specific molecule.
 */
int copie_rates_molecules(vector<type_codage_react>& reaction_list, vector<double>& rate, const int numero_mol = -1,
                          const int nombre_old_rate = 0);

/**
 * Calculates transition rates for all molecules or for a specific molecule if `numero_mol` is provided.
 */
int calcul_rates_molecules(vector<Internal_state>& Level, MC_algorithmes Algorithme_MC,
                           vector<type_codage_react>& reaction_list, vector<double>& rate, const vector<Molecule>& Mol,
                           const Field& fieldB, const Field& fieldE, const vector<Laser>& laser, const double t,
                           const int numero_mol, const int N_Mol, FitParams& params, MatrixXcd& H, MatrixXcd& E0_cm,
                           MatrixXcd& Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[]);

/************************************************************************/
/************************ Reaction Handling Functions *******************/
/************************************************************************/

/**
 * Scans all rates and executes a reaction if the probability matches.
 */
int do_reaction_FastRoughMethod(const MC_algorithmes Algorithme_MC, const gsl_rng* r, const vector<double>& rate,
                                const vector<type_codage_react>& reaction_list, vector<Molecule>& Mol, const int n_reac,
                                const vector<Laser>& laser, const double dt_KMC, ofstream& file_rate,
                                int& number_photons, FitParams& params);

/**
 * Executes a reaction but does not update potentials. Potentials must be updated separately.
 * @return The index of the affected molecule.
 */
int do_reaction(const MC_algorithmes Algorithme_MC, const gsl_rng* r, const vector<double>& rate,
                const vector<type_codage_react>& reaction_list, vector<Molecule>& Mol, const int n_reac,
                const vector<Laser>& my_laser, const double dt_KMC, ofstream& file_rate, bool first_call,
                int& number_photons, FitParams& params);

/**
 * Generates a random unit vector for spontaneous emission based on quantization axis or polarization vector.
 */
Vecteur3D get_unit_vector_spontaneous_emission(const gsl_rng* r, complex<double> e_pol_dipole_transition[3],
                                               Vecteur3D quantization_axis, int delta_M, FitParams& params);

#endif
>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
