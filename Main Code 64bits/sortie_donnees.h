/*
  Name: Data Output Functions
  Description:
  This header file contains functions for saving and exporting simulation data,
  including molecular states, transition rates, and laser properties. These functions
  support detailed statistical analysis and post-processing.

  Key Features:
  - Save random number generator states to ensure reproducibility.
  - Output molecular populations, transition rates, and laser spectra/intensity.
  - Support for time-dependent outputs for dynamic simulations.

  Dependencies:
  - Includes classes and utilities for molecules, lasers, fields, parameters, and statistics.

  Notes:
  - Functions are designed to output data in a structured format, compatible with external analysis tools.
*/

#ifndef Sortie_Donnee_SEEN
#define Sortie_Donnee_SEEN

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <gsl/gsl_rng.h> // For random generator
#include <gsl/gsl_randist.h> // For Gaussian random generator
#include "Molecule.h" // Molecule class
#include "Laser.h" // Laser class
#include "Field.h" // Field class
#include "Stat_Molecule.h" // Statistics on molecules
#include "Transition_rate_calcul.h" // Transition rate calculations
#include "params.h" // Parameter handling
#include "datacards.h" // DataCards for parameter configuration
#include "diagonalization.h" // Energy level diagonalization

using namespace std;

/**
 * Saves the state of the random number generator to a file.
 * Ensures reproducibility of simulations by saving the generator's state.
 * @param r Pointer to the GSL random generator.
 * @param nom_file File path to save the state.
 */
void save_random_generator(gsl_rng *r, const char *nom_file);

/**
 * Example data output for molecular states and transitions.
 * Outputs data such as molecule properties, field configurations, and laser interactions.
 */
void Sortie_donnee_example(ofstream &file_out, vector<Molecule> &Mol,
                           vector<Internal_state> &Levels, const Field &fieldB,
                           const Field &fieldE, const vector<Laser> &laser,
                           const double t, const int nb_mol, FitParams &params,
                           DataCards &data, const int number_photons);

/**
 * Outputs molecular states and transitions.
 * Similar to Sortie_donnee_example but supports customized data export.
 */
void Sortie_donnee(ofstream &file_out, vector<Molecule> &Mol,
                   vector<Internal_state> &Levels, const Field &fieldB,
                   const Field &fieldE, const vector<Laser> &laser,
                   const double t, const int nb_mol, FitParams &params,
                   DataCards &data, const int number_photons);

/**
 * Outputs molecular populations in states (v, J).
 * Data is structured by time for analysis.
 */
void Sortie_donnee_pop_vJ(ofstream &file_out, const vector<Molecule> &Mol,
                          const int nb_Mol, const double t, const int NX,
                          const int N_two_JX, FitParams &params);

/**
 * Outputs molecular populations in states (N, J, 2M).
 * Data is structured by time for analysis.
 */
void Sortie_donnee_pop_NJ2M(ofstream &file_out, const vector<Molecule> &Mol,
                            const int nb_Mol, const double t, const int NX,
                            const int N_two_JX, FitParams &params);

/**
 * Outputs molecular populations in states (N).
 * Data is structured by time for analysis.
 */
void Sortie_donnee_pop_N(ofstream &file_out, const vector<Molecule> &Mol,
                         const int nb_Mol, const double t, const int NX,
                         const int N_two_JX, FitParams &params);

/**
 * Outputs populations in state (v) over time.
 */
void Sortie_donnee_pop_v(ofstream &file_out, const vector<Molecule> &Mol,
                         const int nb_Mol, const double t, const int NX,
                         FitParams &params, int number_photons);

/**
 * Outputs an example of transition rate data.
 * Provides rates, reaction lists, and molecular states for analysis.
 */
void Sortie_rate_example(ofstream &file_rate, const vector<double> &rate,
                         vector<Internal_state> &Level,
                         const vector<type_codage_react> &reaction_list,
                         const vector<Molecule> &Mol, const Field &fieldB,
                         const Field &fieldE, const vector<Laser> &laser,
                         const int N_Mol, const double t, FitParams &params);

/**
 * Outputs transition rate data for the simulation.
 */
void Sortie_rate(ofstream &file_rate, const vector<double> &rate,
                 vector<Internal_state> &Level,
                 const vector<type_codage_react> &reaction_list,
                 const vector<Molecule> &Mol, const Field &fieldB,
                 const Field &fieldE, const vector<Laser> &laser,
                 const int N_Mol, const double t, FitParams &params);

/**
 * Outputs the spectrum of a laser.
 * @param file_out Output stream for saving the spectrum.
 * @param laser Vector of lasers in the simulation.
 * @param params Fit parameters.
 * @param num_laser Laser index (default: 0).
 */
void Sortie_laser_spectrum(ofstream &file_out, const vector<Laser> &laser,
                           FitParams &params, int num_laser = 0);

/**
 * Outputs the time-dependent intensity of a laser.
 * @param file_out Output stream for saving intensity data.
 * @param laser Vector of lasers in the simulation.
 * @param params Fit parameters.
 * @param num_laser Laser index (default: 0).
 */
void Sortie_laser_intensity(ofstream &file_out, const vector<Laser> &laser,
                            FitParams &params, int num_laser = 0);

#endif
