// Initialization of levels, transitions, line forces, and lifetimes
// Laser parameters, molecular positions, and velocities...

#ifndef Initialisation_programme_SEEN
#define Initialisation_programme_SEEN

#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>                   // For random number generator
#include <gsl/gsl_randist.h>               // For Gaussian random distribution
#include <gsl/gsl_sf_lambert.h>            // For Lambert W function
#include <gsl/gsl_cdf.h>                   // For cumulative distribution functions

#include "Vecteur3D.h"
#include "Molecule.h"                      // Molecule class
#include "Laser.h"                         // Laser class
#include "Field.h"                         // Field (Force) class
#include "Transitions_initialisation.h"
#include "shift_molecule.h"
#include "params.h"                        // Parameter handling
#include <sstream>                         // For ostringstream

/************************************************************************/
/********************* PROGRAM INITIALIZATION ***************************/
/************************************************************************/

/**
 * Initializes levels, transitions, line forces, and lifetimes.
 * @param nom_file_Levels Path to the file containing level information.
 * @param nom_file_Lines Path to the file containing line information.
 * @param Level Vector to store internal state levels.
 * @param params Fit parameters for the simulation.
 * @return Number of levels initialized.
 */
int initialisation_trans_mol(const char *nom_file_Levels, const char *nom_file_Lines,
                             vector<Internal_state>& Level, FitParams& params);

/**
 * Initializes the state of the random number generator.
 * @param r Pointer to the random generator.
 * @param nombre_seed Seed value for random number generation.
 * @param nom_file Path to the file containing random state configuration.
 */
void set_random_generator(gsl_rng* r, int nombre_seed, const char* nom_file);

/**
 * Initializes molecule positions and velocities with a Gaussian distribution.
 * @param r Pointer to the random generator.
 * @param Mol Vector of molecules to initialize.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param Nb_of_type_of_different_Mol Number of different molecule types.
 * @param params Fit parameters for the simulation.
 * @param data DataCards object for parameter input.
 */
void Init_Molecule(const gsl_rng* r, vector<Molecule>& Mol, const Field fieldB,
                   const Field fieldE, const int Nb_of_type_of_different_Mol,
                   FitParams& params, DataCards data);

/**
 * Initializes the lasers used in the simulation.
 * @param laser Vector of lasers to initialize.
 * @param Nb_lasers Number of lasers to configure.
 * @param params Fit parameters for the simulation.
 * @param nom_file_spectrum Path to the file containing spectrum data.
 * @param nom_file_intensity Path to the file containing intensity data.
 */
void Init_Laser(vector<Laser>& laser, const int Nb_lasers, FitParams& params,
                const char* nom_file_spectrum, const char* nom_file_intensity);

/**
 * Applies time-dependent modifications to simulation parameters.
 * @param t Current time in the simulation.
 * @param params Parameters to modify.
 */
void Modif_Param(const double t, FitParams& params);

/************************************************************************/
/********************* FIELD INITIALIZATION *****************************/
/************************************************************************/

/**
 * Initializes fields such as magnetic and electric fields.
 * By default, fields are initialized to zero if no files are provided.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param params Fit parameters for the simulation.
 * @param nom_Magn_file Path to the file containing magnetic field data (optional).
 * @param nom_Elec_file Path to the file containing electric field data (optional).
 */
void Init_Field(Field& fieldB, Field& fieldE, FitParams& params,
                const char* nom_Magn_file = nullptr, const char* nom_Elec_file = nullptr);

#endif  // Initialisation_programme_SEEN
