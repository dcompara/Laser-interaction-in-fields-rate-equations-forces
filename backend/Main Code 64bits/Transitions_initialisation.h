<<<<<<< HEAD
/*
  Name: Transition Initialization
  Description:
  This header file contains functions to initialize energy levels, transitions, lifetimes, and populations for molecules.
  It handles:
  - Reading energy levels and transitions from input files.
  - Initializing Zeeman and dipole matrices for diagonalization.
  - Setting initial populations and lifetimes for molecular states.

  Notes:
  - Units for energies and transitions are typically in cm^-1.
  - Includes functionality to compute Wigner 3-j coefficients for angular momentum coupling.
  - Supports parameter reading and matrix initialization for advanced simulations.

  Dependencies:
  - GSL for random number generation and Wigner 3-j calculations.
  - Eigen for matrix manipulations.
*/

#ifndef Transition_init_SEEN
#define Transition_init_SEEN

#include <iostream>
#include <vector>
#include <Eigen/Eigen> // For Hamiltonian matrices
#include <gsl/gsl_rng.h> // For random number generation
#include <gsl/gsl_randist.h> // For Gaussian random generation
#include <gsl/gsl_sf_coupling.h> // For Wigner 3-j coefficients

#include "Internal_state.h"
#include "Molecule.h"
#include "params.h"
#include "datacards.h" // For reading parameter files

using namespace std;
using namespace Eigen;

/************************************************************************/
/************************ Transition Initialization *********************/
/************************************************************************/

/**
 * Reads energy, Zeeman, and dipole matrices for diagonalization.
 * @param E0_cm Matrix to store the zero-field energy levels in cm^-1.
 * @param Zeeman_cm_B Matrix to store the Zeeman energy shifts.
 * @param d0 Array of matrices for dipole moment components.
 */
void Read_Energy_Zeeman_dipole_for_Diagonalization(MatrixXcd& E0_cm, MatrixXcd& Zeeman_cm_B, MatrixXd d0[]);

/**
 * Reads the energy levels from a file and initializes the `Level` vector.
 * @param nom_file Path to the file containing energy level data.
 * @param Level Vector to store the internal states (output).
 * @param params Simulation parameters.
 * @return The number of energy levels read from the file.
 */
int Level_List(const char* nom_file, vector<Internal_state>& Level, FitParams& params);

/**
 * Reads the transition lines from a file and initializes the `Level` vector.
 * @param nom_file Path to the file containing transition line data.
 * @param Level Vector to store the internal states (output).
 * @param params Simulation parameters.
 * @return The number of transitions read from the file.
 */
int Line_List(const char* nom_file, vector<Internal_state>& Level, FitParams& params);

/**
 * Initializes the lifetimes (Gamma) for each level.
 * @param Level Vector of internal states to initialize (output).
 * @param params Simulation parameters.
 */
void initialisation_Gamma(vector<Internal_state>& Level, FitParams& params);

/**
 * Initializes the initial population probabilities for molecules.
 * Populations are randomly distributed across levels using a provided random generator.
 * @param r Random number generator.
 * @param Mol Vector of molecules (output).
 * @param Nb_molecules Number of molecules to initialize.
 * @param Level Vector of internal states representing the energy levels.
 */
void initialisation_proba(const gsl_rng* r, vector<Molecule>& Mol, const int Nb_molecules, const vector<Internal_state>& Level);

#endif

=======
/*
  Name: Transition Initialization
  Description:
  This header file contains functions to initialize energy levels, transitions, lifetimes, and populations for molecules.
  It handles:
  - Reading energy levels and transitions from input files.
  - Initializing Zeeman and dipole matrices for diagonalization.
  - Setting initial populations and lifetimes for molecular states.

  Notes:
  - Units for energies and transitions are typically in cm^-1.
  - Includes functionality to compute Wigner 3-j coefficients for angular momentum coupling.
  - Supports parameter reading and matrix initialization for advanced simulations.

  Dependencies:
  - GSL for random number generation and Wigner 3-j calculations.
  - Eigen for matrix manipulations.
*/

#ifndef Transition_init_SEEN
#define Transition_init_SEEN

#include <iostream>
#include <vector>
#include <Eigen/Eigen> // For Hamiltonian matrices
#include <gsl/gsl_rng.h> // For random number generation
#include <gsl/gsl_randist.h> // For Gaussian random generation
#include <gsl/gsl_sf_coupling.h> // For Wigner 3-j coefficients

#include "Internal_state.h"
#include "Molecule.h"
#include "params.h"
#include "datacards.h" // For reading parameter files

using namespace std;
using namespace Eigen;

/************************************************************************/
/************************ Transition Initialization *********************/
/************************************************************************/

/**
 * Reads energy, Zeeman, and dipole matrices for diagonalization.
 * @param E0_cm Matrix to store the zero-field energy levels in cm^-1.
 * @param Zeeman_cm_B Matrix to store the Zeeman energy shifts.
 * @param d0 Array of matrices for dipole moment components.
 */
void Read_Energy_Zeeman_dipole_for_Diagonalization(MatrixXcd& E0_cm, MatrixXcd& Zeeman_cm_B, MatrixXd d0[]);

/**
 * Reads the energy levels from a file and initializes the `Level` vector.
 * @param nom_file Path to the file containing energy level data.
 * @param Level Vector to store the internal states (output).
 * @param params Simulation parameters.
 * @return The number of energy levels read from the file.
 */
int Level_List(const char* nom_file, vector<Internal_state>& Level, FitParams& params);

/**
 * Reads the transition lines from a file and initializes the `Level` vector.
 * @param nom_file Path to the file containing transition line data.
 * @param Level Vector to store the internal states (output).
 * @param params Simulation parameters.
 * @return The number of transitions read from the file.
 */
int Line_List(const char* nom_file, vector<Internal_state>& Level, FitParams& params);

/**
 * Initializes the lifetimes (Gamma) for each level.
 * @param Level Vector of internal states to initialize (output).
 * @param params Simulation parameters.
 */
void initialisation_Gamma(vector<Internal_state>& Level, FitParams& params);

/**
 * Initializes the initial population probabilities for molecules.
 * Populations are randomly distributed across levels using a provided random generator.
 * @param r Random number generator.
 * @param Mol Vector of molecules (output).
 * @param Nb_molecules Number of molecules to initialize.
 * @param Level Vector of internal states representing the energy levels.
 */
void initialisation_proba(const gsl_rng* r, vector<Molecule>& Mol, const int Nb_molecules, const vector<Internal_state>& Level);

#endif

>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
