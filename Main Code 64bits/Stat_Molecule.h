/*
  Name: Molecular Statistics
  Description:
  This file defines the `Stat` class and functions for calculating statistical parameters
  for a collection of molecules, such as:
    - Mean and standard deviation of velocities and positions.
    - Energy and temperature (various definitions based on velocity distributions).
    - Effective temperature derived from the velocity norm ordering.

  Key Features:
  - Computes statistics for individual and grouped molecular states.
  - Handles energy and trajectory statistics.
  - Useful for analyzing molecular dynamics simulations.

  Notes:
  - Several temperature definitions are supported:
      1. `T_cin`: \( \sigma(v - \langle v \rangle)^2 / k_Boltzmann \)
      2. Based on the standard deviation of velocity (\( \sigma \)).
      3. Effective temperature using the ordered velocity norm (\( 1/2 \, m \, v_i^2 \)).
*/

#ifndef Stat_Mol_SEEN
#define Stat_Mol_SEEN

#include <iostream>
#include <vector>
#include "Molecule.h"
#include "one_body.h"
#include "algorithmes.h" // For sorting (qsort)
#include "params.h"
#include <gsl/gsl_statistics.h> // For statistical calculations (mean, variance, etc.)

using namespace std;

/************************************************************************/
/****************** Molecular Statistics Data Structure *****************/
/************************************************************************/

/**
 * Represents statistical parameters for a collection of molecules.
 * Includes averages, standard deviations, energy, and temperature calculations.
 */
class Stat {
public:
    int population;       // Population of the studied state
    Vecteur3D mean_v;     // Mean velocity
    Vecteur3D mean_pos;   // Mean position
    Vecteur3D sigma_v;    // Standard deviation of velocity
    Vecteur3D sigma_pos;  // Standard deviation of position

    double E_pot;         // Potential energy: \( \sum_i [q V_{elec}(r_i) - m g z_i + \hbar (E_{cm} - E_{cm0}) + 0.5 \sum_{j \neq i} \frac{q_i q_j}{4 \pi \epsilon_0 r_{ij}}] \)
    double E_cin;         // Kinetic energy: \( \frac{1}{2} m (v_x^2 + v_y^2 + v_z^2) \)

    Vecteur3D Temp_cin;   // Temperature based on \( \frac{1}{2} m \langle v - \langle v \rangle \rangle^2 \)
    Vecteur3D Temp_1D_50; // Temperature for the best 50% particles in 1D
    double Temp_3D_50;    // Effective 3D temperature for the best 50% particles

    // Constructor
    Stat() : population(0),
             mean_v(Vecteur3D(0., 0., 0.)),
             mean_pos(Vecteur3D(0., 0., 0.)),
             sigma_v(Vecteur3D(0., 0., 0.)),
             sigma_pos(Vecteur3D(0., 0., 0.)),
             E_pot(0.), E_cin(0.),
             Temp_cin(Vecteur3D(0., 0., 0.)),
             Temp_1D_50(Vecteur3D(0., 0., 0.)),
             Temp_3D_50(0.) {}
};

/************************************************************************/
/****************** Molecular Statistics Functions **********************/
/************************************************************************/

/**
 * Calculates statistics for a collection of molecules.
 * Computes mean, standard deviation, energy, and temperature for each internal state.
 * @param Mol Vector of molecules.
 * @param stat_Mol Vector to store statistics for each internal state.
 * @param Level Vector of internal states.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param Nb_Mol Total number of molecules.
 * @param Nb_state Total number of states.
 * @param params Simulation parameters.
 */
void stat_molecule(const vector<Molecule> &Mol, vector<Stat> &stat_Mol,
                   const vector<Internal_state> &Level, const Field &fieldB,
                   const Field &fieldE, const int Nb_Mol, const int Nb_state,
                   FitParams &params);

/**
 * Calculates statistics for a single internal state.
 * Returns the number of molecules in the studied state.
 * @param Mol Vector of molecules.
 * @param stat Structure to store statistics for the studied state.
 * @param Level Vector of internal states.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param num_niveau_etudie Index of the studied state.
 * @param Nb_Mol Total number of molecules.
 * @param params Simulation parameters.
 * @return Number of molecules in the studied state.
 */
int stat_molecule_un_niveau(const vector<Molecule> &Mol, Stat &stat,
                            const vector<Internal_state> &Level,
                            const Field &fieldB, const Field &fieldE,
                            const int num_niveau_etudie, const int Nb_Mol,
                            FitParams &params);

/**
 * Calculates statistics for trajectories.
 * Computes minimum, mean (sum), and maximum distances between molecules.
 * @param Mol Vector of molecules.
 * @param dist_min Minimum distance (output).
 * @param dist_sum Sum of distances (output).
 * @param dist_max Maximum distance (output).
 * @param nb_call Number of calls to the function (used for statistics).
 * @param Nb_Mol Total number of molecules (default: 1).
 */
void stat_traj_molecule(const vector<Molecule> &Mol, double &dist_min,
                        double &dist_sum, double &dist_max, int &nb_call,
                        const int Nb_Mol = 1);

/**
 * Calculates statistics for molecules in a specific state.
 * Populates a list of molecules in the studied state and computes statistics.
 * @param Mol Vector of molecules.
 * @param liste_Adresses_Mol_dans_niveau List of molecules in the studied state (output).
 * @param stat Structure to store statistics for the studied state.
 * @param Level Vector of internal states.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param Nb_Mol Total number of molecules.
 * @param params Simulation parameters.
 * @return Number of molecules in the studied state.
 */
int stat_molecule_form_list(const vector<Molecule> &Mol,
                            vector<Molecule> &liste_Adresses_Mol_dans_niveau,
                            Stat &stat, const vector<Internal_state> &Level,
                            const Field &fieldB, const Field &fieldE,
                            const int Nb_Mol, FitParams &params);

#endif
