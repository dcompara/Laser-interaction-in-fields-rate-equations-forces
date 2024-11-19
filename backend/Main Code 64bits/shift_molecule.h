/*
  Name: Energy Shift Calculation for Molecules
  Description:
  This file defines functions to calculate and update energy shifts for molecules
  due to various interactions such as Zeeman, Stark, and dipolar effects. These
  calculations are influenced by external magnetic (B) and electric (E) fields, lasers,
  and other parameters. The file also provides functionality to compute accelerations
  and potentials for molecules based on these shifts.

  Notes:
  - Current implementation focuses on Zeeman effects but can be extended for Stark
    and dipolar interactions.
  - Time dependency is considered in the calculations.
  - Dipolar potential is excluded in certain updates to prevent detuning accumulation.
*/

#ifndef Shift_molecule_SEEN
#define Shift_molecule_SEEN

#include "Molecule.h"
#include "Laser.h"
#include "Field.h"
#include "Initialisation_programme.h"     // For field parameters
#include "params.h"
#include <iostream>

using namespace std;

/**
 * Calculates the potential energy of a molecule's internal state
 * (Zeeman, Stark, or dipolar effects) without gravity.
 * This function does not update the molecule's potential.
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param laser Vector of lasers.
 * @param t Current time.
 * @param params Simulation parameters.
 * @return Potential energy of the molecule.
 */
double new_pot(vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
               const Field& fieldE, const vector<Laser>& laser, const double t,
               FitParams& params);

/**
 * Updates the potential energy of a molecule's internal state.
 * Dipolar potential is excluded to avoid detuning accumulation.
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param laser Vector of lasers.
 * @param t Current time.
 * @param params Simulation parameters.
 */
void set_pot_mol(vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                 const Field& fieldE, const vector<Laser>& laser, const double t,
                 FitParams& params);

/**
 * Updates the potential energy of all molecules in the system.
 * Dipolar potential is excluded to avoid detuning accumulation.
 * @param Mol Vector of molecules.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param laser Vector of lasers.
 * @param t Current time.
 * @param nb_mol Number of molecules in the system.
 * @param params Simulation parameters.
 */
void set_pot_all_mol(vector<Molecule>& Mol, const Field& fieldB, const Field& fieldE,
                     const vector<Laser>& laser, const double t, const int nb_mol,
                     FitParams& params);

/**
 * Updates and returns the acceleration of a single molecule.
 * The acceleration is calculated as: acc = -grad(E_pot) / mass.
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param laser Vector of lasers.
 * @param t Current time.
 * @return Acceleration vector of the molecule.
 */
Vecteur3D new_acc(vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                  const Field& fieldE, const vector<Laser>& laser, const double t);

/**
 * Updates the acceleration of all molecules in the system.
 * @param Mol Vector of molecules.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param laser Vector of lasers.
 * @param t Current time.
 * @param nb_mol Number of molecules in the system.
 * @param params Simulation parameters.
 */
void set_acc_all_mol(vector<Molecule>& Mol, const Field& fieldB, const Field& fieldE,
                     const vector<Laser>& laser, const double t, const int nb_mol,
                     FitParams& params);

/**
 * Calculates the energy shift of a molecule due to external fields (Zeeman, Stark, or dipolar).
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param laser Vector of lasers.
 * @param t Current time.
 * @param params Simulation parameters.
 * @return Energy shift of the molecule.
 */
double delta_shift(vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                   const Field& fieldE, const vector<Laser>& laser, const double t,
                   FitParams& params);

/**
 * Calculates the field-induced energy shift in s^-1 due to given B and E fields.
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param Bfield Magnetic field strength.
 * @param Efield Electric field strength.
 * @return Field-induced energy shift in s^-1.
 */
double delta_field_shift_B_E(vector<Molecule>& Mol, const int n_mol, const double Bfield,
                             const double Efield);

/**
 * Calculates the dipolar potential of a molecule's state.
 * This includes summation over transitions and interactions with other states.
 * @param Mol Vector of molecules.
 * @param n_mol Index of the molecule.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @param laser Vector of lasers.
 * @param t Current time.
 * @param params Simulation parameters.
 * @return Dipolar potential of the molecule's state.
 */
double delta_dipolaire(vector<Molecule>& Mol, const int n_mol, const Field& fieldB,
                       const Field& fieldE, const vector<Laser>& laser, const double t,
                       FitParams& params);

/**
 * Diagonalizes energy levels and transition frequencies in the presence of B and E fields.
 * Useful for complex interactions and precise energy calculations.
 * @param Internal_state_in Initial internal state of the molecule.
 * @param Internal_state_out Final internal state of the molecule.
 * @param pol_locale Local polarization vector.
 * @return Result of the diagonalization process.
 */
double diagonalization_new_Levels_Lines(Internal_state const& Internal_state_in,
                                        Internal_state const& Internal_state_out,
                                        Vecteur3D const pol_locale);

#endif
