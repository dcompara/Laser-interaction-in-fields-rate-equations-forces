// CALCUL DU SHIFT EN ENERGIE DE LA MOLECULE
// Zeeman pour  l'instant mais Stark, dipolaire ensuite, ils peuvent dépendre du temps


#include "Molecule.h"
#include "Laser.h"
#include "Field.h"
#include "Initialisation_programme.h"     // Pour les champs
#include "params.h"

#include <iostream>
using namespace std;

#ifndef Shift_molecule_SEEN
#define Shift_molecule_SEEN

// Calcul (sans mise à jour des énergies)  du  potentiel de l'état interne de la molécule (Zeeman, Stark, dipolaire) sans gravité
double  new_pot(vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params);

// Mise à jour du  potentiel de l'état interne de la molécule
// Ne prend pas en compte le potentiel dipolaire pour éviter des accumulation de detuning
void  set_pot_mol(vector <Molecule> &Mol, const int n_mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params);

// Mise à jour du potentiel de l'état interne de toutes les molécules
// Ne prend pas en compte le potentiel dipolaire pour éviter des accumulation de detuning
void  set_pot_all_mol(vector <Molecule> &Mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params);



// Mise à jour  (et retourne) de l'accélération de la molécule
// Acc = -grad E_pot/mass
Vecteur3D  new_acc(vector <Molecule> &Mol, const int n_mol,   const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t);

// Mise à jour de l'accélération de toutes les molécules
void  set_acc_all_mol(vector <Molecule> &Mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params);



// Shift de la particule lié aux champ: zeeman, Stark ou dipolaire
double  delta_shift(vector <Molecule> &Mol, const int n_mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t,  FitParams &params);


// Calcul du shift en champ E ou B in s^-1
double  delta_field_shift_B_E(vector <Molecule> &Mol, const int n_mol, const double Bfield, const double Efield);


// Potentiel dipolaire de l'état j par Σ i>j delta_ij gamma_ij/Gamma_ij - Σ k<j delta_jk gamma_jk/Gamma_jk
double delta_dipolaire(vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t,  FitParams &params);

// For more complex case we want to diagonalize the Energy levels and hte transition Frequncies in B, E fields
double diagonalization_new_Levels_Lines(Internal_state const & Internal_state_in, Internal_state const & Internal_state_out, Vecteur3D const pol_locale);


#endif
