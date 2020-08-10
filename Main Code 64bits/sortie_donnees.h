

#ifndef Sortie_Donnee_SEEN
#define Sortie_Donnee_SEEN

#include <iostream>
using namespace std;


#include <gsl/gsl_rng.h>                   // for random generator
#include <gsl/gsl_randist.h>               // for gaussian random generator

#include  <stdlib.h>
#include   <stdio.h>
#include  <fstream>                        // to read the data and put them in files

#include <iostream>                // for setprecision
#include <iomanip>

#include "Molecule.h"                       // Classe molecule
#include "Laser.h"                       // Classe Laser
#include "Field.h"                       // Classe Laser

#include "Stat_Molecule.h"                       // statistiques sur un ensemble de molecules
#include "Transition_rate_calcul.h"                     // pour le spectra façonné

#include "params.h"                         // Pour les paramètres
#include "datacards.h"                         // Pour les paramètres
#include "diagonalization.h" // Pour les paramètres


//Sauve l'état du générateur de nombre aléatoire pour ne jamais reprendre le même
void save_random_generator(gsl_rng * r,  const char *nom_file);


void Sortie_donnee_example(ofstream & file_out,  vector <Molecule> &Mol, vector <Internal_state> &Levels, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol, FitParams &params,  DataCards &data, const int number_photons);


void Sortie_donnee(ofstream & file_out,  vector <Molecule> &Mol, vector <Internal_state> &Levels, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol, FitParams &params,  DataCards &data, const int number_photons);


// Toutes à la suites en temps
void Sortie_donnee_pop_vJ(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, const int N_two_JX,  FitParams &params);

// Sortie des populations dans l'état vX à la suite les unes des autres en temps
void Sortie_donnee_pop_v(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, FitParams &params, int number_photons);

void Sortie_rate_example(ofstream & file_rate,  const vector <double> &rate, vector <Internal_state> &Level, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int N_Mol, const double t, FitParams &params);

void Sortie_rate(ofstream & file_rate,  const vector <double> &rate, vector <Internal_state> &Level, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int N_Mol, const double t, FitParams &params);

// Sortie du spectre laser
void Sortie_laser_spectrum(ofstream & file_out, const vector <Laser> &laser, FitParams &params, int num_laser=0);


#endif
