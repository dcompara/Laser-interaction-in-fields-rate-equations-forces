// Initialisation des niveaux, des transitions, des forces de raies et des durées de vie
// Des paramètres lasers, des positions, vitesses des molécules ....

#ifndef Initialisation_programme_SEEN
#define Initialisation_programme_SEEN


#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>                   // for random generator
#include <gsl/gsl_randist.h>               // for gaussian random generator
#include <gsl/gsl_cdf.h>

#include "Vecteur3D.h"
#include "Molecule.h"                       // Classe molecule
#include "Laser.h"                       // Classe Laser
#include "Field.h"                       // Classe Champ
#include "Transitions_initialisation.h"
#include "shift_molecule.h"
#include <sstream>  // pour ostringstream




/************************************************************************/
/******************* INITIALISATION DU PROGRAMME ************************/
/************************************************************************/


// Initialisation des niveaux, des transitions, des forces de raies et des durées de vie
// Retroune le nb de niveaux
int   initialisation_trans_mol(const char *nom_file_Levels, const char *nom_file_Lines ,   vector <Internal_state> &Level, FitParams &params);

//Initialise l'état du générateur de nombre aléatoire
void set_random_generator(gsl_rng  * r, int nombre_seed, const char *nom_file);


#include "params.h"
// Initialise la position et la vitesse (gaussienne) des molecules
void Init_Molecule(const gsl_rng * r, vector <Molecule> &Mol, const Field fieldB, const Field fieldE, const int Nb_of_type_of_different_Mol, FitParams &params, DataCards data);

// Initialise les  lasers
void Init_Laser(vector <Laser> &laser, const int Nb_lasers, FitParams &params, const char *nom_file_spectrum, const char *nom_file_intensity);

// Modification temporelle des paramètres
void Modif_Param(const double t ,  FitParams &params);


/************************************************************************/
/******************* INITIALISATION DES CHAMPS ************************/
/************************************************************************/

#include "params.h"
// Initialise les champs, zero par défaut
void Init_Field(Field &fieldB, Field &fieldE, FitParams &params, const char *nom_Magn_file=NULL, const char *nom_Elec_file=NULL);

#endif
