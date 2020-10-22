/*

Program to read Line or Energy list from


LES UNITES ICI SONT CELLES DES FICHIERS (cm^-1 typiquement)

Avant j'utilisaait des fichiers Franck-Condon et Hönl-London
Je les laisse ici mais il ne servent plus
La notation  [NA][NX] indique qu'il y a
NA lignes et NX colonnes dans les tableaux
Il faut donc bien vérifier que dans les fichiers
vA indice les lignes et vX indice les colonnes

En effet les boucles sont vA puis vX


*/


#ifndef Transition_init_SEEN
#define Transition_init_SEEN

#include <iostream>
using namespace std;

#include "Internal_state.h"
#include <gsl/gsl_rng.h>                   // for random generator
#include <gsl/gsl_randist.h>               // for gaussian random generator

#include "Molecule.h"                       // Classe molecule
#include "params.h"
#include "datacards.h"                      // Pour lire le fichier de paramètres

#include <gsl/gsl_sf_coupling.h> // POUR les 3j de Wigner
// double gsl_sf_coupling_3j (int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc)
// gsl_sf_coupling_3j (2*ja,2*jb,2*jc,2*ma,2*mb,2*mc);
// compute the Wigner 3-j coefficient,

// ja jb jc
// ma mb mc

// where the arguments are given in half-integer units, ja = two ja/2, ma = two ma/2
// Gives 0 in wrong cases


/*** Program to read Line or Energy list from  ***/
// Manifold(1 for upper or 0 for lower typically)  M(in Field it is real M)  bound_level  #(number to discriminate the levels) Population J N  Energy0(in 0 field)   Delta   C
// Where for fit in Energy in Fields  Energy +/- Sqrt(Delta^2/4 + C^2*F^2) where F is field
// C= Linear coefficient is in cm-1/T
// Energy and Delta are in cm-1
// Retourne le nb de niveaux (1 si contient seulement Level[0])
 int Level_List(const char *nom_file, vector <Internal_state> &Level, FitParams &params);


/*** LINE LIST *****/
// Les fichiers  contiennent
 // UpperManifold	M'	bound_level'	#'	LowerManifold	M"	bound_level"	#"  Eupper	Elower	dip_Debye
// Read the file containing the position, assignment and intensity of lines in the simulated spectrum.
// Retourne le nb de transitions
int Line_List(const char *nom_file, vector <Internal_state> &Level, FitParams &params);

// Initialisation des durées de vies
void initialisation_Gamma( vector <Internal_state> &Level, FitParams &params);


// initalisation des probabilités de peuplement de départ
void initialisation_proba(const gsl_rng *r, vector <Molecule> &Mol, const int Nb_molecules, const vector <Internal_state> &Level);

// Initi
#endif
