/*

Program to read Line or Energy list from Pgopher


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


/*** Program to read Line or Energy list from Pgopher ***/
// Manifold(1 for upper or 0 for lower typically)  M(in Field it is real M)  Sym(Parity +/-)  #(number to discriminate the levels) Population J N  Energy0(in 0 field)   Delta   C
// Where for fit in Energy in Fields  Energy +/- Sqrt(Delta^2/4 + C^2*F^2) where F is field
// C= Linear coefficient is in cm-1/T
// Energy and Delta are in cm-1
// Retourne le nb de niveaux (1 si contient seulement Level[0])
 int Pgopher_Level_List(const char *nom_file, vector <Internal_state> &Level, FitParams &params);


/*** LINE LIST *****/
// Les fichiers  contiennent
 // UpperManifold	M'	Sym'	#'	LowerManifold	M"	Sym"	#"	Position	Intensity	Eupper	Elower	Spol
// Read the file containing the position, assignment and intensity of lines in the simulated spectrum.
// Retourne le nb de transitions
int Pgopher_Line_List(const char *nom_file, vector <Internal_state> &Level, FitParams &params);


// initalisation des probabilités de peuplement de départ
void initialisation_proba(const gsl_rng *r, vector <Molecule> &Mol, const int Nb_molecules, const vector <Internal_state> &Level);



/************************************************************************/
/************************** INITIALISATION FC, Energie, HL **************/
/************************************************************************/

// Initialisation des durées de vies
void initialisation_Gamma( vector <Internal_state> &Level, FitParams &params);


// Modification des fichiers de niveaux et de transitions en fonction des FC et des Energies
void Modification_FC_trans_mol(const char *nom_file_Levels, const char *nom_file_Lines ,  vector <Internal_state> &Level, FitParams &params, DataCards data);


// Initialization Of FC factors for absorption and emission
void initialisation_FC(double **FCondon, const char *nom_file_FC, const int NvAmax, const int NvXmax);


// Energy of vibrational levels
// The file format is v Ev Bv (en cm-1)
void initialisation_energie(double **Ev, const char *nom_file_E_v, const int Nvmax, const int NJmax, const double E_v0_J0 = 0.);

// Crée le nouveau fichier de données à partir des FC.
// On lit un fichier de base contenant juste vX=0 et vA=0. Ensuite des FC et des Energie Ev et Bv des niveaux et on forme le nouveau ficheir
void New_Pgopher_Level_Line(const char *nom_file_initial_Level, const char *nom_file_initial_Line, vector <Internal_state> &Level, double **FC, double **EvX, double **EvA, int NXout, int NAout, int N_Two_JA_out_max, int N_Two_JX_out_max, FitParams &params);



// Initi
#endif
