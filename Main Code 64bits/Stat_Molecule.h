/*
  Name:
  Copyright:
  Author:
  Date: 23/10/06 10:01
  Description: Paramètres (température, positions, dispersion, moyenne des molécules) (utilisé dans Affichage_Mol.h)


*/

#ifndef Stat_Mol_SEEN
#define Stat_Mol_SEEN

#include <iostream>
using namespace std;


#include "Molecule.h"
#include "one_body.h"
#include <gsl/gsl_statistics.h> // pour les stat (mean, variance, ...)
#include "algorithmes.h"    // Pour qsort
#include "params.h"

/************************************************************************/
/****************** structure des données Stat *************************/
/************************************************************************/

// Paramètres statistiques d'une assemblée de Molecule;
// Moyenne, ecart type (= standard deviation) des vitesses et position
// Energie et température
// Il y a plusieurs température:
// T_cin sigma_(v-<v>)^2 /k_Boltzmann
// Mais afin d'éviter les molécules qui ont "explosées" (ejectées du piège et avec grande vitesses), on prend plutôt m sigma^2 /k_Boltzmann où sigma = standard deviation of the velocity sqrt(v^2)
//  Enfin Temp_ordonnee  est donné par ordonancement de la norme des vitesse des molécules pour extraire une température effective en utilisant 1/2 m v_i^2


class Stat
{
public:
    int population; // population du niveau étudié

    Vecteur3D mean_v;
    Vecteur3D mean_pos;
// ecart type (racine carrée de la variance)
    Vecteur3D sigma_v;
    Vecteur3D sigma_pos;

    double E_pot;  // sum_i [q V_elec(r_i) - m g z_i + hbar [E_cm-Ecm0]   + 0.5 sum_j not i qi qj /(4 pi epsilon_0 r_ij) ]
    double E_cin; // 1/2 m (v_x^2 +v_y^2 + v_z^2)

    Vecteur3D Temp_cin;  // from 1/2 m <v-<v>>^2
    Vecteur3D Temp_1D_50; // for 50% of the best particles
    double Temp_3D_50; //  / T_3D^2 = (Tx_50^2+Ty_50^2+Tz_50^2)/3

//constructeur
    Stat(): population(0), mean_v(Vecteur3D(0.,0.,0.)), mean_pos(Vecteur3D(0.,0.,0.)), sigma_v(Vecteur3D(0.,0.,0.)), sigma_pos(Vecteur3D(0.,0.,0.)), E_pot(0.), E_cin(0.), Temp_cin(Vecteur3D(0.,0.,0.)),  Temp_1D_50(Vecteur3D(0.,0.,0.)), Temp_3D_50(0.) {};
};



// Paramètres statistiques d'une assemblée de Molecule;
// Le dernier niveau donne la statistique globale (tout niveaux confondus)
void stat_molecule(const vector <Molecule> &Mol,  vector <Stat> &stat_Mol, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const int Nb_Mol, const int Nb_state, FitParams &params);

// Paramètres statistiques d'une assemblée de Molecule; Retourne le nb de molécules n dans l'état étudié
int stat_molecule_un_niveau(const vector <Molecule> &Mol, Stat &stat, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const int num_niveau_etudie, const int Nb_Mol, FitParams &params);
// Statistique sur les trajectoires
// distance Min, moyenne (somme) et max

int stat_molecule_form_list(const vector <Molecule> &Mol,  vector <Molecule> &liste_Adresses_Mol_dans_niveau, Stat &stat, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const int Nb_Mol, FitParams &params);

void stat_traj_molecule(const vector <Molecule> &Mol, double & dist_min, double & dist_sum, double & dist_max, int & nb_call, const int Nb_Mol=1);


#endif
