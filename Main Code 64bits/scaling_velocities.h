/*
Scaling velocities: to make mean temperature of the clouds converge to the initial value that is entered in liste_param.h
*/


//#ifndef scaling_velocities_SEEN
//#define scaling_velocities_SEEN





#include "Atome.h"
#include "Internal_state.h"
#include "Vecteur3D.h"
#include "Laser.h"                       // Classe Laser
#include "Molecule.h"
#include "one_body.h"
#include "Stat_Molecule.h"


#include "algorithmes.h"    // Pour qsort
#include "params.h"


#include "Field.h"

// fonction rescale velocities - applique Lambda aux vitesses, à un même temps t.
void scaling_velocities( const double t, vector <Molecule> &Mol, const int nb_mol,const vector <Internal_state> &Level, const Field &fieldB,const Field &fieldE, FitParams &params,double coupling_efficiency, ofstream & file_scal);
// fonction qui calcule le coefficient lambda de Berendsen, à un tmeps t.
// Briefly v(t+dt) = v(t) * sqrt[1+coef(T0/T(t)-1)]
// So T(t+dt) ~  T(t)+coef(T0-T(t))
Vecteur3D calcul_lambda( const double t, vector <Molecule> &Mol, const int nb_mol,const vector <Internal_state> &Level, const Field &fieldB,const Field &fieldE, FitParams &params,double coupling_efficiency, ofstream & file_scal);
//fonction qui applique Le coefficient Lambda aux vitesses, à un temps t. (il faut l'appliquer à t+dt).
void rescaling_velocities_after_dt ( const double t, vector <Molecule> &Mol, const int nb_mol,const vector <Internal_state> &Level, const Field &fieldB,const Field &fieldE, FitParams &params,double coupling_efficiency, ofstream & file_scal, Vecteur3D Lambda);
