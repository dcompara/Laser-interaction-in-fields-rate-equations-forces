// CALCUL DU SHIFT EN ENERGIE DE LA MOLECULE
// Zeeman- Stark, dipolaire ensuite

#include "shift_molecule.h"
#include "one_body.h"
#include "algorithmes.h"
#include "Transition_rate_calcul.h"


#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>



// Calcul (sans mise à jour des énergies) du  potentiel de l'état interne de la molécule juste avec les champs
// Il n'y a pas la gravité elle est rajouté directement à la main dans l'accélération
double  new_pot(vector <Molecule> &Mol, const int n_mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params)
{
    // Vecteur3D pos=my_mol.get_pos();
    // double delta =delta_field_shift_B_E(my_mol, fieldB.get_Field(pos).mag(), fieldE.get_Field(pos).mag()); // Mettre si on ne veux pas le potentiel dipolaire
    double delta = delta_shift(Mol, n_mol,  fieldB, fieldE, laser, t,  params);
    // my_mol.Energy_cm = my_mol.Energy0_cm + delta/Conv_Ecm_delta;
    return HBAR*delta;
}


// Mise à jour du  potentiel de l'état interne de la molécule
// Ne prend pas en compte le potentiel dipolaire pour éviter des accumulation de detuning
void  set_pot_mol(vector <Molecule> &Mol, const int n_mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t,  FitParams &params)
{
    Molecule my_mol = Mol[n_mol];
    double pot;
    Vecteur3D pos=my_mol.get_pos();
    double B=fieldB.get_Field(pos).mag(); // Magnétique
    double E=fieldE.get_Field(pos).mag(); // Electrique
    pot =  HBAR * delta_field_shift_B_E(Mol, n_mol, B, E); // Recalcul le shift (dipolaire, magn, électrique) de cette molécule
    my_mol.Energy_cm = my_mol.Energy0_cm + (pot/HBAR)/Conv_Ecm_delta;
    return;
}

// Mise à jour du potentiel de l'état interne de toutes les molécules
// Ne prend pas en compte le potentiel dipolaire pour éviter des accumulation de detuning
void  set_pot_all_mol(vector <Molecule> &Mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params)
{
    for (int nb = 0; nb < nb_mol; nb ++)
        set_pot_mol(Mol, nb, fieldB, fieldE, laser, t,  params);

    return;
}



// Mise à jour (et retourne)  de l'accélération de la molécule
Vecteur3D  new_acc(vector <Molecule> &Mol, const int n_mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t)
{
    Molecule my_mol = Mol[n_mol];

    Vecteur3D grad_pot_B,grad_pot_E;
    Vecteur3D pos=my_mol.get_pos();
    double B=fieldB.get_Field(pos).mag(); // Magnétique
    double E=fieldE.get_Field(pos).mag(); // Electrique
    Vecteur3D grad_field_B2=fieldB.get_grad_field_F2(pos); // gradient de ||B||^2 (Magnétique)
    Vecteur3D grad_field_E2=fieldE.get_grad_field_F2(pos); // Electrique

    grad_pot_E = my_mol.Grad_Energy_Shift_E(E, grad_field_E2);
    grad_pot_B = my_mol.Grad_Energy_Shift_B(B, grad_field_B2);

    Vecteur3D acc;

    acc = gravity - grad_pot_E/my_mol.get_mass()  - grad_pot_B/my_mol.get_mass() ;

    my_mol.set_acc(acc);

    return acc;
}


// Mise à jour de l'accélération de toutes les molécules
void  set_acc_all_mol(vector <Molecule> &Mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params)
{
    for (int nb = 0; nb < nb_mol; nb ++)
        new_acc(Mol, nb, fieldB, fieldE, laser,  t );

    return;
}



// Shift of particl i
// Due tio fields: zeeman, Stark, dipolar
double  delta_shift(vector <Molecule> &Mol, const int n_mol,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params)
{
    Molecule my_mol = Mol[n_mol];

    // double delta= delta_dipolaire_shift(my_mol, my_laser, params); // Dipolaire
    double delta=0.; // Pas de dipolaire

    Vecteur3D pos=my_mol.get_pos();
    double B=fieldB.get_Field(pos).mag(); // Magnétique
    double E=fieldE.get_Field(pos).mag(); // Electrique
    delta +=delta_field_shift_B_E(Mol, n_mol, B, E);
    delta += delta_dipolaire(Mol, n_mol, fieldB, fieldE, laser, t, params);
    return delta;
}


// Calcul du shift en champ E ou B in s^-1
double  delta_field_shift_B_E(vector <Molecule> &Mol, const int n_mol, const double Bfield, const double Efield)
{
//    cout << "B " << Bfield << endl;
//    cout << " my_mol.Energy_Shift_B_cm " << my_mol.Energy_Shift_B_cm(Bfield) << endl;
    Molecule my_mol = Mol[n_mol];
    double delta =  (my_mol.Energy_Shift_B_cm(Bfield) + my_mol.Energy_Shift_E_cm(Efield))*Conv_Ecm_delta ;//  cm^(-1) -> s^-1
    return  delta;
}


/*** Calcul du shift dipolaire
// Gamma (force de raie Coeff Einstein de Cette raie précisement) et GammaTot = Gamma + tout autres élergissement (laser, collision, desexcitations vers d'autres niveaux ...)
// on prend par définition hbar delta_dip = delta Gamma'/Gamma_Tot si l'état est bas en énergie et - cela si l'inverse
// où Gamma'=taux d'absorpion = taux d'emission stimulée = Gamma_Tot OmegaRabi^2/(Gamma_tot^2 + 4 delta^2)
// De plus pour les lasers façonnés cela peut poser problème.

// Evidemment si il y a 2 lasers à 3 niveaux déjà cela pose problème en général si s>>1. Il faut habiller proprement les niveaux.
En s'inspirant de l'article "an atom faucet" nous cherchons une formule ne dépassant pas l'effet à haute intensité et correcte à basse.
Ainsi le décalage delta_j ln(1+s_j) crée par le laser j devient
[(Sum_j delta_j s_j) ln (1+sum_j s_j)]/sum_j s_j qui est correct à basse intensité (s_j <<1) et aussi si on fait I= I/2+I/2
EN fait nous n'utilisons pas cette formule encore
***/

/**
 Même si il n'y a pas de lumière à la raie correspondante on calcul, comme si il y en avait, le détuning
Cela semble correct pour des spectres avec peu de trous, mais si on a des peignes il faut prendre I_reel et non I_tot avant façonnage
**/


// Potentiel dipolaire(/Hbar) de l'état j par Σ i>j delta_ij gamma_ij/Gamma_ij - Σ k<j delta_jk gamma_jk/Gamma_jk
// Is calculated in the sur la fonction rates_molecule
double delta_dipolaire(vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params)
{
    double delta_pot_dipolaire =0.;
    vector <double> rate_local;
    vector <type_codage_react> reaction_list_local;
    vector <Internal_state> Level; // Require to call the function

// TODO (dc#5#): To be modified if we want to have the correct dipolar force in diagonalized case. Here just OK for non diagonalized case so I put fake d0 and d matrix

    MatrixXcd H,E0_cm, Zeeman_cm_B, d[3];
    MatrixXd d0[3];

    rates_molecule(Level, reaction_list_local, rate_local, Mol, n_mol, fieldB, fieldE, laser, t,  delta_pot_dipolaire,  params, H,E0_cm, Zeeman_cm_B, d0, d, false);

    return delta_pot_dipolaire;
}






