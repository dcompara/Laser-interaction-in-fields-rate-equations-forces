/*
  Name: Algorithme de calcul des taux
  Author: Daniel Comparat
  Date: 16/12/08


REMARQUES:
La prédissociation est traitée comme un parametre phénoménologiue
Les transition a deux photons aussi
la desexcitation dans le continuuum peut être prise en compte car parfois Sum_FC != 1 la partie non 1 est le continuum


ATTENTION
Unités des fichier (cm^-1)
Les lasers sont CW.
 Dans le cas d'un laser femto il ne faut pas forcément mettre Imoyen cf thèse d'Antoine Monmayrant (IV.5)

*/


#ifndef Transition_rate_SEEN
#define Transition_rate_SEEN

#include <iostream>
using namespace std;

#include  <iostream>                       // to include cout, cin
#include  <cmath>                          // to include sqrt(), etc.
#include "laser.h"                       // Classe Laser
#include  "shift_molecule.h"              // Pour le shift en énergie delta
#include  "sortie_donnees.h"              // Pour le shift en énergie delta
#include  "Kinetic_Monte_Carlo.h"       // to include the algorithm
#include  "algorithmes.h"               // to include Is_switch laser
#include <complex>                      // Pour les calculs d'émission spontanée

#include <Eigen/Eigen>  // For matrix manipulation
using namespace Eigen;

#include "params.h"     // Pour mettre des paramètres dedans

/************************************************************************/
/************************** Interaction laser - 2 niveaux *************/
/************************************************************************/

// Taux d'emission spontanée d'une transition (sans prendre en compte l'ordre des états:
double Gamma_spon(const double dip_Debye, const double energie_cm);


// Forme gaussienne Gamma = FWHM  = 2 \sqrt{2 \ln(2)} \sigma  = 2.3548 \sigma
inline double Gauss(const double I_tot, const double delta, const double Gamma);

// Forme laser plate
// retourne l'intensité (spectrale) au décalage delta souhaité
inline double intensity_flat(const double I_tot, const double delta, const double Gamma);

// Forme  lorentzienne
// I*Gamma/(2.*pi)/(Gamma^2/4 + delta^2)
inline double Lorentz(const double I_tot, const double delta, const double Gamma);

// Approximation of the Voigt profile [Lorentz*Gaussian](omega) where Lorentz and Gaussian FWHM values of G_L and G_G and centered on 0.
double Pseudo_Voigt(double delta, const double G_L, const double G_G);


// comb spectum = gaussian laser but with individual Lorentzian comb lines
// the comb line are positioned at nu_offset + n*nu_repetition
// The intensity is taken at the detuning of the nearest (lorentizian) comb line
// The linewidth is due to the spontaneous emission of the 2 levels system drived by the comb laser
double comb_shape(double I,double delta, const double G_L, const double G_G, const int n_las, const double Energy_transition_cm, FitParams &params);


// intensity of a pseudo_BBR (pseudo_Black Body Radiation) that is a BBR spectrum but with a gaussian beam
inline double pseudo_BBR_intensity(const double I, const double Energy_transition_cm, const double T);



// Calcul de l'intensité locale [Lorentz*I](omega) with a natural Cauchy-Lorents linewitdh with FWHM Gamma_atomic.
// delta = omega_Laser - omega_atomic
double intensity_Convolution_linewidth(const double I, const double delta, const double Gamma_atomic, const double Gamma_Laser, const int laser_type, const int n_las, const double Energy_transition_cm, const double Energy_transition_laser_cm, FitParams &params);

// Taux à deux niveau Gamma'
// is_bound_transition is here to detcte transition (ionization) to continuum (0 = false that is ioniozation ; 1 =  true that is bound - bound transition
double rate_two_level(const double I_loc, const double dip_Debye);

// Rate for (photo-)ionization
double rate_ionization(const double I_loc, const double dip_Debye, const double delta, const double Energy_transition_cm);

// Field ionization rate
double rate_field_ionization(const Laser& my_laser,  const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE);

// Calculation of the excitation rate ofr a laser (or a class of coherent laser) for a given transition (bound or ionization)
// Add this rate only is is_rate_calculated is true
double rate_excitation(vector <type_codage_react> &reaction_list, vector <double> &rate, const int n_las, const vector <Molecule> &Mol, const int n_mol, const Vecteur3D k, const double Itot_loc,  const double Itot,
                       const double dipole_debye, const double delta,  const double Energy_transition_cm,  const Field &fieldB, const Field &fieldE, const int is_bound_transition, const Internal_state Internal_state_out, const bool is_rate_calculated ); // calculate rate and reaction_list for this transition

/************************************************************************/
/************************** Fonctions transition Molecule ***************/
/************************************************************************/

// Calcul de tous les taux d'emission spontanée de la molécule
int rates_molecule_spon(vector <Internal_state> &Level, vector <type_codage_react> &reaction_list, vector <double> &rate, const Molecule &my_mol, const Field &fieldB, const Field &fieldE, const int num_mol, FitParams &params);


// Calcul of rate (if no interferece between laser)  between level in and out for a given laser and for a given molecule. This add the light shift effect to the dipolar potential (delta_pot_dipolaire)
// update the detuning (delta), the polarization vector also updated), Gamma_spot_total, sqrt of the local intensity of the laser
int rates_single_molecule_laser_level(const int n_las, double dipole, double &delta, double &eps_pol, double &Gamma_spon_tot, double sqrt_intensity_loc[],
                                      vector <Internal_state> &Level, Internal_state &Internal_state_in, Internal_state &Internal_state_out, vector <type_codage_react> &reaction_list, vector <double> &rate,
                                      const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const Laser &my_laser, const double t, double &delta_pot_dipolaire,
                                      FitParams &params, bool is_rate_calculated, int is_bound_transition=1, const int n_level_in =0, const int n_level_out =1, const double Gamma_in=0., MatrixXd d[]= {});


// Calcul de tous les taux de la molécule  pour tous les lasers
// On les ajoutes aux taux des autres
//  return nb_rate
// It calculate also the dipole detuning which is just rate * delta/Gamma
// In this last case we do not necessarily to calculate the rate, as for the dipolar shift, so is_rate_calculated would be false
int rates_molecule(vector <Internal_state> &Level, vector <type_codage_react> &reaction_list, vector <double> &rate,
                   const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, double &delta_pot_dipolaire, FitParams &params, bool is_rate_calculed = true);






// return the decay rate from Level i = sum_j<i Gamma_ij
double Gamma_Level_from_diagonalized_dipole(vector <Internal_state> &Level, MatrixXd d[], const int i);



// Copie tous les taux sauf celui de la molécule numero_mol
// Retourne le nombre de taux de ces molécules
int copie_rates_molecules(vector <type_codage_react> &reaction_list, vector <double> &rate, const int numero_mol = -1, const int nombre_old_rate = 0);


// Calcul de tous les taux de toutes les molécules si numero_mol = aucune (-1)
// Sinon on ne recalcule que celui de la molécule numero_mol
int calcul_rates_molecules(vector <Internal_state> &Level, MC_algorithmes Algorithme_MC, vector <type_codage_react> &reaction_list, vector <double> &rate, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
                           const double t, const int numero_mol, const int N_Mol, FitParams &params );



/***********************************************/
/********** réalisation d'une réaction *********/
/***********************************************/


// Scan all rates and realize the reaction if the rate*random < dt=0.1/rate_max
// We scan all molecules and make at most one rate per molecule but
// Typically  a tens of Molecules will be affected
int do_reaction_FastRoughMethod(const MC_algorithmes Algorithme_MC, const gsl_rng * r, const vector <double> &rate, const vector <type_codage_react> &reaction_list, vector <Molecule> &Mol,
                                const int n_reac, const vector <Laser> &laser, const double dt_KMC, ofstream & file_rate, int &number_photons, FitParams &params);

// Effectue la réaction. Ne met pas les potentiels à jour ensuite. Il faut le faire à part.
// Retourne le numéro de la molécule affectée
int do_reaction(const MC_algorithmes Algorithme_MC, const gsl_rng * r,  const vector <double> &rate, const vector <type_codage_react> &reaction_list,
                vector <Molecule> &Mol, const int n_reac, const  vector <Laser> &my_laser, const double dt_KMC, ofstream & file_rate, bool first_call, int &number_photons, FitParams &params);

// Gives a random unit vector (in the lab frame) for the spontaneous emission for a transition delta_M=-1,0,+1 and a quantization axis
// Or (if diagonalized) for a polarization vector e_pol:
//  e_pol[q] = normalized dipole transition <i|d_(q)|j>; q=-1,0,1. So in the quantification axis e_pol = sum_q epol_q E^q
// the probability distribution linked with f(r)= (3/8π)[1-|r.e_pol|^2]
Vecteur3D get_unit_vector_spontaneous_emission(const gsl_rng * r, Vecteur3D e_pol_dipole_transition, Vecteur3D quantization_axis, int delta_M, FitParams &params);


#endif
