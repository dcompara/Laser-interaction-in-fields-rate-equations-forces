/*
  Name: Algorithme de calcul Monte Carlo cinétique
  Author: Daniel Comparat
  Date: 17/09/06 08:14

Résout les équations de taux

Il y a plusieurs algorithmes
cf
Kinetic Monte Carlo modeling of dipole blockade in Rydberg excitation experiment
Amodsen Chotia, Matthieu Viteau, Thibault Vogt, Daniel Comparat and Pierre Pillet
New Journal of Physics 10 (2008) 045031

+ Fast Rough Method = several Random choice to evolve typically N_Mol/10 times faster.

*/





#ifndef KMC_SEEN
#define KMC_SEEN

using namespace std;

#include <vector>
#include  <iostream>                       // to include cout, cin
#include  "algorithmes.h"                  // to include binarySearch
#include  <cmath>                          // to include sqrt(), etc.
#include <gsl/gsl_rng.h>                   // for random generator
#include <gsl/gsl_randist.h>               // for gaussian random generator
#include  "constantes_SI.h"


enum MC_algorithmes { Aucun_MC = -1,
                      Kinetic_Monte_Carlo = 0,
                      Random_Selection_Method = 1,
                      First_Reaction_Method = 2,
                      Fast_Rough_Method = 3
                    };



// Permet de tirer avec une proba uniforme dans une liste rate[]
// Il faut utiliser gsl_ran_discrete_preproc si on doit tirer plusieur fois ?

int tirage_random_discrete(const gsl_rng *r, const vector <double> &array, const int nb_rate, double &sum_rate_tot);




/************************************************************************
************************************************************************
 FONCTION PRINCIPALE POUR LE  Monte Carlo


Les fonctions (selon l'algorithme choisi) retournent le numéro de la réaction et changent le temps.
*/



/*****
THE KMC ALGORITHM IS THE FOLLOWING
• Initializing the system to its given state called k at the actual time t.
• Creating a new rate list Gamma_lk for the system, l = 1, . . . , N.
• Choosing a unit-interval uniform random number generator r : 0 < r < 1 and calculating
the first reaction rate time t0 by solving Int_t^t_0 Sum_l Gamma_kl(t) = - ln r
(for constant rate it is Delta_t * Sum_l Gamma_kl = - ln r
• Choosing the new state l from the unifor probability of the rate Gamm_kl(t_0)

*/


// Voici l'algorithme pour des taux INDEPENDANT du temps (c'est généralement notre cas).
// Un algorithme plus poussé avec évolution N_corps existe dans le programme Rydberg
// retourne le numéro de la réaction
int KMC_step(const gsl_rng * r, const  vector <double> &rate,  double &dt);


/************************************************************************
************************************************************************
 FONCTION PRINCIPALE POUR LE Kinetic Monte Carlo avec First Reaction Method :

ALGORITHME AVEC TAUX DEPENDANT DU TEMPS
An Introduction To Monte Carlo Simulations Of Surface Reactions
A. P. J. Jansen: cond-mat/0303028


Evolution temporelle First Reaction Method (FRM)
appellée aussi Discrete Event Simulation


According to this method, when the system is in a given configuration a,
the set of all possible reactions is determined,
and a time of occurrence tab is generated for EACH reaction
Then, the reaction with the SMALLEST tab  is selected,
the configuration is changed accordingly, and the time t is incremented in tab
Finally, the set of possible reactions is generated according
to the new configuration b, and the procedure is repeated

CE QUI CHANGE comparé au KMC précédent est que l'on calcule
integralle_{t_0}^{t_0+\Delta t_i} rate_i[t] dt =  -\log u_i.
Et on choisi le plus petit \Delta t_i, i.e. la première réaction.


************************************************************************/


// Cherche le temps minimal pour l'évolution FRM
// Retourne le numéro de la réaction ayant ce temps minimal (la FRM !)

int FRM_step(const gsl_rng * r, const  vector <double> &rate,  double &dt);

/*
Evolution temporelle Random Selection Method (RSM)
Choix aléatoire de l'atome et de sa réaction et du temps.
Parfois il n'y a pas de réaction
*/

int RSM_step(const gsl_rng * r, const  vector <double> &rate, double &dt);


// Calculate the max rate and then give the time for an evolution 1/rate_Max
int Fast_Rough_Method_step(const gsl_rng * r, const  vector <double> &rate, double &dt);

// Trouve la réaction (son numéro dans la liste des taux rate[]
// Trouve le temps où elle doit être effectuée
int find_reaction(MC_algorithmes Algorithme_MC , const gsl_rng * r, const vector <double> &rate, double & dt);


#endif

