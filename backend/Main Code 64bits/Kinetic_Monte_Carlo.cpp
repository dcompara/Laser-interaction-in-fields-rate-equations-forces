/*
  Name: Algorithme de calcul Monte Carlo cinétique
  Author: Daniel Comparat
  Date: 17/09/06 08:14

Résoult les équations de taux

Il y a plusieurs algorithmes
cf
Kinetic Monte Carlo modeling of dipole blockade in
Rydberg excitation experiment
Amodsen Chotia1, Matthieu Viteau, Thibault Vogt, Daniel Comparat and Pierre Pillet
New Journal of Physics 10 (2008) 045031
ou le ReadMe

+  Fast Rough Method . That evolves several molecules during the same time step.


*/



#include  "Kinetic_Monte_Carlo.h"


// Permet de tirer avec une proba uniforme dans une liste rate[]
// Il faut utiliser gsl_ran_discrete_preproc si on doit tirer plusieur fois ?

int tirage_random_discrete(const gsl_rng *r, const vector <double> &array, const int nb_rate, double &sum_rate_tot)
{

    double *Sum_array = new double[(nb_rate+1)]; // cannot used vector because binarySearch
    Sum_array[0]=0.;
    for (int i = 0; i < nb_rate; i++)
        Sum_array[i+1] = Sum_array[i] + array[i];   // On crée la liste croissante de taux intégré en temps

    double u=1.-gsl_rng_uniform(r);     // pour tirer au hasard dans (0,1] car gsl_rng_uniform(r) dans [0,1)
    int sol = binarySearch(Sum_array, 0,  nb_rate, Sum_array[nb_rate]*u) ;
    // la réaction qui a été choisi est telle que R_{i-1} < u R <= R_i, i.e. rate_tot[nat] < u*R <= rate_tot[nat+1]

    sum_rate_tot = Sum_array[nb_rate];

    delete [] Sum_array;
    return(sol);
}






/************************************************************************
************************************************************************
 FONCTION PRINCIPALE POUR LE  Monte Carlo


 Remarque 1: Pour se retrouver dans le cas des théorèmes d'évolution indiquant que le KMC est correct
           il faut que TOUTES les réactions soient prises en compte par le KMC.

 Remarque 2: Les taux varient car les positions varient mais on néglige cela ici.

Les fonctions (selon l'algorithme choisi) retournent le numéro de la réaction et changent le temps.
*/



/*****
THE KMC ALGORITHM IS THE FOLLOWING
• Initializing the system to its given state called k at the actual time t.
• Creating a new rate list Gamma_lk for the system, l = 1, . . . , N.
• Choosing a unit-interval uniform random number generator r : 0 < r < 1 and calculating
the first reaction rate time t0 by solving Int_t^t_0 Sum_l Gamma_kl(t) = - ln r
(for constant rate it is Delta_t * Sum_l Gamma_kl = - ln r
• Choosing the new state l from the uniform probability of the rate Gamma_kl(t_0)

*/



// Voici l'algorithme pour des taux INDEPENDANT du temps (c'est généralement notre cas).
// Un algorithme plus poussé avec évolution N_corps existe dans le programme Rydberg
// retourne le numéro de la réaction
int KMC_step(const gsl_rng * r, const vector <double> &rate, double &dt)
{

    double sum_rate; // somme des taux
    int nb_rate = rate.size();
    int n_reac = tirage_random_discrete(r, rate, nb_rate, sum_rate); // Choix de la réaction (et calcul de la somme)

    // Choix du temps Delta_t = - ln random/ Sum_l Gamma_kl;
    double random = 1.-gsl_rng_uniform(r); //  uniform random number u \in (0,1]
    dt = - log(random)/sum_rate;  // Incrémentation aléatoire du temps pour rendre compte du caractère Poissonien du processus

    return n_reac;
}





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
integrall_{t_0}^{t_0+\Delta t_i} rate_i[t] dt =  -\log u_i.
Et on choisi le plus petit \Delta t_i, i.e. la première réaction.


************************************************************************/


// Cherche le temps minimal pour l'évolution FRM
// Retourne le numéro de la réaction ayant ce temps minimal (la FRM !)

int FRM_step(const gsl_rng * r, const  vector <double> &rate,  double &dt)
{
    int  nb_rate = rate.size();
    //Calcul, via  l'intégrale temporelle int_{t_0}^(t_0+t_KMC_estimé) des taux de transition, du temps t_KMC
    double *random_number = new double[nb_rate];
    for (int i=0; i<nb_rate; i++)
        random_number[i] = 1.-gsl_rng_uniform(r);   // Choix des nombres aléatoires pour l'algorithme Monte Carlo.

    int i_time_min = -1;
    dt =  VERY_LARGE_NUMBER; // On initialise le temps minimum à une valeur "infinie"

    double dtime;

    for (int i = 0; i < nb_rate; i++)
    {
        if (rate[i] > 0) dtime = (-log(random_number[i]))/rate[i];
        else dtime= VERY_LARGE_NUMBER;
        // On calcul aléatoirement un temps pour chaque transition

        if (dtime < dt)
        {
            i_time_min = i;   // On cherche le temps minimum
            dt = dtime;
        }
    }
    delete[] random_number;

    return  i_time_min;
}




/*
Evolution temporelle Random Selection Method (RSM)
Choix aléatoire de l'atome et de sa réaction et du temps.
Parfois il n'y a pas de réaction
*/

int RSM_step(const gsl_rng * r, const  vector <double> &rate, double &dt)
{
    int  nb_rate = rate.size();
    double rate_max = 0. ;   // Initialisation pour le taux de transition maximum
    for (int i = 0; i < nb_rate; i++) // Recherche des taux de transition maximum
    {
        if (rate[i] > rate_max)
            rate_max = rate[i];  // Taux de transition maximum
    }

    dt = 1/(nb_rate*rate_max) ; // Calcul d'un intervalle de temps "aléatoire" où une seule réaction au maximum est possible

    int n_reac = (int) (gsl_rng_uniform(r)*nb_rate); // Tirage au hasard d'un état entre 0 et nb_rate-1

    double Delta_t = gsl_rng_uniform(r)/rate[n_reac];    // Tirage au hasard du temps de transition pour cet état
    if (Delta_t < dt)
        return n_reac; // L'état change si son "temps" de transition est plus court que celui prévu
    else
        return -1; // Aucun  changement d'état


}



// Calculate a small time step  dt to have Gamma dt ~0.1
// So  in the futur will realize only the reaction having a rate Gamma that verifies
// Gamma dt  *(random) < 1 we do it
// Very efficient for a large number of molecules because typically N_Mol/10 will evolve.
int Fast_Rough_Method_step(const gsl_rng * r, const  vector <double> &rate, double &dt)
{
    int  nb_rate = rate.size();
    double rate_max = 0. ;   // Initialisation pour le taux de transition maximum
    for (int i = 0; i < nb_rate; i++) // Recherche des taux de transition maximum
    {
        if (rate[i] > rate_max)
            rate_max = rate[i];  // Taux de transition maximum
    }

    dt = 0.1/rate_max; // This will be the dt_KMC the 0.1 is arbitrary. Here a tens of the rate will typically be performed

    return 0; // arbitrary (but I put n_reac= 0 in order to avoid problems in rate[n_reac] after). Because several reaction will then be performed in do_reaction
}




// Trouve la réaction (son numéro dans la liste des taux rate[]
// Trouve le temps où elle doit être effectuée
int find_reaction(MC_algorithmes Algorithme_MC , const gsl_rng * r, const vector <double> &rate, double & dt)
{
    int    nb_rate = rate.size();
    int n_reac=0; // numéro de la réaction dans la liste reaction_list[]
    // Fait une évolution N-corps et s'il y a lieu effectue la réaction KMC.
    // Met à jour les intégration temporelle Trate et sum_rate_t et les potentiels et champs
    switch (Algorithme_MC)
    {
    case Aucun_MC : // Il n'y a aucune réaction, on ne sort jamais des boucles
        n_reac = -1;
        dt = VERY_LARGE_NUMBER; // No reaction so infinit time to realize it
        break;


    case Kinetic_Monte_Carlo :
        n_reac = KMC_step(r, rate, dt);
        break;

    case  Random_Selection_Method :
        n_reac = RSM_step(r, rate,  dt);
        break;

    case First_Reaction_Method :
        n_reac = FRM_step(r, rate,  dt);
        break;

    case Fast_Rough_Method :
        n_reac = Fast_Rough_Method_step(r, rate,  dt);

        break;

    }
    return n_reac;
}



