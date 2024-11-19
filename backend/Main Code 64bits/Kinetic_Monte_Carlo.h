<<<<<<< HEAD
/*
  Name: Kinetic Monte Carlo Algorithm
  Author: Daniel Comparat
  Date: 17/09/06 08:14

  Description:
  Implements various algorithms for solving rate equations using Kinetic Monte Carlo (KMC) methods.

  References:
  - "Kinetic Monte Carlo modeling of dipole blockade in Rydberg excitation experiments,"
    Amodsen Chotia, Matthieu Viteau, Thibault Vogt, Daniel Comparat, and Pierre Pillet,
    New Journal of Physics 10 (2008) 045031.

  Algorithms:
  - Kinetic Monte Carlo (KMC)
  - Random Selection Method (RSM)
  - First Reaction Method (FRM)
  - Fast Rough Method (FRM), optimized for speed (evolves approximately N_Mol/10 times faster).

*/

#ifndef KMC_SEEN
#define KMC_SEEN

#include <vector>
#include <iostream>           // For standard I/O (cout, cin)
#include <cmath>              // For mathematical operations (sqrt, etc.)
#include <gsl/gsl_rng.h>      // For random number generator
#include <gsl/gsl_randist.h>  // For Gaussian random number generator
#include "constantes_SI.h"    // SI constants and conversions
#include "algorithmes.h"      // Includes binary search

using namespace std;

// Enumeration of Monte Carlo algorithms
enum MC_algorithmes {
    Aucun_MC = -1,                // No Monte Carlo algorithm
    Kinetic_Monte_Carlo = 0,      // Standard Kinetic Monte Carlo
    Random_Selection_Method = 1, // Random Selection Method
    First_Reaction_Method = 2,   // First Reaction Method
    Fast_Rough_Method = 3        // Fast Rough Method
};

/**
 * Selects an index from a discrete probability distribution defined by `array`.
 * @param r Random number generator.
 * @param array Array of probabilities (rates).
 * @param nb_rate Number of rates in the array.
 * @param sum_rate_tot Sum of all rates (output parameter).
 * @return Index of the selected reaction.
 */
int tirage_random_discrete(const gsl_rng* r, const vector<double>& array, const int nb_rate, double& sum_rate_tot);

/************************************************************************
************************************************************************
 MAIN FUNCTION FOR KINETIC MONTE CARLO (KMC)

 The functions (depending on the chosen algorithm) return the index of the reaction and update the time.
************************************************************************/

/**
 * Kinetic Monte Carlo Algorithm:
 * - Initializes the system in a given state at time t.
 * - Creates a new rate list for the system.
 * - Uses a uniform random number to calculate the time to the first reaction.
 * - Determines the new state based on the calculated time and probabilities.
 *
 * Assumes time-independent rates (most common case).
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int KMC_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
************************************************************************
 FIRST REACTION METHOD (FRM)

 Time-dependent rates:
 - Determines all possible reactions and calculates a time of occurrence for each.
 - Selects the reaction with the smallest time and updates the system state.
 - Repeats for the new configuration.

 Algorithm:
 - For each reaction, calculate:
   ∫_{t_0}^{t_0+Δt_i} rate_i[t] dt = -log(u_i),
   where u_i is a uniform random number.
 - Select the reaction with the smallest Δt_i.

************************************************************************/

/**
 * First Reaction Method:
 * Finds the reaction with the minimum evolution time and updates the time.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int FRM_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
 RANDOM SELECTION METHOD (RSM)

 Selects an atom, its reaction, and a time at random.
 Sometimes no reaction occurs.
************************************************************************/

/**
 * Random Selection Method:
 * Chooses a random reaction and time step.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction (or no reaction).
 */
int RSM_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
 FAST ROUGH METHOD (FRM)

 Optimized method:
 - Calculates the maximum rate and determines the time for evolution as 1/rate_Max.
************************************************************************/

/**
 * Fast Rough Method:
 * Computes the maximum rate and determines the evolution time.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int Fast_Rough_Method_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
 FIND REACTION

 Determines the reaction index and time based on the chosen Monte Carlo algorithm.
************************************************************************/

/**
 * Finds the reaction index and time for execution.
 * @param Algorithme_MC Chosen Monte Carlo algorithm.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int find_reaction(MC_algorithmes Algorithme_MC, const gsl_rng* r, const vector<double>& rate, double& dt);

#endif  // KMC_SEEN
=======
/*
  Name: Kinetic Monte Carlo Algorithm
  Author: Daniel Comparat
  Date: 17/09/06 08:14

  Description:
  Implements various algorithms for solving rate equations using Kinetic Monte Carlo (KMC) methods.

  References:
  - "Kinetic Monte Carlo modeling of dipole blockade in Rydberg excitation experiments,"
    Amodsen Chotia, Matthieu Viteau, Thibault Vogt, Daniel Comparat, and Pierre Pillet,
    New Journal of Physics 10 (2008) 045031.

  Algorithms:
  - Kinetic Monte Carlo (KMC)
  - Random Selection Method (RSM)
  - First Reaction Method (FRM)
  - Fast Rough Method (FRM), optimized for speed (evolves approximately N_Mol/10 times faster).

*/

#ifndef KMC_SEEN
#define KMC_SEEN

#include <vector>
#include <iostream>           // For standard I/O (cout, cin)
#include <cmath>              // For mathematical operations (sqrt, etc.)
#include <gsl/gsl_rng.h>      // For random number generator
#include <gsl/gsl_randist.h>  // For Gaussian random number generator
#include "constantes_SI.h"    // SI constants and conversions
#include "algorithmes.h"      // Includes binary search

using namespace std;

// Enumeration of Monte Carlo algorithms
enum MC_algorithmes {
    Aucun_MC = -1,                // No Monte Carlo algorithm
    Kinetic_Monte_Carlo = 0,      // Standard Kinetic Monte Carlo
    Random_Selection_Method = 1, // Random Selection Method
    First_Reaction_Method = 2,   // First Reaction Method
    Fast_Rough_Method = 3        // Fast Rough Method
};

/**
 * Selects an index from a discrete probability distribution defined by `array`.
 * @param r Random number generator.
 * @param array Array of probabilities (rates).
 * @param nb_rate Number of rates in the array.
 * @param sum_rate_tot Sum of all rates (output parameter).
 * @return Index of the selected reaction.
 */
int tirage_random_discrete(const gsl_rng* r, const vector<double>& array, const int nb_rate, double& sum_rate_tot);

/************************************************************************
************************************************************************
 MAIN FUNCTION FOR KINETIC MONTE CARLO (KMC)

 The functions (depending on the chosen algorithm) return the index of the reaction and update the time.
************************************************************************/

/**
 * Kinetic Monte Carlo Algorithm:
 * - Initializes the system in a given state at time t.
 * - Creates a new rate list for the system.
 * - Uses a uniform random number to calculate the time to the first reaction.
 * - Determines the new state based on the calculated time and probabilities.
 *
 * Assumes time-independent rates (most common case).
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int KMC_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
************************************************************************
 FIRST REACTION METHOD (FRM)

 Time-dependent rates:
 - Determines all possible reactions and calculates a time of occurrence for each.
 - Selects the reaction with the smallest time and updates the system state.
 - Repeats for the new configuration.

 Algorithm:
 - For each reaction, calculate:
   ∫_{t_0}^{t_0+Δt_i} rate_i[t] dt = -log(u_i),
   where u_i is a uniform random number.
 - Select the reaction with the smallest Δt_i.

************************************************************************/

/**
 * First Reaction Method:
 * Finds the reaction with the minimum evolution time and updates the time.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int FRM_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
 RANDOM SELECTION METHOD (RSM)

 Selects an atom, its reaction, and a time at random.
 Sometimes no reaction occurs.
************************************************************************/

/**
 * Random Selection Method:
 * Chooses a random reaction and time step.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction (or no reaction).
 */
int RSM_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
 FAST ROUGH METHOD (FRM)

 Optimized method:
 - Calculates the maximum rate and determines the time for evolution as 1/rate_Max.
************************************************************************/

/**
 * Fast Rough Method:
 * Computes the maximum rate and determines the evolution time.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int Fast_Rough_Method_step(const gsl_rng* r, const vector<double>& rate, double& dt);

/************************************************************************
 FIND REACTION

 Determines the reaction index and time based on the chosen Monte Carlo algorithm.
************************************************************************/

/**
 * Finds the reaction index and time for execution.
 * @param Algorithme_MC Chosen Monte Carlo algorithm.
 * @param r Random number generator.
 * @param rate Array of reaction rates.
 * @param dt Time increment for the next reaction (output parameter).
 * @return Index of the selected reaction.
 */
int find_reaction(MC_algorithmes Algorithme_MC, const gsl_rng* r, const vector<double>& rate, double& dt);

#endif  // KMC_SEEN
>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
