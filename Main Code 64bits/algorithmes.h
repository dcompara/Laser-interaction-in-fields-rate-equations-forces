/*********************
fonction signe
binary search
wait en seconde
*********************/

#ifndef Algorithmes_SEEN
#define Algorithmes_SEEN

#include <time.h>                // For clock()
#include  <fstream>                        // to read the data and put them in files
#include <iostream>
#include <math.h>       /* floor */
using namespace std;


// Fonction signe
template <typename T> inline int sgn(T t)
{
    return t > 0 ? 1 : t < 0 ? -1 : 0;
}


// Macro to Switch laser (based on the rule 0 = false ; 1 =  true) Gives
// 0 if t is between 0 and T1 (modulo T1+T2)
// 1 if t is between T1, T1+T2 (modulo T1+T2)
// Thus we can switch between lasers by nlas + Nb_laser * Is_Switch
int Is_Switch(double T1, double T2, double t);


/*
Fonction de recherche binaire par dichotomie
Retourne position tel que sortedArray[position] < key <= sortedArray[position+1]
-1 si rien n'est trouvé

Provient de http://www.fredosaurus.com/notes-cpp/algorithms/searching/binarysearch.html
Voir aussi  #include <algorithm>                     // Pour le binary search
Existe aussi dans gsl c'est gsl_histogram_find
lower_bound existe aussi
*/

int binarySearch(double sortedArray[], int  first, int  last, double key);

// Retourne le nb de lignes d'un fichier
int number_line_file(const char *nom_file);

// Attend le nb de seconds (précis à la milliseconde)
void wait(double sec);

// Compare des doubles pour qsort
// if they match in ranking, the function shall return zero; if elem1 goes before elem2, it shall return a negative value; and if it goes after, a positive value.
int compare_doubles(const void *x, const void *y);


#endif


