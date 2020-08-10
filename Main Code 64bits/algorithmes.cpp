/*********************
Modulo
binary searchr
File number of line
wait en seconde
*********************/

#include "algorithmes.h"



// Macro to Switch laser (based on the rule 0 = false ; 1 =  true) Gives
// 0 if t is between 0 and T1 (modulo T1+T2)
// 1 if t is between T1, T1+T2 (modulo T1+T2)
// Thus we can switch between lasers by nlas + Nb_laser * Is_Switch
int Is_Switch(double T1, double T2, double t)
{
    double reste = t - (T1 + T2) * floor(t/(T1 + T2)); // is t modulo T1+T2
    if (reste < T1)
        return 0;
    else
        return 1;
}



/*
Fonction de recherche binaire par dichotomie
Retourne position tel que sortedArray[position] < key <= sortedArray[position+1]
-1 si rien n'est trouvé

Provient de http://www.fredosaurus.com/notes-cpp/algorithms/searching/binarysearch.html
Voir aussi  #include <algorithm>                     // Pour le binary search
Existe aussi dans gsl c'est gsl_histogram_find
lower_bound existe aussi
*/
int binarySearch(double sortedArray[], int  first, int  last, double key)
{
    int llast = last ;
    int ffirst = first;
    int mid;
    while (ffirst <= llast)
    {
        mid = (int) (ffirst + llast) / 2;  // compute mid point.
        if (key > sortedArray[mid])
            ffirst = mid + 1;  // repeat search in top half.
        else if (key < sortedArray[mid])
            llast = mid - 1; // repeat search in bottom half.
        else
            return mid; // Trouver la position exacte key = sortedArray[mid]
    }
    return llast; // sortedArray[llast] < key < sortedArray[first] =  sortedArray[llast+1]
}



// Retourne le nb de lignes d'un fichier
int number_line_file(const char *nom_file)
{
    ifstream file(nom_file);
    int n = 0;
    string s;
    while( getline( file, s ) )
    {
        n++;
    }
    file.close();
    return n;
}



// Attend le nb de seconds (précis à la milliseconde)
void wait(double sec)
{
    clock_t start;
    start = clock();


    while ((clock()-start)/double(CLOCKS_PER_SEC) < sec)
    {
    }
}


// Compare des doubles pour qsort
// if they match in ranking, the function shall return zero; if elem1 goes before elem2, it shall return a negative value; and if it goes after, a positive value.
int compare_doubles (const void *x, const void *y)
{

    // x and y are pointers to doubles.

    // Returns -1 if x < y
    //          0 if x == y
    //         +1 if x > y

    double dx, dy;

    dx = *(double *)x;
    dy = *(double *)y;

    if (dx < dy)
    {
        return -1;
    }
    else if (dx > dy)
    {
        return +1;
    }
    return 0;
}
