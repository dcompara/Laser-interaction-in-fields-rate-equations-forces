/*********************
Function sign
Binary search
Wait in seconds
*********************/

#ifndef Algorithmes_SEEN
#define Algorithmes_SEEN

#include <time.h>          // For clock()
#include <fstream>         // For reading data and writing to files
#include <iostream>
#include <math.h>          // For floor()
using namespace std;

// Function to get the sign of a number
template <typename T>
inline int sgn(T t)
{
    return t > 0 ? 1 : t < 0 ? -1 : 0;
}

/*
Binary search function using dichotomy.
Returns position such that sortedArray[position] < key <= sortedArray[position+1].
Returns -1 if nothing is found.

Source: http://www.fredosaurus.com/notes-cpp/algorithms/searching/binarysearch.html
See also #include <algorithm>  // For binary search
Also exists in GSL as gsl_histogram_find.
lower_bound exists as well.
*/
int binarySearch(double sortedArray[], int first, int last, double key);

// Returns the number of lines in a file
int number_line_file(const char *filename);

// Waits for a specified number of seconds (accurate to milliseconds)
void wait(double sec);

// Compares doubles for qsort
// If they match in ranking, the function shall return zero;
// if elem1 goes before elem2, it shall return a negative value;
// if it goes after, a positive value.
int compare_doubles(const void *x, const void *y);

#endif

