/*********************
Modulo
Binary search
File number of lines
Wait in seconds
*********************/

#include "algorithmes.h"
#include <chrono>
#include <thread>

/**
 * Binary search function using dichotomy.
 *
 * @param sortedArray - Array of doubles, must be sorted in ascending order.
 * @param first - Index of the first element to consider.
 * @param last - Index of the last element to consider.
 * @param key - The value to search for.
 * @return Index `position` such that sortedArray[position] < key <= sortedArray[position+1],
 *         or -1 if not found.
 */
int binarySearch(double sortedArray[], int first, int last, double key)
{
    int llast = last;
    int ffirst = first;
    int mid;
    while (ffirst <= llast)
    {
        mid = (ffirst + llast) / 2;  // Compute mid-point.
        if (key > sortedArray[mid])
            ffirst = mid + 1;  // Search in the top half.
        else if (key < sortedArray[mid])
            llast = mid - 1;  // Search in the bottom half.
        else
            return mid;  // Exact position found.
    }
    return llast;  // sortedArray[llast] < key < sortedArray[llast + 1].
}

/**
 * Counts the number of lines in a given file.
 *
 * @param filename - Path to the file.
 * @return Number of lines in the file.
 */
int number_line_file(const char *filename)
{
    std::ifstream file(filename);
    int n = 0;
    std::string s;
    while (std::getline(file, s))
    {
        n++;
    }
    return n;
}

/**
 * Pauses execution for the specified number of seconds.
 * Accurate to milliseconds.
 *
 * @param sec - Duration to wait, in seconds.
 */
void wait(double sec)
{
    std::this_thread::sleep_for(std::chrono::duration<double>(sec));
}

/**
 * Comparison function for doubles used in qsort.
 *
 * @param x - Pointer to the first element.
 * @param y - Pointer to the second element.
 * @return Negative if x < y, zero if x == y, positive if x > y.
 */
int compare_doubles(const void *x, const void *y)
{
    double dx = *(const double *)x;
    double dy = *(const double *)y;

    return (dx < dy) ? -1 : (dx > dy) ? 1 : 0;
}
