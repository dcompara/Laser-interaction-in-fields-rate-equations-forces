#include "diagonalization.h"
#include "laser.h"                          // for the Euler_angles

#include  <iostream>                       // to include cout, cin


// diagonalized the Hamiltionian for the current molecule its field etc.. and give the eigenvectors and eigenvalues and dipoles (in Debye) update all Level[n].Energy_cm.
void Diagonalization(vector <Internal_state> &Level, const Molecule &my_mol, const Field &fieldB, const Field &fieldE,
                     FitParams &params,  SelfAdjointEigenSolver<MatrixXcd> &es,  MatrixXcd d[])
{
    const int nb_levels=23; //   Level.size();
    // I tried a dynamical size (or vector) but may be not enough and was not easily compatible with the matrix and speed. But should be tried again

    const int nb_incoherent_low_levels=3; //  These are the "dead" levels number 0, 1, ... nb_incoherent_low_levels-1  for annihilation so for incoherent dipole sum
    const int nb_incoherent_high_levels=0; //  These are the "continuum" levels number nb_levels-nb_incoherent_high_levels, nb_levels-2, nb_levels-1  for photoionization so for incoherent dipole sum



    /******* ORDER OF LEVELS (the n=0 and n=1 manifold, the one for spontaneous emission, should be order in Energy) ****************

    The annihilation is treated like a spontaneous emission down to a dead level. BEcause we have states with |m|=0,1,2 we have created 3 dead levels
    with |m|=0,1 name dead_-1 dead_0 dead_1. We have create (abritrary) a pi spontaneous emission with all m=-1,0,1 and only the 3P2-2   and 3P22 are simga+ and sigma-

    For convenience we have ordered the level in energy such that the Zeeman effect keep this ordering (for very small field)
    And for degenerate levels that stay degenerate even with Zeeman effect we have added a small shift to ensure the energy ordering

    The matrix are calculated using a Mathematica code which use n S L J M_J ordering that is the most natural orderning. So we use it here.


    However the diagonalization then gives an ordering in energy that is

    dead_-1 dead_0 dead_1    n=1[  1S00     3S1-1    3S11    3S10],   N=2[  1S00    3P00   3P1-1    3P11     3P10    1P10   1P1-1    1P11   3P2-2    3P22    3P2-1   3P21   3P20     3S1-1   3S11   3S10
    [i] 0        1       2           3         4        5       6              7       8       9       10      11     12     13       14      15      16     17      18      19      20       21     22

    *************************************/


    double E0_cm[nb_levels][nb_levels] = // With 1e-5 offset
    {
        {-10000.00001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, -10000.00000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, -9999.99999, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, -6.81758, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, -1e-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1e-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 41147.81261, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 41148.05609, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.23870, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.23871, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.23872, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.29957, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.29958, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.29959, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.38478, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.38479, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.38480, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.38481, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.38482, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.66480, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.66481, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.66482}
    };

    // With 1e-3 offset
//    {
//        {-10000.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, -10000.000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, -9999.999, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, -6.817, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, -1e-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 1e-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 41147.812, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 41148.056, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.237, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.238, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.239, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.298, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.299, 0, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.300, 0, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.382, 0, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.383, 0, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.384, 0, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.385, 0, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.386, 0, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.663, 0, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.664, 0},
//        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41148.665}
//    };

    /** For the Zeeman effect in cm^-1  ****/

    // Ruggero = Dermer = mine = Pauline (if (-1)^(l+s-j)(-1)^(lp+sp-jp) included)
    double Zeeman_cm_B[nb_levels][nb_levels] =
    {
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0.9337307964640155, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0.9337307964640155, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.9337307964640155},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.5390897266891428, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.6602473779824212, 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.6602473779824212, 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., -0.5390897266891428, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.7623880028197907, 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., -0.6602473779824212, 0., 0., 0., 0., 0., 0., 0., 0.6602473779824212, 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.6602473779824212, 0., 0., 0., 0., 0., 0., 0., 0.6602473779824212, 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.6602473779824212, 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.6602473779824212, 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.7623880028197907, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
        {0., 0., 0., 0., 0., 0., 0., 0.9337307964640155, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}
    };


    /***  <i|d_q|j> in Debye is coded here as d[q+1][i][j] (real), i = line, j = column ; that is for a i<-->j transition (with E_i> E_j) ****/

    double dipole[3][nb_levels][nb_levels]=
    {
        {   {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.007877448760964387, 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3.7868764758423907, 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., -2.1863541527154653, 0., 0., -2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., -1.5459858474604744, 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3.7868764758423907, 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., -2.6777260355839694, 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -15.250486562036826, 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 2.1863541527154653, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 8.80487252186473, 0.},
            {0., 0., 0., 0., 0., 0., -2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -10.783722464410557},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., -2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -10.783722464410557, 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 3.7868764758423907, 0., 0., 0., 15.250486562036826, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0.007877448760964387, 0., 0., 0., 3.7868764758423907, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 15.250486562036826, 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 1.5459858474604744, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 6.2259850676936495, 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., -8.80487252186473, 0., 0., -10.783722464410557, 0., 0., 0., 0., 0., 0., 0., -6.2259850676936495, 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -15.250486562036826, 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -10.783722464410557, 0., 0., 0., 0., 0., 0., 0., -10.783722464410557, 0., 0., 0., 0.}
        },
        {   {0., 0., 0., 0., 4.7220552656855155, 0., 0., 0., 0., 0.00048814973718447113, 0., 0., 0., 0.0026750393062167006, 0., 0., 0., 0.007877448760964387, 0., 0., 0.14482997883476675, 0., 0.},
            {0., 0., 0., 159.87787069042116, 0., 0., 4.7220552656855155, 4.88155837172815, 0.015436732757477556, 0., 0., 0.00048814973718447113, 0.0026750393062167006, 0., 0., 0., 0., 0., 0., 0.007877448760964387, 0., 0., 0.14482997883476675},
            {0., 0., 0., 0., 0., 4.7220552656855155, 0., 0., 0., 0., 0.00048814973718447113, 0., 0., 0., 0.0026750393062167006, 0., 0., 0., 0.007877448760964387, 0., 0., 0.14482997883476675, 0.},
            {0., 159.87787069042116, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3.7868764758423907, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {4.7220552656855155, 0., 0., 0., 0., 0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0.},
            {0., 0., 4.7220552656855155, 0., 0., 0., 0., 0., 0., 0., -2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0.},
            {0., 4.7220552656855155, 0., 0., 0., 0., 0., 0., -2.1863541527154653, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 3.091971694920949, 0., 0., 0.},
            {0., 4.88155837172815, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 15.250486562036826, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0.015436732757477556, 0., 0., 0., 0., -2.1863541527154653, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -8.80487252186473},
            {0.00048814973718447113, 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0., 0.},
            {0., 0., 0.00048814973718447113, 0., 0., -2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -10.783722464410557, 0.},
            {0., 0.00048814973718447113, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0.0026750393062167006, 0., 3.7868764758423907, 0., 0., 0., 15.250486562036826, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0.0026750393062167006, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0.0026750393062167006, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0.007877448760964387, 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0., 0.},
            {0., 0., 0.007877448760964387, 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0.},
            {0., 0.007877448760964387, 0., 0., 0., 0., 3.091971694920949, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 12.451970135387299},
            {0.14482997883476675, 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0., 0., 0., 0., 0.},
            {0., 0., 0.14482997883476675, 0., 0., 0., 0., 0., 0., 0., -10.783722464410557, 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0., 0., 0., 0.},
            {0., 0.14482997883476675, 0., 0., 0., 0., 0., 0., -8.80487252186473, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 12.451970135387299, 0., 0., 0.}
        }
        ,
        {
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -0.007877448760964387, 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3.7868764758423907, 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -3.7868764758423907, 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., -2.1863541527154653, 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., -1.5459858474604744, 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., -2.6777260355839694, 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -15.250486562036826, 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 2.1863541527154653, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 8.80487252186473, 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557},
            {0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 3.7868764758423907, 0., 0., 0., 15.250486562036826, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0.007877448760964387, 0., 0., 3.7868764758423907, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 15.250486562036826, 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 2.6777260355839694, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557},
            {0., 0., 0., 0., 1.5459858474604744, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 6.2259850676936495, 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -15.250486562036826, 0., 0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., -8.80487252186473, 0., 0., 10.783722464410557, 0., 0., 0., 0., 0., 0., 0., -6.2259850676936495, 0., 0., 0.},
            {0., 0., 0., 0., 0., 0., 0., 0., 0., 10.783722464410557, 0., 0., 0., 0., 0., 0., 0., -10.783722464410557, 0., 0., 0., 0., 0.}
        }
    };

// For small sizes, especially for sizes smaller than (roughly) 16, using fixed sizes is hugely beneficial to performance, as it allows Eigen to avoid dynamic
// So here I use fixed size
// Furthermore it seems that using dynamical and MatrixXcd (not with MatrixXd for H) create a memory problem!!

    /*** dipole matrix element ***/
    MatrixXcd d0[3]; // d0 = matrix dipole in zero field
    d0[0] = MatrixXcd(nb_levels,nb_levels);
    d0[1] = MatrixXcd(nb_levels,nb_levels);
    d0[2] = MatrixXcd(nb_levels,nb_levels);
    for (int i=0; i<nb_levels; i++)
        for (int j=0; j<nb_levels; j++)
        {
            for(int n_polar = -1; n_polar <= 1; n_polar++)
            {
                d0[n_polar+1](i,j)= dipole[n_polar+1][i][j]; // d0[q+1]_ij = 0_<i | d_q | j>_0   (i = out and j = in)
                // cout << "polar " << n_polar << " i  "  << i  << " j " << j  << " d " <<  d0[n_polar+1](i,j) << endl;
            }
        }


    /*** Hamiltonian matrix element: Zeeman + Stark (External field + Dynamical Stark) ***/
    MatrixXcd H(nb_levels,nb_levels); // Hamiltonian Matrix. It is an hermitian matrix so I use complex not MatrixXcd

    Vecteur3D r,v,B,F,dip_ij;
    r = my_mol.get_pos();
    v = my_mol.get_vel();
    B = fieldB.get_Field(r); // I neglect the dynamical modification of B field for v << c cf https://en.wikipedia.org/wiki/Classical_electromagnetism_and_special_relativity
    F = fieldE.get_Field(r) + v.cross(B); // External field + Dynamical Stark, noted F to avoid confusion with E = Energy
    complex<double> d_dot_F;

    // Hamiltonian. We assume here an adiabatic following of the state and so that their quantization axis is the B axis
    for (int i=0; i<nb_levels; i++)
        for (int j=0; j<nb_levels; j++)
        {
            dip_ij = Vecteur3D(dipole[0][i][j],dipole[1][i][j],dipole[2][i][j])*Debye/(100.*hPlanck*C); // (to have the dipole such that d.F in cm^-1)
            dipole_dot_Electric_Ffield(dip_ij, B, F, d_dot_F); // Hamiltonian Stark = - d.F = +e r.F
            /**** SIGN TO BE CHECKED - d.F = +e r.F ***/
            H(i,j)=  E0_cm[i][j] + B.mag()*Zeeman_cm_B[i][j] + d_dot_F; // H_ij = 0_<i | H | j>_0 = 0_<j | H | i>_0  (cf Eq (3) of Dermer PRA 40, 5526 (1989)

//           cout << " i,j " << i <<  " " << j << "   " << E0_cm[i][j] << " " << Zeeman_cm_B[i][j] << "  "  << Stark_cm_Bv[i][j] << endl;
        }

    /*** diagonalization: gives new eigen energies (in Level[n].Energy_cm) and eigen vectors  ***/

    es.compute(H);  // calculate the new eigenvectors |i> (i start from 0)  and new eingen_Energies
    // E_i = es.eigenvalues()(i) (in incresing order)
    // es.eigenvectors()(i,j) = 0<i | j>  gives (i=line, j = column index) the new (column) vector |j> in function of the old |i>_0



//
//    for( int j0 = 0; j0 < nb_levels; j0++ ) // Ligne
//    {
//        for( int j = 0; j < nb_levels; j++ )
//        {
//            // cout <<  es.eigenvectors()(j0,j).real() << " "; // 0<j0 | j>
//
//            cout <<  (round(100.* es.eigenvectors()(j0,j).real()))/100  << "+i" << (round(100.* es.eigenvectors()(j0,j).imag()))/100 << "   ";
//        }
//        cout << endl;
//    }
//    cout << endl << endl;





    for( int n = 0; n < nb_levels; n++ )
    {
        Level[n].Energy_cm = es.eigenvalues()(n);
    }


    /**** calcul of the new dipoles (in Debye) in d[q]

    evec = es.eigenvectors() verifie evec(j0,j) = 0<j0 | j>  gives (j0=line, j = column index) the new (column) vector |j> in function of the old |j0>_0
    d0[q+1]_i0 j0 = 0_<i0 | d_q | j0>_0 . So d[q+1]_ij = <i | d_q | j> = Sum i0,j0   <i |i0>0 0<i0| d_q | j0>0 0<j0|j> =  Sum i0,j0  evec^dag (i,i0)  d_q(i0,j0)  evec(j0,j)

    IN CONCLUSION: the new dipole are given by d[polar] = evec^dag.d0[polar].evec =  <i | d_q | j> with evec_j = |j> = sum_|j>_0   0_<j| j>.


    But for some transition (low to dead levels and high to continuum we use an incoherent sum. To do this we use the Hadamar product A�B = (A.cwiseProduct(B)) defined by (A�B)__ij = A_ij B_ij
    d_incoh[q+1]_ij^2 = <i | d_incoh_q | j> = Sum i0,j0   | <i |i0>0 0<i0| d_q | j0>0 0<j0|j> |^2 =  Sum i0,j0  |evec^dag (i,i0)  d_q(i0,j0)  evec(j0,j)|^2 = Sum i0,j0  evec^dag (i,i0) evec^dag(i,i0)^*   d_q(i0,j0) d_q(i0,j0)^*  evec(j0,j) evec(j0,j)^*
    d_incoh[q+1]_ij^2 = Sum i0,j0  evec^dag�evec^dag* (i,i0)    d_q�d_q* (i0,j0)  evec�evec* (j0,j)  =  (evec^dag�evec^dag* .   d_q�d_q* .  evec�evec*)_ij
    with for such state typically the |j> levels does not change so |j>=|j0>_0 and the last sum is Sum i0  |<i |i0>0 0<i0| d_q | j=j0>0 |^2
    But we keep in case of order changing in the levels files ...

    IN CONCLUSION dSQUARE_incoh = evec^dag�evec^dag* .   d_q�d_q* .  evec�evec*
    ***/

    for(int n_polar = -1; n_polar <= 1; n_polar++)
    {
        // We first calculate all dipoles as incoherent sum and then we replace the central part with the proper coherent dipoles. Remark: We do twice the work (so it is slower, but what is long is probably the diagonalization) but is is very simple to write

        // incoherent sum for all dipoles
        d[n_polar+1] = ( ( es.eigenvectors().adjoint()).cwiseProduct(es.eigenvectors().adjoint().conjugate()) ) * ( d0[n_polar+1].cwiseProduct(d0[n_polar+1].conjugate()) ) * ( (es.eigenvectors()).cwiseProduct(es.eigenvectors().conjugate()) );
        // this is not yet the matrix of the dipole but the one of dipole squared d_incoh[q+1]_ij^2
        for (int i=0; i<nb_levels; i++)
            for (int j=0; j<nb_levels; j++)
                d[n_polar+1](i,j) = sqrt(d[n_polar+1](i,j)); // To go back to d_incoh[q+1]_ij

        // We replace the central part with the proper block matrix of coherent dipole. Using matrix.block<p,q>(i,j) = Block of size (p,q), starting at (i,j)
        const int i = nb_incoherent_low_levels;
        const int p = nb_levels - nb_incoherent_high_levels - nb_incoherent_low_levels;  // size of coherent block
        d[n_polar+1].block<p,p>(i,i) = ( (es.eigenvectors().adjoint())*d0[n_polar+1]*(es.eigenvectors()) ).block<p,p>(i,i);
    }

    H.resize(0,0);
    d0[0].resize(0,0);
    d0[1].resize(0,0);
    d0[2].resize(0,0);
}


//   d.E (for the Stark effect that is -d.E)
// where the dipole  vector d= sum_p d_p e^p is given in the local quantification axis (typically the magnetic field)
void dipole_dot_Electric_Ffield (const Vecteur3D& dipole, const Vecteur3D& axe_quant,  const Vecteur3D& E, complex<double> &d_dot_E)
{
    double dp,d0,dm;
    dm = dipole(0);
    d0 = dipole(1);
    dp = dipole(2);

    Vecteur3D Euler_angles_axe_quant = Euler_angles(axe_quant);
    Vecteur3D Euler_angles_axe_electric_field = Euler_angles(E);

// The link between the polar angles (theta, phi) and the Euler angle (alpha, beta, gamma) in ZXZ convention as we used them (Wikipedia)  are
    // alpha= phi +pi/2; beta = theta; gamma = -pi/2
    double phi_F = Euler_angles_axe_quant(0) - pi/2.; // polar angle for the quantization axis
    double theta_F = Euler_angles_axe_quant(1); // polar angle for the quantization axis
    double phi_k = Euler_angles_axe_electric_field(0) - pi/2.; // polar angle for the electric field axis (along the vector E)
    double theta_k = Euler_angles_axe_electric_field(1); // polar angle for the electric field axis (along the vector E)


    complex<double> twice_dip; // Twice the dipole
    const complex<double> i(0., 1.);

    /*** We try to optimize the calcul so to reduce at maximmum the numebr of operations ****/

    double sin_theta_k = sin(theta_k); // The notation k is used here because this is the notation for a field as if it was a laser propagating along k with a polarization along k by a laser (cf effectif_dipole_local() function)
    double sin_theta_F = sin(theta_F);
    double cos_theta_k = cos(theta_k);
    double cos_theta_F = cos(theta_F);
    double cos_phi_F_minus_phi_k = cos(phi_F-phi_k);
    double sin_phi_F_minus_phi_k = sin(phi_F-phi_k);

    double sqrt2 = sqrt(2.);

    double dp_minus_dm = dp-dm;
    double dp_plus_dm = dp+dm;

    twice_dip  = cos_theta_k*(2.*d0*cos_theta_F + sqrt2*dp_minus_dm*sin_theta_F)+
                 sin_theta_k*(-sqrt2*dp_minus_dm*cos_theta_F*cos_phi_F_minus_phi_k+
                              2.*d0*cos_phi_F_minus_phi_k*sin_theta_F - i*sqrt2*dp_plus_dm*sin_phi_F_minus_phi_k);

    d_dot_E = E.mag() * twice_dip/2.;
}

