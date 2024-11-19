/*
  Name: Class "Field"
  Copyright:
  Author: Daniel Comparat
  Date: 15/9/2012

Class Field.

A field \( F \) is defined by analytical formulas in the vector form \( F = F_0 + F_1 + F_2 + \dots + F_n \) for a given \( n \).
For the electric field, we also include the electric potential.

Initialization starts at zero.

There are several types of possible fields, as defined by `type_field_read`:

0: Defined up to the second order: Fields along x, y, and z are decomposed by component. Example for \( F_x \): \( F_0x + F_1x \cdot x + F_2x \cdot x^2 + \dots + F_nx \cdot x^n \).
1: Fields in Helmholtz coils (typically used with magnetic fields).
2: 3D field maps derived from 2D cylindrical symmetry: 4 columns (r, z, \( F_r(r, z) \), \( F_z(r, z) \)).
3: 3D field maps from 2D cylindrical symmetry, including first and second derivatives (\( \frac{\partial F_r}{\partial r} \), \( \frac{\partial F_z}{\partial z} \), and \( \frac{\partial^2}{\partial r \partial z} \)).
4: 3D field maps with six components: x, y, z, \( F_x \), \( F_y \), \( F_z \). (TO BE IMPLEMENTED)
*/

#ifndef Field_SEEN
#define Field_SEEN

#include "Vecteur3D.h"
#include <vector>
#include "BicubicInterpolation_module.h"
#include "TricubicInterpolation_module.h"

using namespace std;


enum // Différents types de champs
{
    Magnetic_Field = 0,
    Electric_Field = 10
};


class Field
{
protected:
    Vecteur3D F0;             // Bias (zero order Field)
    Vecteur3D F1;             // Gradient (first order)
    Vecteur3D F2;             // quadratic term
    Vecteur3D Fn;             // power n term


public:
    int n;                     //power for the last term. Will be 3 by default
    int type_field; //  Magnetic_Field or Electric_Field
    int type_field_read;

// DEFAULT: 0: 2nd order plus a nth order
// 1:  Field in Helmoltz coils (so usualy goes with a field that is magnetic)
// 2: Field map 3D from 2D cylindrical symmetry: 4 columns r,z, F_r(r,z); F_z(r,z)
// 3: Field map 3D from 2D cylindrical symmetry: F_r(r,z); F_z(r,z)+ derivative d/dr; d/dz and d^2/drdz
// 4: Field map 3D: 6x,y,z, Bx, By, Bz (TO BE DONE)

    int Nb_bobines;         // Nombre de bobines
    double gap_bobines;         // Ecart entre les bobines
    double courant_bobines;         // Courant dans les bobines
    double rayon_bobines;         // rayon les bobines
    int is_Helmholtz ;      // Si true on ajoute à chaque bobine la même décalée du rayon
    double potential; // 0 for magnetic (may be used for electric potential)
    vector< vector<double> > matrix_field; // gives the field on a regular grid. Depend on the type_field_read. This table can be  r z Br Bz or x y z Bx By Bz or with derivative etc ... NO MORE THAN 100 columns!
    // **matrix  would be another choice cf http://www.cplusplus.com/forum/articles/7459/

    // BE CAREFUL matrix_field[j][i]  is the element [i,j] of the line i, column j of the file
    // THUS matrix_field[j] is the column number j of the file.

    Bicubic_Interpolator Interpolator2D; // Class Bicubic Interpolator. Contains matrix size, spacing in grids and algorithm for interpolation
    // LM_Interpolator Interpolator3D; // Class Triubic Interpolator

    // Forme canonique d'une classe
    Field();                        //constructeur
    Field(const Field&);            // Constructeur de copie
    Field& operator = (const Field&);       // Affectation par recopie
    virtual ~Field() {};                       // Destructeur par defaut


    // Get components
    Vecteur3D  get_F0()  const
    {
        return F0;
    }
    Vecteur3D  get_F1()  const
    {
        return F1;
    }
    Vecteur3D  get_F2()  const
    {
        return F2;
    }
    int  get_n()  const
    {
        return n;
    }

    Vecteur3D  get_Fn()  const
    {
        return Fn;
    }



    // Set components
    void  set_F0(const Vecteur3D& new_F0)
    {
        F0 = new_F0;
    }

    void  set_F1(const Vecteur3D& new_F1)
    {
        F1 = new_F1;
    }
    void  set_F2(const Vecteur3D& new_F2)
    {
        F2 = new_F2;
    }
    void  set_n(const int& new_n)
    {
        n = new_n;
    }
    void  set_Fn(const Vecteur3D& new_Fn)
    {
        Fn = new_Fn;
    }



    friend Field operator +(const Field , const Field);    // surcharge de l'opérateur +


// I do not update it. So it is incomplete!
    virtual void write(ostream & flux)
    {
        if (&flux == &cout)
            cout << "Bias : " << F0 << "\t";
        else
            flux << F0 << "\t";
        if (&flux == &cout)
            cout << "Grad : " << F1 << "\t";
        else
            flux << F1 << "\t";
        if (&flux == &cout)
            cout << "Grad Grad : " << F2 << "\t";
        else
            flux << F2 << "\t";
        if (&flux == &cout)
            cout << "potential : " << potential << "\t";
        else
            flux << potential << "\t";

    }

    // read data
    // Same incomplete

    virtual void read(istream & flux)
    {
        if (&flux == &cin)
            cout << "entrez Bias field : ";
        flux >> F0;
        if (&flux == &cin)
            cout << "entrez Grad : ";
        flux >> F1;
        if (&flux == &cin)
            cout << "entrez Grad Grad : ";
        flux >> F2;
        if (&flux == &cin)
            cout << "entrez potential : ";
        flux >> potential;
    }


// Get field at a given position from the grid
    Vecteur3D  get_Field_interpolation2D(const Vecteur3D pos, int type_field_read = 0) const;

    // Get Grad (F.F) = gradient of the norm of the field squared at a given position from the grid
    // Also give the field F
    void  get_norm2_grad_Field_interpolation2D(const Vecteur3D pos, Vecteur3D &field, Vecteur3D &grad_F2, int type_field_read = 0) const;


// Get field at a given position
    Vecteur3D  get_Field(const Vecteur3D &pos) const;


// Get potential at a given position
    // Be careful give the potential as if it was an electric potential (so it integrates the field)
    double  get_electric_potential(const Vecteur3D pos)  const;

    // Calcul du gradient de la norme (au carré) du champ local F
// Grad (F.F) = 2 (F.grad) F
    Vecteur3D get_grad_field_F2(const Vecteur3D& pos ) const;



    /***************  BOBINE  *****************/


// CALCUL APPROCHE. Les formules analytiques existes mais sont complexes.
// Ici ce sont des formules plus simples mais très proche de la réalité (cf bz2, brho2 Mathematica ficher calcul_Sisyphus.nb
// Champ   pour bobine rayon r centrée en z=0 axe Oz avec courant de I_coil ampère
    // ATTENTION ICI LES CALCULS SONT fait en supposant les bobines d'axes Oz alors qu'elle seront d'axe Ox en réalité.
    // Il faudra donc appeler ces fonctions non pas en (x,y,z) mais en (y,z,x)

// Pour info le champ au centre est Mu I0/ (2 r)
    // Champ crée par une bobine positionée en pos_bobine


    Vecteur3D get_field_coil(const Vecteur3D pos, double r, double I_coil, Vecteur3D pob_bobine, double &Brho_sur_rho, double &Bz ) const;
    Vecteur3D get_field_coil(const Vecteur3D pos, double r, double I_coil, Vecteur3D pob_bobine ) const;

// Champ pour Nb_coils bobines de rayon r la première centrée en z=0 les autres séparées de  z_gap
    Vecteur3D get_field_coils(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz, double &Brho_sur_rho, double &Bz) const;
    Vecteur3D get_field_coils(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz) const;


// Gradient du champ   pour une bobines rayon r centrée en z=0 axe Oz
// En fait donne   dBz_drho, dBz_dz,  dBrho_drho, dBrho_dz
    void get_grad_field_coil_4(const Vecteur3D pos, double r, double I_coil, double &dBz_drho_sur_rho, double &dBz_dz,  double &dBrho_drho, double &dBrho_dz_sur_rho) const;
    // Gradient du champ   pour des bobines rayon r axe Oz
    // En fait donne   dBz_drho, dBz_dz,  dBrho_drho, dBrho_dz
    void get_grad_field_coils_4(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz, double &dBz_drho_sur_rho, double &dBz_dz,  double &dBrho_drho, double &dBrho_dz_sur_rho) const;


// Gradient de la Norme  du champ   pour des bobines rayon r centrée en z=0 axe Oz
    Vecteur3D get_grad_normfield2_coils(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz) const;


// Erase and Read the matrix of discrete value from a file and return the number of elements
// The  file depends on the type of field but should be ordered x By Bz.
// It should be first along the n1 axis [0,0,0], [1,0,0], ..., [n1-1,0,0], [0,1,0], ...
// and  equally spaced along each axis,
    int Read_Matrix( const char * Fieldfilename);

// Calculation and output of partial derivation of the field
// For instance for the type_field 2 (that does not contains the derivatives. I calculate them)
    void Calculate_Derivative_Matrix(int number_matrix_lines, const char * Fieldfilename) const;


// Using the Matrix of the field (that should have been read before!) we extract the grid parameters and initialize the Interpolators depending on the type of field used
    void Init_Interpolator();

};

ostream& operator << (ostream &flux, Field field);

istream& operator >> (istream &flux,  Field & field); //at est modifié!


#endif
