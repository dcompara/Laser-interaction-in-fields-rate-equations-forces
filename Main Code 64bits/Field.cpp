#include "Field.h"
#include <math.h>
#include "constantes_SI.h"
#include "algorithmes.h"
// ------------
// Constructors
// ------------


Field::Field()                       // Constructeur par défaut
{
    F0.set(0.,0.,0.);
    F1.set(0.,0.,0.);
    F2.set(0.,0.,0.);
    Fn.set(0.,0.,0.);
    n = 3;
    type_field = Magnetic_Field; // By default we look for Zeemnan effect
    type_field_read=0;   // Par défaut le champ est donné par F0,F1 et F2
    Nb_bobines=0;
    gap_bobines=0.;
    courant_bobines = 0.;
    rayon_bobines = 1.;         //Pour éviter des divisions par zéro (mais comme le courant est nul cela importe peu)
    is_Helmholtz = 0;
    potential = 0.;
    // We start with empty matrix memory
    matrix_field =  vector< vector<double> > ();

};

//Field::~Field()   // Destructeur
//{
//    matrix_field.clear();
//}
//;

// NOT CORRECT BECAUSE THE MATRIX IS NOT COPY
Field::Field(const Field & field)               // Constructeur de (re)copie
{
    F0 = field.F0;
    F1 = field.F1;
    F2 = field.F2;
    n = field.n;
    Fn = field.Fn;
    type_field = field.type_field;
    type_field_read=field.type_field_read;
    Nb_bobines=field.Nb_bobines;
    gap_bobines=field.gap_bobines;
    courant_bobines=field.courant_bobines;
    rayon_bobines=field.rayon_bobines;
    is_Helmholtz = field.is_Helmholtz;
    potential = field.potential;

    // Bicubic_Interpolator_r_z = field.Bicubic_Interpolator_r_z; // Does not exist so !!
};

Field & Field::operator = (const Field & field)          // Affectation par recopie
{
    if (this != &field) // On vérifie que les objets ne sont pas les mêmes !
    {
        F0 = field.F0;
        F1 = field.F1;
        F2 = field.F2;
        n = field.n;
        Fn = field.Fn;
        type_field = field.type_field;
        type_field_read=field.type_field_read;
        Nb_bobines=field.Nb_bobines;
        gap_bobines=field.gap_bobines;
        courant_bobines=field.courant_bobines;
        rayon_bobines=field.rayon_bobines;
        is_Helmholtz = field.is_Helmholtz;
        potential = field.potential;
    }
    return *this;
}

//--------------------------------------------
// Surcharge des opérateurs +,-,*,/
// + (-) = addition (soustraction) membre à membre.
//--------------------------------------------

Field operator +(const Field field1, const Field field2)
{
    Field sum;

    sum.set_F0(field1.get_F0()+field2.get_F0());
    sum.set_F1(field1.get_F1()+field2.get_F1());
    sum.set_F2(field1.get_F2()+field2.get_F2());
    sum.set_n(field1.get_n()+field2.get_n());
    sum.set_Fn(field1.get_Fn()+field2.get_Fn());
    sum.type_field = field1.type_field + field2.type_field;
    sum.type_field_read = field1.type_field_read + field2.type_field_read;
    sum.Nb_bobines = field1.Nb_bobines + field2.Nb_bobines;
    sum.gap_bobines= field1.gap_bobines + field2.gap_bobines;
    sum.courant_bobines = field1.courant_bobines + field1.courant_bobines;
    sum.rayon_bobines = field1.rayon_bobines + field2.rayon_bobines;
    sum.is_Helmholtz = field1.is_Helmholtz; // may be also field1.is_Helmholtz&&field2.is_Helmholtz
    sum.potential = field1.potential + field2.potential;

    return sum;
}


//----------------------------------
// Surdéfinition des entrées sorties
//----------------------------------

ostream& operator << (ostream &flux,  Field field)
{
    field.write(flux);
    return(flux);
}

istream& operator >> (istream &flux,  Field & field) //at est modifié!
{
    field.read(flux);
    return(flux);
}


// Get field at a given position from the grid
Vecteur3D  Field::get_Field_interpolation2D(const Vecteur3D pos, int type_field_read) const
{
    const double x= pos.x(), y = pos.y(), z=pos.z();
    const double r = sqrt(x*x + y*y); // cylindrical coordinates;
    double Fz,Fr;
    Fr = Interpolator2D.GetValue(matrix_field, r, z, 2, type_field_read ); //Return evaluation of Br(r,z). Br is on column 2 of the matrix
    Fz = Interpolator2D.GetValue(matrix_field, r, z, 3, type_field_read ); //Return evaluation of Bz(r,z). Bz is on column 3 of the matrix


    Vecteur3D field;

    if (r<=SMALL_NUMBER) // r==0 is too restrictive and we are safer like that
    {
        field = Vecteur3D(0,0,Fz);
    }
    else
    {
        field = Vecteur3D(Fr*x/r, Fr*y/r, Fz); // Fx=Fr Cos(theta); Fy= Fr Sin(theta)
    }
    return field;
}


// Get Grad (F.F) = gradient of the norm of the field squared at a given position from the grid
// Also give the field F
void  Field::get_norm2_grad_Field_interpolation2D(const Vecteur3D pos, Vecteur3D  &field, Vecteur3D &grad_F2, int type_field_read) const
{
    const double x= pos.x(), y = pos.y(), z=pos.z();
    const double r = sqrt(x*x + y*y); // cylindrical coordinates;

    double Fr,dFr_dr,dFr_dz;
    double Fz,dFz_dr,dFz_dz;


    Interpolator2D.GetValue(matrix_field, r, z, Fr, dFr_dr, dFr_dz, 2, type_field_read); // Fr is on column 2 of the matrix
    Interpolator2D.GetValue(matrix_field, r, z, Fz, dFz_dr, dFz_dz, 3, type_field_read);

    double grad_F2_z = 2.*(Fr*dFr_dz + Fz*dFz_dz); // Grad (F.F) = 2 (F.grad) F

    if (r<=SMALL_NUMBER)
    {
        field = Vecteur3D(0,0,Fz);
        grad_F2 = Vecteur3D(0,0,grad_F2_z);
    }
    else
    {
        field = Vecteur3D(Fr*x/r, Fr*y/r, Fz); // Fx=Fr Cos(theta); Fy= Fr Sin(theta)
        double grad_F2_r = 2.*(Fr*dFr_dr + Fz*dFz_dr);
        grad_F2 = Vecteur3D(grad_F2_r*x/r, grad_F2_r*y/r, grad_F2_z);
    }
    return;
}

// Get field at a given position
Vecteur3D  Field::get_Field(const Vecteur3D &pos) const
{
    Vecteur3D Field_value;

    if (this->type_field_read == 0 )
    {
        Field_value = Vecteur3D(F0.x()+F1.x()*pos.x() + F2.x()*pos.x()*pos.x() + Fn.x()*pow(pos.x(),n), F0.y()+F1.y()*pos.y() +  F2.y()*pos.y()*pos.y() + Fn.y()*pow(pos.y(),n), F0.z()+F1.z()*pos.z()+  F2.z()*pos.z()*pos.z() + Fn.z()*pow(pos.z(),n));

        return Field_value;
    }

    if (this->type_field_read == 1) // N-coils of Ox axis

    {
        Field_value =  get_field_coils(Vecteur3D(pos.y(),pos.z(),pos.x()), this->rayon_bobines, this->courant_bobines, this->Nb_bobines, this->gap_bobines, this->is_Helmholtz);
        // les bobines sont calculées comme étant d'axes Oz alors qu'elle seront d'axe Ox en réalité. Il faudra donc appeler ces fonctions non pas en (x,y,z) mais en (y,z,x)
        return Vecteur3D(Field_value.z(),Field_value.x(),Field_value.y());
    }


    if (this->type_field_read == 2 || this->type_field_read == 3) // 2: Field map 3D but from 2D cylindrical symmetry: r,z
    {
        Field_value = get_Field_interpolation2D(pos,type_field_read);
    }

    if (this->type_field_read >= 4)
    {
        cout << "need to be done " << endl;
    }


    return Field_value;

}


// Get potential at a given position
// Be careful give the potential as if it was an electric potential (so it integrates the field)
double  Field::get_electric_potential(const Vecteur3D pos)  const
{
    double pot = 0.;

    if (this->type_field_read == 0)
    {
        double x= pos.x(), y= pos.y(), z = pos.z();
        pot +=  - F0.x()* x - F1.x()* x*x/2. -  F2.x()* x*x*x/3. - Fn.x() * pow(x,n+1)/(n+1.);  // Ex = - d/dx  pot
        pot +=  - F0.y()* y - F1.y()* y*y/2. -  F2.y()* y*y*y/3. - Fn.y() * pow(y,n+1)/(n+1.);
        pot +=  - F0.z()* z - F1.z()* z*z/2. -  F2.z()* z*z*z/3. - Fn.z() * pow(z,n+1)/(n+1.);
        return pot;
    }
    else
    {
        cout << " electric potential not implemented for this type of field " << endl;
    }

    return pot;

}


// Calcul du gradient de la norme (au carré) du champ local F
// Grad (F.F) = 2 (F.grad) F
Vecteur3D Field::get_grad_field_F2(const Vecteur3D& pos) const
{
    Vecteur3D grad_F2;

    if (this->type_field_read == 0)
    {
        double along_x = 2.*( F0.x()+F1.x()*pos.x() + F2.x()*pos.x()*pos.x() +  Fn.x()*pow(pos.x(),n) )*( F1.x()+2.*F2.x()*pos.x() +  n*Fn.x()*pow(pos.x(),n-1));
        double along_y = 2.*( F0.y()+F1.y()*pos.y() + F2.y()*pos.y()*pos.y() +  Fn.y()*pow(pos.y(),n) )*( F1.y()+2.*F2.y()*pos.y() +  n*Fn.y()*pow(pos.y(),n-1));
        double along_z = 2.*( F0.z()+F1.z()*pos.z() + F2.z()*pos.z()*pos.z() +  Fn.z()*pow(pos.z(),n) )*( F1.z()+2.*F2.z()*pos.z() +  n*Fn.z()*pow(pos.z(),n-1));
        grad_F2 = Vecteur3D(along_x, along_y, along_z);
        return grad_F2;
    }

    if (this->type_field_read == 1)
    {
        grad_F2 = get_grad_normfield2_coils(Vecteur3D(pos.y(),pos.z(),pos.x()), this->rayon_bobines, this->courant_bobines, this->Nb_bobines, this->gap_bobines, this->is_Helmholtz);
        // les bobines sont calculées comme étant d'axes Oz alors qu'elle seront d'axe Ox en réalité. Il faudra donc appeler ces fonctions non pas en (x,y,z) mais en (y,z,x)
        return Vecteur3D(grad_F2.z(),grad_F2.x(),grad_F2.y());
    }

    if (this->type_field_read == 2 || this->type_field_read == 3)
    {
        Vecteur3D  field;
        get_norm2_grad_Field_interpolation2D(pos, field, grad_F2,type_field_read);
    }


    if (this->type_field_read >= 4)
    {
        cout << "need to be done " << endl;
    }

    return grad_F2;
}



/***************  BOBINE  *****************/

// CALCUL APPROCHE. Les formules analytiques existes mais sont complexes.
// Ici ce sont des formules plus simples mais très proche de la réalité (cf bz2, brho2 Mathematica ficher calcul_Sisyphus.nb (ou plutot magnetiqe.nb)
// Champ   pour bobine rayon r centrée en z=0 axe Oz avec courant de I_coil ampère

// Pour info le champ au centre est µ0 I/ (2 r)
// Champ crée par une bobine positionée en pos_bobine. Donne aussi les valeurs  Brho_sur_rho et Bz
Vecteur3D Field::get_field_coil(const Vecteur3D pos, double r, double I_coil, const Vecteur3D pos_bobine,   double &Brho_sur_rho, double &Bz) const
{

    double x=pos.x() - pos_bobine.x();
    double y=pos.y() - pos_bobine.y();
    double z=pos.z() - pos_bobine.z();


    double rho = sqrt(x*x+y*y);
    double intermediate = (r + rho)*(r + rho) + z*z;
    double denominator =  2.*(z*z + (r - rho)*(r - rho) )*intermediate*sqrt(intermediate); // 2 (z^2 + (r -rho)^2) (z^2 + (r +rho)^2)^(3/2)
    Brho_sur_rho = I_coil * MU0*(3./2.)*r*r*z/denominator; // j'ai changé 3 en 3/2
 // cerr << " attention error sur gradient " << endl;
    double Bx =  x*Brho_sur_rho;// Bx = brho Cos theta avec Cos theta = x/rho *)
    double By =  y*Brho_sur_rho; // By = brho Sin theta avec Sin theta = y/rho *)
    Bz = I_coil * MU0*r*r*(r*r + r*rho - 2*rho*rho + z*z)/denominator;

    return  Vecteur3D(Bx,By,Bz);
}


Vecteur3D Field::get_field_coil(const Vecteur3D pos, double r, double I_coil, const Vecteur3D pos_bobine) const
{

    double Brho_sur_rho=0.;
    double Bz=0.;

    return  get_field_coil(pos, r, I_coil, pos_bobine, Brho_sur_rho, Bz);
}



// Champ pour Nb_coils bobines de rayon r la première centrée en z=0 les autres séparées de  z_gap
Vecteur3D Field::get_field_coils(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz, double &Brho_sur_rho_tot, double &Bz_tot) const
{
    Vecteur3D B(0.,0.,0.);

    double x=pos.x();
    double y=pos.y();
    double z=pos.z();

    double Brho_sur_rho=0.;
    double Bz=0.;
    Brho_sur_rho_tot=0.;
    Bz_tot=0.;

    for (int N = 0; N < Nb_coils; N++)
    {
        double z_bobine =  (N-0.5) * z_gap;
        B += get_field_coil(Vecteur3D(x, y,  z),  r, I_coil, Vecteur3D(0., 0., z_bobine), Brho_sur_rho, Bz);
        Brho_sur_rho_tot += Brho_sur_rho; // Rho ne change pas lorsqu'on scan z!
        Bz_tot += Bz;
    }

    if (is_Helmholtz == 1)
        for (int N = 0; N < Nb_coils; N++)
        {
            double z_bobine =  r + (N-0.5) * z_gap; // On ajoute ue bobine décalée du rayon (configuration Helmholtz)
            B += get_field_coil(Vecteur3D(x, y,  z),  r, I_coil, Vecteur3D(0., 0., z_bobine), Brho_sur_rho, Bz);
            Brho_sur_rho_tot += Brho_sur_rho; // Rho ne change pas lorsqu'on scan z!
            Bz_tot += Bz;
        }

    if (is_Helmholtz == -1) //anti-Helmholtz
        for (int N = 0; N < Nb_coils; N++)
        {
            double z_bobine =  r + (N-0.5) * z_gap; // On ajoute ue bobine décalée du rayon (configuration Helmholtz)
            B += get_field_coil(Vecteur3D(x, y,  z),  r, -I_coil, Vecteur3D(0., 0., z_bobine), Brho_sur_rho, Bz);
            Brho_sur_rho_tot += Brho_sur_rho; // Rho ne change pas lorsqu'on scan z!
            Bz_tot += Bz;
        }
    return B;
}


Vecteur3D Field::get_field_coils(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz) const
{
    double Brho_sur_rho_tot=0.;
    double Bz_tot=0.;
    return  Field::get_field_coils(pos,  r, I_coil, Nb_coils, z_gap, is_Helmholtz, Brho_sur_rho_tot, Bz_tot);
}

// Gradient du champ   pour une bobines rayon r centrée en z=0 axe Oz
// En fait donne   dBz_drho, dBz_dz,  dBrho_drho, dBrho_dz
void Field::get_grad_field_coil_4(const Vecteur3D pos, double r, double I_coil, double &dBz_drho_sur_rho, double &dBz_dz,  double &dBrho_drho, double &dBrho_dz_sur_rho) const
{
    double x=pos.x();
    double y=pos.y();
    double z=pos.z();

    double rho = sqrt(x*x+y*y);

    double d1= z*z + (r - rho)*(r - rho);
    double d2 = z*z + (r + rho)*(r + rho) ;
    double denominator =  2.*d1*d1 * d2*d2*sqrt(d2);

    dBz_drho_sur_rho =  3.*r*r * (-3.*z*z* (r*r + z*z) +   2.* r*r*r*rho - (2.* r*r + z*z) * rho*rho - 2.* r*rho*rho*rho + 2.*rho*rho*rho*rho);
    dBz_dz = - 3. *r*r* z * ((r*r + z*z)*(r*r + z*z) +   r * (r*r + z*z)*rho - (r*r + 3.*z*z)*rho*rho + 3.* r*rho*rho*rho -   4.*rho*rho*rho*rho);
    dBrho_drho =  3.* r*r* z * ((r*r + z*z)*(r*r + z*z) - r * (r*r + z*z)* rho +  3. * (r - z) * (r + z) *rho*rho + r*rho*rho*rho - 4.*rho*rho*rho*rho);
    dBrho_dz_sur_rho = 3. * r*r * (r*r*r*r - 4.* z*z*z*z + 2. * r * z*z*rho - 3. * z*z*rho*rho +rho*rho*rho*rho - r*r* (3.* z*z + 2.*rho*rho));

    dBz_drho_sur_rho =  dBz_drho_sur_rho * I_coil * MU0/denominator;
    dBz_dz = dBz_dz*I_coil * MU0/denominator;
    dBrho_drho = dBrho_drho *I_coil * MU0/denominator;
    dBrho_dz_sur_rho = dBrho_dz_sur_rho * I_coil * MU0/denominator;

    return;
}

// Gradient du champ   pour des bobines rayon r centrée en z=0 axe Oz
// En fait donne   dBz_drho, dBz_dz,  dBrho_drho, dBrho_dz
void Field::get_grad_field_coils_4(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz, double &dBz_drho_sur_rho_tot, double &dBz_dz_tot,  double &dBrho_drho_tot, double &dBrho_dz_sur_rho_tot) const
{
    double x=pos.x();
    double y=pos.y();
    double z=pos.z();

    double dBz_drho_sur_rho=0.;
    double dBz_dz=0.;
    double dBrho_drho=0.;
    double dBrho_dz_sur_rho=0.;

    dBz_drho_sur_rho_tot=0.;
    dBz_dz_tot=0.;
    dBrho_drho_tot=0.;
    dBrho_dz_sur_rho_tot=0.;

    for (int N = 0; N < Nb_coils; N++)
    {
        double z_bobine =  (N-0.5) * z_gap;
        Field::get_grad_field_coil_4(Vecteur3D(x,y,z-z_bobine), r, I_coil, dBz_drho_sur_rho, dBz_dz,  dBrho_drho, dBrho_dz_sur_rho);

        dBz_drho_sur_rho_tot += dBz_drho_sur_rho;
        dBz_dz_tot += dBz_dz;
        dBrho_drho_tot += dBrho_drho;
        dBrho_dz_sur_rho_tot += dBrho_dz_sur_rho;
    }

    if (is_Helmholtz == 1)
        for (int N = 0; N < Nb_coils; N++)
        {
            double z_bobine =  r+ (N-0.5) * z_gap; // On ajoute ue bobine décalée du rayon (configuration Helmholtz)
            Field::get_grad_field_coil_4(Vecteur3D(x,y,z-z_bobine), r, I_coil, dBz_drho_sur_rho, dBz_dz,  dBrho_drho, dBrho_dz_sur_rho);

            dBz_drho_sur_rho_tot += dBz_drho_sur_rho;
            dBz_dz_tot += dBz_dz;
            dBrho_drho_tot += dBrho_drho;
            dBrho_dz_sur_rho_tot += dBrho_dz_sur_rho;
        }


    if (is_Helmholtz == -1) //anti-Helmholtz
        for (int N = 0; N < Nb_coils; N++)
        {
            double z_bobine =  r+ (N-0.5) * z_gap; // On ajoute ue bobine décalée du rayon (configuration anti-Helmholtz)
            Field::get_grad_field_coil_4(Vecteur3D(x,y,z-z_bobine), r, -I_coil, dBz_drho_sur_rho, dBz_dz,  dBrho_drho, dBrho_dz_sur_rho);

            dBz_drho_sur_rho_tot += dBz_drho_sur_rho;
            dBz_dz_tot += dBz_dz;
            dBrho_drho_tot += dBrho_drho;
            dBrho_dz_sur_rho_tot += dBrho_dz_sur_rho;
        }

    return;
}


// Gradient de la Norme  du champ   pour des bobines rayon r centrée en z=0 axe Oz
// Grad B.B = 2 (B.Grad) B
Vecteur3D Field::get_grad_normfield2_coils(const Vecteur3D pos, double r, double I_coil, int Nb_coils, double z_gap, int is_Helmholtz) const
{
    Vecteur3D B(0.,0.,0.);
    Vecteur3D gradB2(0.,0.,0.);

    double x=pos.x();
    double y=pos.y();
    double rho = sqrt(x*x+y*y);

    double Brho_sur_rho=0.;
    double Bz=0.;
    B = Field::get_field_coils(pos, r, I_coil, Nb_coils, z_gap, is_Helmholtz, Brho_sur_rho, Bz);

    double dBz_drho_sur_rho=0.;
    double dBz_dz=0.;
    double dBrho_drho=0.;
    double dBrho_dz_sur_rho=0.;
    Field::get_grad_field_coils_4(pos, r, I_coil, Nb_coils, z_gap, is_Helmholtz, dBz_drho_sur_rho, dBz_dz, dBrho_drho, dBrho_dz_sur_rho);


    double GradB2rho_sur_rho = 2. *(Brho_sur_rho*dBrho_drho + Bz*dBz_drho_sur_rho); // Mettre le sur rho permet d'éviter les x/rho lorsque rho=0
    double GradB2z = 2. *(Brho_sur_rho*rho*dBrho_dz_sur_rho*rho + Bz*dBz_dz);


    gradB2 +=  Vecteur3D(GradB2rho_sur_rho *x, GradB2rho_sur_rho*y, GradB2z);

    return gradB2;
}



/** READ THE MATRIX ****/

// Erase and Read the matrix of discrete value from a file and return the number of elements
// The  file depends on the type of field but should be ordered x By Bz.
// It should be first along the n1 axis [0,0,0], [1,0,0], ..., [n1-1,0,0], [0,1,0], ...
// and  equally spaced along each axis,
// THE FIRST LINE CONTAINS ANY COMMENTS you want!
// THE SECOND CONTAINS NAMES like x,y,z, Bx, By, Bz, dBx_dx ... SEPARATED BY TABS
// Return number of lines;
int Field::Read_Matrix(const char* nom_file)
{
    ifstream file(nom_file);

    cout << "name file " << nom_file << endl << endl;

    if ( !file || nom_file == NULL)
    {
        cerr << "No able to open the file " << nom_file << endl;  // Better to not put because sometimes their is no file (and so if this line is here, we will have all the time a message) and we just as the defautl values
        file.close();
        return 0;
    }

    cout << "Reading field map from ASCII file" << endl;

// READ MATRIX SIZE!
    int number_matrix_lines = number_line_file(nom_file) - 2; // To remove the first lines of comments
    int number_matrix_columns = 1;
    string line;
    getline(file, line);// First line
    getline(file, line); // second line


    for (int i=0; i< (int) line.size(); i++)
    {
        if (line.data()[i]== '\t')  // tab = '\t'
        {
            number_matrix_columns++;
        }
    }


    matrix_field.clear(); // To avoid to insert the default file at the begigning

    vector<double> column[100]; //  column
    // We assume that there is less than 100 columns x,y,z, BX By Bz dBx/dx etc ... !
    if (number_matrix_columns > 100)
    {
        cerr << "too many column in field map file " << endl;
        exit(1);
    }

    for (int i = 0; i < number_matrix_lines; ++i)
        for (int j = 0; j < number_matrix_columns; ++j)
        {
            double current_value;
            file >> current_value;
            column[j].push_back(current_value);
        }

    for (int j = 0; j < number_matrix_columns; ++j)
        matrix_field.push_back(column[j]); // matrix_field[j] = column number j of the file

    // SO BE CAREFUL matrix_field[j][i]  is the element [i,j] of the line i, column j
    file.close();

    return number_matrix_lines;
}

// Calculation and output of partial derivation of the field
// For  instance for the type_field 2 (that does not contains the derivatives. I calculate them)
void Field::Calculate_Derivative_Matrix(int number_matrix_lines, const char * output_filename) const
{

    ofstream file_out(output_filename);

    for (int i = 0; i < number_matrix_lines; ++i)
    {
        for (int j = 0; j < 4; ++j) // THe first 4 columns are x,y, Fx,Fy
        {
            file_out << matrix_field[j][i]  << " " ;
        } // Then we add the derivates
        file_out << Interpolator2D.XDerivative(matrix_field, i,2)<< " " ;
        file_out << Interpolator2D.XDerivative(matrix_field, i,3)<< " " ;

        file_out << Interpolator2D.YDerivative(matrix_field, i,2)<< " " ;
        file_out << Interpolator2D.YDerivative(matrix_field, i,3)<< " " ;

        file_out << Interpolator2D.XYDerivative(matrix_field, i,2)<< " " ;
        file_out << Interpolator2D.XYDerivative(matrix_field, i,3);

        file_out << endl;
    }
}


// Using the Matrix of the field (that should have been read before!)  we extract the grid parameters and initialize the Interpolators
// Field[m_xindex][m_yindex] is given by matrix_field[j][i] where j is the column corresponding to the field and i = m_index*m_xdimension + m_yindex
void Field::Init_Interpolator()
{
    int number_matrix_lines = matrix_field[0].size();
    if (this->type_field_read == 2  || this->type_field_read == 3) // Field map 3D from 2D cylindrical symmetry: r,z,
    {
        double r0 =  matrix_field[0][0]; // Recall matrix_field[j][i]  is the element [i,j] of the line i, column j;
        double z0 = matrix_field[1][0];
        double delta_z = matrix_field[1][1] - matrix_field[1][0]  ; // The file is ordered first in 1st column then in second!

        double rfinal = matrix_field[0][number_matrix_lines-1];
        double zfinal = matrix_field[1][number_matrix_lines-1];

        int zDimension = 1.001 + (zfinal-z0)/delta_z;  // Number of different  z values: z0, z0 +delta_z, ..., zfinal =  z0 + (Zdimension-1) * delta_z
        int rDimension =  0.001 + number_matrix_lines/zDimension; // Number of different  r values: r0, r0 +delta r, ..., rfinal = r0 + (rdimension-1) * delta_r
// To avoid round errors (seen!) we add 0.001
        double delta_r = (rfinal - r0)/(rDimension-1);

        cout << " init Interpolator 2D " << endl;
        Interpolator2D.Init(rDimension,  zDimension, delta_r, delta_z, r0, z0); //Interpolator2D.Init(XDimension,  YDimension, delta_x, delta_y, X0, Y0);
    }

    if (this->type_field_read >= 4)
    {

        cout << " init Interpolator : NOT DONE YET " << endl;
    }
}
