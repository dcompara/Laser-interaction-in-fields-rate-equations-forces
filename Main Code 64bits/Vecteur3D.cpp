// Classe Vecteur à 3 dimensions
// Basée sur ThreeVector.h,
// as part of the CLHEP - a Class Library for High Energy Physics.
// http://wwwinfo.cern.ch/asd/lhc++/clhep/index.html*

// Cette classe était assez incomplète par exemple l'operator / qui manque ...
// Je l'ai en partie completée
// Il y a en fait beaucoup de fonction en ligne qui ne sont pas définie ensuite


// Voir aussi la classe « vecteur_ » (exemple du cours C++)
//         Patrick TRAU ULP-IPST Strasbourg novembre 04
// Pour les stream


// En gros
// Une fois défini un Vecteur3D point=Vecteur3D(0,1,2); ou
// Vecteur3D point(0,1,2);
// point.x() pour prendre la coordonnée X

#include "Vecteur3D.h"


extern const Vecteur3D HepXHat, HepYHat, HepZHat;

ostream& operator << (ostream &f,Vecteur3D v)
{
    v.affiche(f);
    return(f);
}

istream& operator >> (istream &f,Vecteur3D &v) //v est modifié!
{
    v.saisie(f);
    return(f);
}

typedef Vecteur3D HepThreeVectorD;
typedef Vecteur3D HepThreeVectorF;


// --------------
// Global methods
// --------------

Vecteur3D operator / (const Vecteur3D & p, double a)
{
    return Vecteur3D(p.x()/a, p.y()/a, p.z()/a);
}

Vecteur3D operator + (const Vecteur3D & a, const Vecteur3D & b)
{
    return Vecteur3D(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
}

Vecteur3D operator - (const Vecteur3D & a, const Vecteur3D & b)
{
    return Vecteur3D(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
}

Vecteur3D operator * (const Vecteur3D & p, double a)
{
    return Vecteur3D(a*p.x(), a*p.y(), a*p.z());
}

Vecteur3D operator * (double a, const Vecteur3D & p)
{
    return Vecteur3D(a*p.x(), a*p.y(), a*p.z());
}

double operator * (const Vecteur3D & a, const Vecteur3D & b)
{
    return a.dot(b);
}



//Add by PAULINE 26/03/2015

Vecteur3D operator / (const Vecteur3D & a, const Vecteur3D & b)
{
    return Vecteur3D(a.x()/b.x(), a.y()/b.y(), a.z()/b.z());
}

Vecteur3D operator - ( const Vecteur3D & p, double a)
{
    return Vecteur3D(p.x() - a, p.y()-a, p.z()-a);
}

Vecteur3D operator + ( const Vecteur3D & p, double a)
{
    return Vecteur3D(p.x() + a, p.y()+ a, p.z()+ a);
}

 Vecteur3D racine(const Vecteur3D & a) // function square root of a 3D vector
{
    return Vecteur3D(sqrt(a.x()),sqrt(a.y()),sqrt(a.z()));
}
 Vecteur3D Hadamard(const Vecteur3D & a, const Vecteur3D & b) // special multiplication: component by component. Multiplication d'Hadamard
{
    return Vecteur3D(a.x()*b.x(), a.y()*b.y(),a.z()*b.z());
}


 Vecteur3D abso(const Vecteur3D & a) //absolute value
{
    return Vecteur3D(abs(a.x()), abs(a.y()),abs(a.z()));
}


// --------------------------
// Set in various coordinates
// --------------------------



void Vecteur3D::gaussian_initialisation (const gsl_rng * r, const Vecteur3D sigma)
{
    dx = gsl_ran_gaussian (r,sigma.getX());
    dy = gsl_ran_gaussian (r,sigma.getY());
    dz = gsl_ran_gaussian (r,sigma.getZ());
}


// Laplace initilisation utilise
// This function returns a random variate from the Laplace distribution with width a.
// p(x)dx = (1/2a)exp(−|x/a|)dx.
void Vecteur3D::laplace_initialisation (const gsl_rng * r, const Vecteur3D sigma)
{
    dx = gsl_ran_laplace (r,sigma.getX());
    dy = gsl_ran_laplace (r,sigma.getY());
    dz = gsl_ran_laplace (r,sigma.getZ());
}

// x in Laplace, y in Laplace and z in gaussian
void Vecteur3D::laplace_laplace_gauss_initialisation (const gsl_rng * r, const Vecteur3D sigma)
{
    dx = gsl_ran_laplace (r,sigma.getX());
    dy = gsl_ran_laplace (r,sigma.getY());
    dz = gsl_ran_gaussian  (r,sigma.getZ());
}

