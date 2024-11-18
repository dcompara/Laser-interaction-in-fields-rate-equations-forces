
// Class: 3D Vector
// Based on ThreeVector.h,
// as part of the CLHEP - a Class Library for High Energy Physics.
// http://wwwinfo.cern.ch/asd/lhc++/clhep/index.html*

// This class was originally quite incomplete, for example, it lacked the operator `/`.
// I have partially completed it.
// There are, in fact, many inline functions that are not defined elsewhere.

// See also the class "vecteur_" (example from the C++ course)
//         Patrick TRAU, ULP-IPST Strasbourg, November 2004
// For stream manipulations.

// In summary:
// Once a `Vecteur3D` is defined, for example:
// `Vecteur3D point = Vecteur3D(0, 1, 2);` or
// `Vecteur3D point(0, 1, 2);`
// You can access the X coordinate with `point.x()`.




#ifndef VECTOR3D_SEEN
#define VECTOR3D_SEEN

#include <cmath>
#include <iostream>
#include <fstream>                         // to read the data and put them in files
#include <gsl/gsl_rng.h>                    // To have random_number generator
#include <gsl/gsl_randist.h>

#include <iostream>
using namespace std;

class Vecteur3D
{

public:

    // Basic properties and operations on 3-vectors:

    enum
    {
        X=0, Y=1, Z=2, NUM_COORDINATES=3, SIZE=NUM_COORDINATES
    };
    // Safe indexing of the coordinates when using with matrices, arrays, etc.
    // (BaBar)

    inline Vecteur3D(double x = 0.0, double y = 0.0, double z = 0.0);
    // The constructor.

    inline Vecteur3D(const Vecteur3D &);
    // The copy constructor.

    inline ~Vecteur3D();
    // The destructor.  Not virtual - inheritance from this class is dangerous.

    void affiche(ostream & pr)
    {
        // pr<<"["<<dx<<","<<dy<<","<<dz<<"]";
        pr<<" "<<dx<<" "<<dy<<" "<<dz<<" ";
        return;
    }

    void saisie(istream & f)
    {
        if (&f == &cin)
            cout<<"entrez x : ";
        f>>dx;
        if (&f==&cin)
            cout<<"entrez y : ";
        f>>dy;
        if (&f==&cin)
            cout<<"entrez z : ";
        f>>dz;
        return;
    }

    double operator () (int i) const
// Get components by index -- 0-based (Geant4)
    {
        if (i == 0)
        {
            return dx;
        }
        else if (i == 1)
        {
            return dy;
        }
        else if (i == 2)
        {
            return dz;
        }
        else
        {
            // cerr << "Vector3D::operator(): bad index" << endl;
            return 0.0;
        }
    }

    inline double operator [] (int) const;
    // Get components by index -- 0-based (Geant4)

    double & operator () (int i)
    // Get components by index.  0-based.
    {
        if (i == 0)
        {
            return dx;
        }
        else if (i == 1)
        {
            return dy;
        }
        else if (i == 2)
        {
            return dz;
        }
        cerr << "use double & operator () (int i)  with inot 0,1,2 " << endl;
        return dx; // To avoid warning in the compiler that "reaches end of non-void function [-Wreturn-type]"

    }




    inline double & operator [] (int);
    // Get components by index.  0-based.

    inline double x() const;
    inline double y() const;
    inline double z() const;
    // The components in cartesian coordinate system.  Same as getX() etc.

    // Get by each coordinate int =0 for x, 1 for y and 2 for z
    inline double get(int);


    inline void setX(double);
    inline void setY(double);
    inline void setZ(double);
    // Set the components in cartesian coordinate system.

    // set by each coordinate int =0 for x, 1 for y and 2 for z
    inline void set(int,double);

    inline void set( double x, double y, double z);
    // Set all three components in cartesian coordinate system.

    inline double phi() const;
    // The azimuth angle.

    inline double theta() const;
    // The polar angle.

    inline double cosTheta() const;
    // Cosine of the polar angle.

    inline double cos2Theta() const;
    // Cosine squared of the polar angle - faster than cosTheta(). (ZOOM)

    inline double mag2() const;
    // The magnitude squared (r^2 in spherical coordinate system).

    inline Vecteur3D square() const;
    // The squared (x^2,y^2,z^2) of each component .

    inline double mag() const;
    // The magnitude (r in spherical coordinate system).

    inline void setPhi(double);
    // Set phi keeping mag and theta constant (BaBar).

    inline void setTheta(double);
    // Set theta keeping mag and phi constant (BaBar).

    void setMag(double);
    // Set magnitude keeping theta and phi constant (BaBar).

    inline double perp2() const;
    // The transverse component squared (rho^2 in cylindrical coordinate system).

    inline double perp() const;
    // The transverse component (rho in cylindrical coordinate system).

    inline void setPerp(double);
    // Set the transverse component keeping phi and z constant.

    void setCylTheta(double);
    // Set theta while keeping transvers component and phi fixed

    inline double perp2(const Vecteur3D &) const;
    // The transverse component w.r.t. given axis squared.

    inline double perp(const Vecteur3D &) const;
    // The transverse component w.r.t. given axis.

    inline Vecteur3D & operator = (const Vecteur3D &);
    // Assignment.

    inline bool operator == (const Vecteur3D &) const;
    inline bool operator != (const Vecteur3D &) const;
    // Comparisons (Geant4).

    bool isNear (const Vecteur3D &, double epsilon=tolerance) const;
    // Check for equality within RELATIVE tolerance (default 2.2E-14). (ZOOM)
    // |v1 - v2|**2 <= epsilon**2 * |v1.dot(v2)|

    double howNear(const Vecteur3D & v ) const;
    // sqrt ( |v1-v2|**2 / v1.dot(v2) ) with a maximum of 1.
    // If v1.dot(v2) is negative, will return 1.

    double deltaR(const Vecteur3D & v) const;
    // sqrt( pseudorapity_difference**2 + deltaPhi **2 )

    inline Vecteur3D & operator += (const Vecteur3D &);
    // Addition.

    inline Vecteur3D & operator -= (const Vecteur3D &);
    // Subtraction.

    inline Vecteur3D operator - () const;
    // Unary minus.

    inline Vecteur3D & operator *= (double);
    // Scaling with real numbers.

    inline Vecteur3D & operator /= (double);
    // Division by (non-zero) real number.

    inline Vecteur3D unit() const;
    // Vector parallel to this, but of length 1.

    inline Vecteur3D orthogonal() const;
    // Vector orthogonal to this (Geant4).

    inline double dot(const Vecteur3D &) const;
    // double product.

    inline Vecteur3D cross(const Vecteur3D &) const;
    // Cross product.

    double angle(const Vecteur3D &) const;
    // The angle w.r.t. another 3-vector.

    double pseudoRapidity() const;
    // Returns the pseudo-rapidity, i.e. -ln(tan(theta/2))

    void setEta  ( double p );
    // Set pseudo-rapidity, keeping magnitude and phi fixed.  (ZOOM)

    void setCylEta  ( double p );
    // Set pseudo-rapidity, keeping transverse component and phi fixed.  (ZOOM)


    // = = = = = = = = = = = = = = = = = = = = = = = =
    //
    // Esoteric properties and operations on 3-vectors:
    //
    // 1 - Set vectors in various coordinate systems
    // 2 - Synonyms for accessing coordinates and properties
    // 3 - Comparisions (dictionary, near-ness, and geometric)
    // 4 - Intrinsic properties
    // 5 - Properties releative to z axis and arbitrary directions
    // 6 - Polar and azimuthal angle decomposition and deltaPhi
    //
    // = = = = = = = = = = = = = = = = = = = = = = = =

    // 1 - Set vectors in various coordinate systems

    inline void setRThetaPhi  (double r, double theta, double phi);
    // Set in spherical coordinates:  Angles are measured in RADIANS

    inline void setREtaPhi  ( double r, double eta,  double phi );
    // Set in spherical coordinates, but specify peudorapidiy to determine theta.

    inline void setRhoPhiZ   (double rho, double phi, double z);
    // Set in cylindrical coordinates:  Phi angle is measured in RADIANS

    void setRhoPhiTheta ( double rho, double phi, double theta);
    // Set in cylindrical coordinates, but specify theta to determine z.

    void setRhoPhiEta ( double rho, double phi, double eta);
    // Set in cylindrical coordinates, but specify pseudorapidity to determine z.

    // 2 - Synonyms for accessing coordinates and properties

    inline double getX() const;
    inline double getY() const;
    inline double getZ() const;
    // x(), y(), and z()

    inline double getR    () const;
    inline double getTheta() const;
    inline double getPhi  () const;
    // mag(), theta(), and phi()

    inline double r       () const;
    // mag()

    inline double rho     () const;
    inline double getRho  () const;
    // perp()

    double eta     () const;
    double getEta  () const;
    // pseudoRapidity()

    inline void setR ( double s );
    // setMag()

    inline void setRho ( double s );
    // setPerp()

    // 3 - Comparisions (dictionary, near-ness, and geometric)

    int compare (const Vecteur3D & v) const;
    bool operator > (const Vecteur3D & v) const;
    bool operator < (const Vecteur3D & v) const;
    bool operator>= (const Vecteur3D & v) const;
    bool operator<= (const Vecteur3D & v) const;
    // dictionary ordering according to z, then y, then x component

    inline double diff2 (const Vecteur3D & v) const;
    // |v1-v2|**2

    static double setTolerance (double tol);
    static inline double getTolerance ();
    // Set the tolerance used in isNear() for Vecteur3Ds

    bool isParallel (const Vecteur3D & v, double epsilon=tolerance) const;
    // Are the vectors parallel, within the given tolerance?

    bool isOrthogonal (const Vecteur3D & v, double epsilon=tolerance) const;
    // Are the vectors orthogonal, within the given tolerance?

    double howParallel   (const Vecteur3D & v) const;
    // | v1.cross(v2) / v1.dot(v2) |, to a maximum of 1.

    double howOrthogonal (const Vecteur3D & v) const;
    // | v1.dot(v2) / v1.cross(v2) |, to a maximum of 1.

    enum
    {
        ToleranceTicks = 100
    };

    // 4 - Intrinsic properties

    double beta    () const;
    // relativistic beta (considering v as a velocity vector with c=1)
    // Same as mag() but will object if >= 1

    double gamma() const;
    // relativistic gamma (considering v as a velocity vector with c=1)

    double coLinearRapidity() const;
    // inverse tanh (beta)

    // 5 - Properties relative to Z axis and to an arbitrary direction

    // Note that the non-esoteric CLHEP provides
    // theta(), cosTheta(), cos2Theta, and angle(const Vecteur3D&)

    inline double angle() const;
    // angle against the Z axis -- synonym for theta()

    inline double theta(const Vecteur3D & v2) const;
    // synonym for angle(v2)

    inline double cosTheta (const Vecteur3D & v2) const;
    inline double cos2Theta(const Vecteur3D & v2) const;
    // cos and cos^2 of the angle between two vectors

    inline Vecteur3D project () const;
    Vecteur3D project (const Vecteur3D & v2) const;
    // projection of a vector along a direction.

    inline Vecteur3D perpPart() const;
    inline Vecteur3D perpPart (const Vecteur3D & v2) const;
    // vector minus its projection along a direction.

    double rapidity () const;
    // inverse tanh(v.z())

    double rapidity (const Vecteur3D & v2) const;
    // rapidity with respect to specified direction:
    // inverse tanh (v.dot(u)) where u is a unit in the direction of v2

    double eta(const Vecteur3D & v2) const;
    // - ln tan of the angle beween the vector and the ref direction.

    // 6 - Polar and azimuthal angle decomposition and deltaPhi

    // Decomposition of an angle within reference defined by a direction:

    double polarAngle (const Vecteur3D & v2) const;
    // The reference direction is Z: the polarAngle is abs(v.theta()-v2.theta()).

    double deltaPhi (const Vecteur3D & v2) const;
    // v.phi()-v2.phi(), brought into the range (-PI,PI]

    double azimAngle  (const Vecteur3D & v2) const;
    // The reference direction is Z: the azimAngle is the same as deltaPhi

    double polarAngle (const Vecteur3D & v2,
                       const Vecteur3D & ref) const;
    // For arbitrary reference direction,
    //  polarAngle is abs(v.angle(ref) - v2.angle(ref)).

    double azimAngle  (const Vecteur3D & v2,
                       const Vecteur3D & ref) const;
    // To compute azimangle, project v and v2 into the plane normal to
    // the reference direction.  Then in that plane take the angle going
    // clockwise around the direction from projection of v to that of v2.


// Random  initiliastion

    void gaussian_initialisation (const gsl_rng * r, const Vecteur3D sigma);

    void laplace_initialisation (const gsl_rng * r, const Vecteur3D sigma);

    void laplace_laplace_gauss_initialisation (const gsl_rng * r, const Vecteur3D sigma);


protected:
    void setSpherical (double r, double theta, double phi);
    void setCylindrical (double r, double phi, double z);
    double negativeInfinity() const;

protected:

    double dx;
    double dy;
    double dz;
    // The components.

    static double tolerance;
    // default tolerance criterion for isNear() to return true.
}
;  // Vecteur3D


extern const Vecteur3D HepXHat, HepYHat, HepZHat;

ostream& operator << (ostream &f,Vecteur3D v);

istream& operator >> (istream &f,Vecteur3D &v); //v est modifié!

typedef Vecteur3D HepThreeVectorD;
typedef Vecteur3D HepThreeVectorF;

Vecteur3D operator / (const Vecteur3D &, double a);
// Division of 3-vectors by non-zero real number

Vecteur3D operator + (const Vecteur3D &, const Vecteur3D &);
// Addition of 3-vectors.

Vecteur3D operator - (const Vecteur3D &, const Vecteur3D &);
// Subtraction of 3-vectors.

double operator * (const Vecteur3D &, const Vecteur3D &);
// double product of 3-vectors.

Vecteur3D operator * (const Vecteur3D &, double a);
Vecteur3D operator * (double a, const Vecteur3D &);
// Scaling of 3-vectors with a real number


//ADD by PAULINE 26/03/2015
Vecteur3D operator / (const Vecteur3D & a, const Vecteur3D & b);
Vecteur3D operator - ( const Vecteur3D & p, double a);
Vecteur3D operator + ( const Vecteur3D & p, double a);

Vecteur3D racine(const Vecteur3D & a); // sqrt of a 3D vector
Vecteur3D Hadamard(const Vecteur3D & a, const Vecteur3D & b);
//modi pauline 8/04
Vecteur3D abso(const Vecteur3D & a); //absolute value
//Vecteur3D operator racine(const Vecteur3D & a);
// Scaling of 3-vectors with a real number



// ------------------
// Access to elements
// ------------------

// x, y, z

inline double & Vecteur3D::operator[] (int i)
{
    return operator()(i);
}

inline double   Vecteur3D::operator[] (int i) const
{
    return operator()(i);
}

inline double Vecteur3D::x() const
{
    return dx;
}
inline double Vecteur3D::y() const
{
    return dy;
}
inline double Vecteur3D::z() const
{
    return dz;
}

inline double Vecteur3D::getX() const
{
    return dx;
}
inline double Vecteur3D::getY() const
{
    return dy;
}
inline double Vecteur3D::getZ() const
{
    return dz;
}


// Get by each coordinate int =0 for x, 1 for y and 2 for zz
inline double Vecteur3D::get(int i)
{
    if (i == 0)
    {
        return dx;
    }
    else if (i == 1)
    {
        return dy;
    }
    else if (i == 2)
    {
        return dz;
    }
    return 0.;
}



inline void Vecteur3D::setX(double x)
{
    dx = x;
}
inline void Vecteur3D::setY(double y)
{
    dy = y;
}
inline void Vecteur3D::setZ(double z)
{
    dz = z;
}


// set by each coordinate int =0 for x, 1 for y and 2 for z
inline void Vecteur3D::set(int i,double value)
{
    if (i == 0)
    {
        dx = value;
    }
    else if (i == 1)
    {
        dy = value;
    }
    else if (i == 2)
    {
        dz = value;
    }
}


inline void Vecteur3D::set(double x, double y, double z)
{
    dx = x;
    dy = y;
    dz = z;
}

// --------------
// Global methods
// --------------

// --------------------------
// Set in various coordinates
// --------------------------

inline void Vecteur3D::setRThetaPhi
( double r, double theta, double phi )
{
    setSpherical(r, theta, phi);
}

inline void Vecteur3D::setREtaPhi
( double r, double eta,  double phi )
{
    setSpherical(r, 2*atan(exp(-eta)), phi);
}

inline void Vecteur3D::setRhoPhiZ
( double rho, double phi, double z)
{
    setCylindrical(rho, phi, z);
}

// ------------
// Constructors
// ------------

inline Vecteur3D::Vecteur3D(double x, double y, double z)
    : dx(x), dy(y), dz(z)
{}

inline Vecteur3D::Vecteur3D(const Vecteur3D & p)
    : dx(p.dx), dy(p.dy), dz(p.dz)
{}

inline Vecteur3D::~Vecteur3D()
{}

inline Vecteur3D & Vecteur3D::operator = (const Vecteur3D & p)
{
    dx = p.dx;
    dy = p.dy;
    dz = p.dz;
    return *this;
}

// ------------------
// Access to elements
// ------------------

// r, theta, phi

inline double Vecteur3D::mag2() const
{
    return dx*dx + dy*dy + dz*dz;
}
inline double Vecteur3D::mag()  const
{
    return sqrt(mag2());
}
inline double Vecteur3D::r()    const
{
    return mag();
}

inline double Vecteur3D::theta()  const
{
    return dx == 0.0 && dy == 0.0 && dz == 0.0 ? 0.0 : atan2(perp(),dz);
}
inline double Vecteur3D::phi() const
{
    return dx == 0.0 && dy == 0.0 ? 0.0 : atan2(dy,dx);
}

inline double Vecteur3D::getR()     const
{
    return mag();
}
inline double Vecteur3D::getTheta() const
{
    return theta();
}
inline double Vecteur3D::getPhi()   const
{
    return phi();
}
inline double Vecteur3D::angle()    const
{
    return theta();
}

inline double Vecteur3D::cosTheta() const
{
    double ptot = mag();
    return ptot == 0.0 ? 1.0 : dz/ptot;
}

inline double Vecteur3D::cos2Theta() const
{
    double ptot2 = mag2();
    return ptot2 == 0.0 ? 1.0 : dz*dz/ptot2;
}

inline void Vecteur3D::setR(double r)
{
    setMag(r);
}

inline void Vecteur3D::setTheta(double th)
{
    double ma   = mag();
    double ph   = phi();
    setX(ma*sin(th)*cos(ph));
    setY(ma*sin(th)*sin(ph));
    setZ(ma*cos(th));
}

inline void Vecteur3D::setPhi(double ph)
{
    double xy   = perp();
    setX(xy*cos(ph));
    setY(xy*sin(ph));
}

// perp, eta,

inline double Vecteur3D::perp2()  const
{
    return dx*dx + dy*dy;
}
inline double Vecteur3D::perp()   const
{
    return sqrt(perp2());
}
inline double Vecteur3D::rho()    const
{
    return perp();
}
inline double Vecteur3D::eta()    const
{
    return pseudoRapidity();
}

inline double Vecteur3D::getRho() const
{
    return perp();
}
inline double Vecteur3D::getEta() const
{
    return pseudoRapidity();
}

inline void Vecteur3D::setPerp(double r)
{
    double p = perp();
    if (p != 0.0)
    {
        dx *= r/p;
        dy *= r/p;
    }
}
inline void Vecteur3D::setRho(double rho)
{
    setPerp (rho);
}


// ----------
// Comparison
// ----------

inline bool Vecteur3D::operator == (const Vecteur3D& v) const
{
    return (v.x()==x() && v.y()==y() && v.z()==z()) ? true : false;
}

inline bool Vecteur3D::operator != (const Vecteur3D& v) const
{
    return (v.x()!=x() || v.y()!=y() || v.z()!=z()) ? true : false;
}

inline double Vecteur3D::getTolerance ()
{
    return tolerance;
}

// ----------
// Arithmetic
// ----------

inline Vecteur3D& Vecteur3D::operator += (const Vecteur3D & p)
{
    dx += p.x();
    dy += p.y();
    dz += p.z();
    return *this;
}

inline Vecteur3D& Vecteur3D::operator -= (const Vecteur3D & p)
{
    dx -= p.x();
    dy -= p.y();
    dz -= p.z();
    return *this;
}

inline Vecteur3D Vecteur3D::operator - () const
{
    return Vecteur3D(-dx, -dy, -dz);
}

inline Vecteur3D& Vecteur3D::operator *= (double a)
{
    dx *= a;
    dy *= a;
    dz *= a;
    return *this;
}

inline Vecteur3D & Vecteur3D::operator /= (double a)
{
    dx /= a;
    dy /= a;
    dz /= a;
    return *this;
}

// -------------------
// Combine two Vectors
// -------------------

inline Vecteur3D Vecteur3D::square() const
{
    return Vecteur3D(dx*dx, dy*dy, dz*dz); // The squared (x^2,y^2,z^2) of each component .
}


inline double Vecteur3D::diff2(const Vecteur3D & p) const
{
    return (*this-p).mag2();
}

inline double Vecteur3D::dot(const Vecteur3D & p) const
{
    return dx*p.x() + dy*p.y() + dz*p.z();
}

inline Vecteur3D Vecteur3D::cross(const Vecteur3D & p) const
{
    return Vecteur3D(dy*p.z()-p.y()*dz, dz*p.x()-p.z()*dx, dx*p.y()-p.x()*dy);
}

inline double Vecteur3D::perp2(const Vecteur3D & p)  const
{
    double tot = p.mag2();
    double ss  = dot(p);
    return tot > 0.0 ? mag2()-ss*ss/tot : mag2();
}

inline double Vecteur3D::perp(const Vecteur3D & p) const
{
    return sqrt(perp2(p));
}

inline Vecteur3D Vecteur3D::perpPart () const
{
    return Vecteur3D (dx, dy, 0);
}
inline Vecteur3D Vecteur3D::project () const
{
    return Vecteur3D (0, 0, dz);
}

inline Vecteur3D Vecteur3D::perpPart (const Vecteur3D & v2) const
{
    return ( *this - project(v2) );
}


// angle entre les 2 vecteurs // FAUX car c'est l'angle entre Oz et q!!
// A MON AVIS TOUTES LES FORMULES CI-DESSOUS SONT FAUSSES
inline double Vecteur3D::angle(const Vecteur3D & q) const
{
    return acos(cosTheta(q));
}
//
//inline double Vecteur3D::angle_2_vec(const Vecteur3D & v2) const
//{
//    return acos(v2.dot(*this)/(mag(v2)*(mag(*this))); // v1 v2 cos(theta) = v1.v2
//}


inline double Vecteur3D::azimAngle(const Vecteur3D & v2) const
{
    return deltaPhi(v2);
}

inline double Vecteur3D::cosTheta(const Vecteur3D & v2) const
{
    return((*this-v2).cosTheta());
}
inline double Vecteur3D::cos2Theta(const Vecteur3D & v2) const
{
    return((*this-v2).cos2Theta());
}
// cos and cos^2 of the angle between MM2 and axe Oz (axe Oz est theta=0).

inline double Vecteur3D::theta(const Vecteur3D & v2) const
{
    return((*this-v2).theta());
}





// ----------
// Properties
// ----------

inline Vecteur3D Vecteur3D::unit() const
{
    double  tot = mag2();
    Vecteur3D p(x(),y(),z());
    return tot > 0.0 ? p *= (1.0/sqrt(tot)) : p;
}

inline Vecteur3D Vecteur3D::orthogonal() const
{
    double x = dx < 0.0 ? -dx : dx;
    double y = dy < 0.0 ? -dy : dy;
    double z = dz < 0.0 ? -dz : dz;
    if (x < y)
    {
        return x < z ? Vecteur3D(0,dz,-dy) : Vecteur3D(dy,-dx,0);
    }
    else
    {
        return y < z ? Vecteur3D(-dz,0,dx) : Vecteur3D(dy,-dx,0);
    }
}



#endif
