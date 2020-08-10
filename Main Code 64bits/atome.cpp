#include "atome.h"
// ------------
// Constructors
// ------------


Atome::Atome()                       // Constructeur par défaut
{
    pos.set(0.,0.,0.);
    vel.set(0.,0.,0.);
    acc.set(0.,0.,0.);
    mass = 0.;
    charge = 0.;
    name ="";
};

//Atome::~Atome()
//{}
//;                    // Destructeur

Atome::Atome(const Atome & at)               // Constructeur de (re)copie
{
    pos = at.pos;
    vel = at.vel;
    acc = at.acc;
    mass = at.mass;
    charge = at.charge;
    name = at.name;
};

Atome & Atome::operator = (const Atome & at)          // Affectation par recopie
{
    if (this != &at) // On vérifie que les objets ne sont pas les mêmes !
    {
        pos = at.pos;
        vel = at.vel;
        acc = at.acc;
        mass = at.mass;
        charge = at.charge;
        name = at.name;
    }
    return *this;
}

//--------------------------------------------
// Surcharge des opérateurs +,-,*,/
// + (-) = addition (soustraction) membre à membre.
// * (/) = multiplication (division) membre à membre (et même sous membres à sous membre. pos_x*pos_x
// On peut aussi le faire avec un réel
//--------------------------------------------

Atome operator +(const Atome at1, const Atome at2)
{
    Atome sum;

    sum.set_pos(at1.get_pos()+at2.get_pos());
    sum.set_vel(at1.get_vel()+at2.get_vel());
    sum.set_acc(at1.get_acc()+at2.get_acc());
    sum.set_mass(at1.get_mass()+at2.get_mass());
    sum.set_charge(at1.get_charge()+at2.get_charge());
    sum.set_name(at1.get_name()+at2.get_name());

    return sum;
}




// ----------
// Comparison
// ----------

bool Atome::operator == (const Atome & at) const
{
    return (at.pos==pos && at.vel==vel && at.acc==acc && at.mass==mass && at.charge==charge && at.name == name) ? true : false;
}

bool Atome::operator != (const Atome & at) const
{
    return (at.pos!=pos || at.vel!=vel || at.acc!=acc || at.mass!=mass || at.charge!=charge || at.name != name) ? true : false;
}

//----------------------------------
// Surdéfinition des entrées sorties
//----------------------------------

ostream& operator << (ostream &flux,Atome at)
{
    at.write(flux);
    return(flux);
}

istream& operator >> (istream &flux,Atome & at) //at est modifié!
{
    at.read(flux);
    return(flux);
}
//----------
// Distances, angles
//----------

// Distance carrée entre deux Atomes
inline double Atome::dist2(const Atome & at) const
{
    double R2 = (pos - at.get_pos()).mag2();
    return(R2);
}

// Distance entre deux Atomes
inline double Atome::dist(const Atome & at) const
{
    return(sqrt(this->dist2(at))); // ou return(sqrt(*this.dist2(at)));
}

// Distance carrée entre deux vitesses d'Atomes
inline double Atome::dist2_vel(const Atome & at) const
{
    double R2 = (vel - at.get_vel()).mag2();
    return(R2);
}

// Distance entre deux deux vitesses d'Atomes
inline double Atome::dist_vel(const Atome & at) const
{
    return(sqrt(this->dist2_vel(at)));
}

inline double Atome::cosTheta (const Atome & at) const
{
    return((pos-at.get_pos()).Vecteur3D::cosTheta());
}
inline double Atome::cos2Theta(const Atome & at) const
{
    return((pos-at.get_pos()).Vecteur3D::cos2Theta());
}
// cos and cos^2 of the angle between the two points and axe Oz (axe Oz est theta=0).

inline double Atome::theta(const Atome & at) const
{
    return((pos-at.get_pos()).Vecteur3D::theta());
}


// ----------------------------------
//   FIN DES FONCTIONS EN LIGNES
// ----------------------------------


