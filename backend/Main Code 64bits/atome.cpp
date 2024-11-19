#include "atome.h"

// ------------
// Constructors
// ------------

Atome::Atome()  // Default constructor
{
    pos.set(0., 0., 0.);
    vel.set(0., 0., 0.);
    acc.set(0., 0., 0.);
    mass = 0.;
    charge = 0.;
    name = "";
}

Atome::Atome(const Atome &at)  // Copy constructor
{
    pos = at.pos;
    vel = at.vel;
    acc = at.acc;
    mass = at.mass;
    charge = at.charge;
    name = at.name;
}

Atome &Atome::operator=(const Atome &at)  // Copy assignment operator
{
    if (this != &at)  // Check for self-assignment
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
// Overloaded operators
// + and - = member-wise addition and subtraction.
// * and / = member-wise multiplication and division.
//--------------------------------------------

Atome operator+(const Atome at1, const Atome at2)
{
    Atome sum;

    sum.set_pos(at1.get_pos() + at2.get_pos());
    sum.set_vel(at1.get_vel() + at2.get_vel());
    sum.set_acc(at1.get_acc() + at2.get_acc());
    sum.set_mass(at1.get_mass() + at2.get_mass());
    sum.set_charge(at1.get_charge() + at2.get_charge());
    sum.set_name(at1.get_name() + at2.get_name());

    return sum;
}

// ----------
// Comparison
// ----------

bool Atome::operator==(const Atome &at) const
{
    return (at.pos == pos && at.vel == vel && at.acc == acc && at.mass == mass && at.charge == charge && at.name == name);
}

bool Atome::operator!=(const Atome &at) const
{
    return !(at == *this);
}

//----------------------------------
// Input and Output Overloading
//----------------------------------

ostream &operator<<(ostream &flux, Atome at)
{
    at.write(flux);
    return flux;
}

istream &operator>>(istream &flux, Atome &at)
{
    at.read(flux);
    return flux;
}

// ----------
// Distances and Angles
// ----------

// Squared distance between two atoms
inline double Atome::dist2(const Atome &at) const
{
    return (pos - at.get_pos()).mag2();
}

// Distance between two atoms
inline double Atome::dist(const Atome &at) const
{
    return sqrt(this->dist2(at));
}

// Squared distance between velocities of two atoms
inline double Atome::dist2_vel(const Atome &at) const
{
    return (vel - at.get_vel()).mag2();
}

// Distance between velocities of two atoms
inline double Atome::dist_vel(const Atome &at) const
{
    return sqrt(this->dist2_vel(at));
}

// Cosine of the angle between positions and z-axis
inline double Atome::cosTheta(const Atome &at) const
{
    return (pos - at.get_pos()).Vecteur3D::cosTheta();
}

// Cosine squared of the angle between positions and z-axis
inline double Atome::cos2Theta(const Atome &at) const
{
    return (pos - at.get_pos()).Vecteur3D::cos2Theta();
}

// Angle between positions and z-axis
inline double Atome::theta(const Atome &at) const
{
    return (pos - at.get_pos()).Vecteur3D::theta();
}


