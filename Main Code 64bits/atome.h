/*
  Name:  Class "Atome"
  Author: Daniel Comparat
  Date: 15/10/06 11:01
  Description:

This class represents the external state of a particle.
Class `Atome` contains a 3D vector for position, velocity, and acceleration,
as well as properties like mass, charge, and name.
 * Properties are initialized to 0.
 * Values can be set using, e.g., `at.set_pos(new_pos)`.
 * Values can be read using, e.g., `at.get_pos()`.
 * Values can be incremented using, e.g., `at.inc_pos(d_pos)`.
 * The `==` and `!=` operators are overridden but note that they compare doubles,
   which may cause inaccuracies due to rounding errors.
 * The `write` and `read` methods allow writing to and reading from a stream
   (e.g., `cout` or a file).

FUNCTIONS:
 - `at1.dist2(at2)` returns the squared distance between `at1` and `at2`.
 - `at1.dist(at2)` returns the distance between `at1` and `at2`.
 - `at1.cosTheta(at2)` computes the cosine of the angle between MM2 and the z-axis.
 - `at1.cos2Theta(at2)` computes the square of the cosine of the angle between MM2 and the z-axis.
*/

#ifndef Atom_SEEN
#define Atom_SEEN

#include <iostream>
#include "Vecteur3D.h"  // Handles 3D vectors

class Atome
{
protected:
    Vecteur3D pos;             // (X, Y, Z) position
    Vecteur3D vel;             // (X, Y, Z) velocity
    Vecteur3D acc;             // (X, Y, Z) acceleration
    double mass;               // Mass of the particle
    double charge;             // Coulomb charge
    std::string name;          // Particle name

public:
    // Canonical form of a class
    Atome();                               // Constructor
    Atome(const Atome&);                   // Copy constructor
    Atome& operator=(const Atome&);        // Copy assignment operator
    virtual ~Atome() {};                   // Virtual destructor

    // Overloaded operators
    friend Atome operator+(const Atome, const Atome);  // Overloaded addition operator

    // Getters
    Vecteur3D get_pos() const { return pos; }
    Vecteur3D get_vel() const { return vel; }
    Vecteur3D get_acc() const { return acc; }
    double get_mass() const { return mass; }
    double get_charge() const { return charge; }
    std::string get_name() const { return name; }

    // Setters
    void set_pos(const Vecteur3D& new_pos) { pos = new_pos; }
    void set_vel(const Vecteur3D& new_vel) { vel = new_vel; }
    void set_acc(const Vecteur3D& new_acc) { acc = new_acc; }
    void set_mass(const double& new_mass) { mass = new_mass; }
    void set_charge(const double& new_charge) { charge = new_charge; }
    void set_name(const std::string& new_name) { name = new_name; }

    // Clear properties
    void clear_pos() { pos = Vecteur3D(0., 0., 0.); }
    void clear_vel() { vel = Vecteur3D(0., 0., 0.); }
    void clear_acc() { acc = Vecteur3D(0., 0., 0.); }
    void clear_mass() { mass = 0.; }
    void clear_charge() { charge = 0.; }

    // Increment properties
    void inc_pos(const Vecteur3D& d_pos) { pos += d_pos; }
    void inc_vel(const Vecteur3D& d_vel) { vel += d_vel; }
    void inc_acc(const Vecteur3D& d_acc) { acc += d_acc; }

    // Comparison operators
    bool operator==(const Atome& at) const;
    bool operator!=(const Atome& at) const;

    // Output the atom's data
    virtual void write(std::ostream& flux)
    {
        flux << "mass: " << mass << "\t"
             << "charge: " << charge << "\t"
             << "position: " << pos << "\t"
             << "velocity: " << vel << "\t"
             << "acceleration: " << acc << "\t"
             << "name: " << name << "\t";
    }

    // Read atom data from a stream
    virtual void read(std::istream& flux)
    {
        if (&flux == &std::cin) std::cout << "Enter mass: ";
        flux >> mass;

        if (&flux == &std::cin) std::cout << "Enter charge: ";
        flux >> charge;

        if (&flux == &std::cin) std::cout << "Enter position: ";
        flux >> pos;

        if (&flux == &std::cin) std::cout << "Enter velocity: ";
        flux >> vel;

        if (&flux == &std::cin) std::cout << "Enter acceleration: ";
        flux >> acc;

        if (&flux == &std::cin) std::cout << "Enter name: ";
        flux >> name;
    }

    // Distances and angles
    inline double dist2(const Atome& at) const;
    inline double dist(const Atome& at) const;
    inline double dist2_vel(const Atome& at) const;
    inline double dist_vel(const Atome& at) const;
    inline double cosTheta(const Atome& at) const;
    inline double cos2Theta(const Atome& at) const;
    inline double theta(const Atome& at) const;  // Cosine and squared cosine of angles
};

#endif
