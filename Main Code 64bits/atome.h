/*
  Name:  classe « Atome »
  Copyright:
  Author: Daniel Comparat
  Date: 15/10/06 11:01
  Description:


EN FAIT C'EST L'ETAT EXTERNE D'UNE PARTICULE
 Classe Atome at contenant un vecteur3D pos, vel, acc and mass + charge  + name
 * Elles sont initialisée à 0.
 * On peut leur mettre des valeurs par (exemple) at.set_pos(new_pos)
 * On peut lire les valeurs par (exemple) at.get_pos()
 * On peut incrémenter les valeurs par  at.inc_pos(d_pos)
 * == et != sont surdéfinis ATTENTION ILS COMPARENT DES DOUBLES (avec les erreurs d'arrondis cela doit être faux)
 * write et read permettent d'écrire et de lire dans un flux (cout (pas cerr ou clog) ou fichier)
 * FONCTIONS
            at1.dist2(at2) est distance carrée entre at1 et at2
            at1.dist(at2) est distance  entre at1 et at2
            at1.cosTheta(at2)
            at1.cos2Theta(at2) cos and cos^2 of the angle between MM2 and axe Oz (axe Oz est theta=0).

  */




#ifndef Atom_SEEN
#define Atom_SEEN

#include <iostream>
using namespace std;

#include "Vecteur3D.h"

class Atome
{
protected:
    Vecteur3D pos;             // (X,Y,Z) position
    Vecteur3D vel;             // (X,Y,Z) velocity
    Vecteur3D acc;             // (X,Y,Z) acceleration
    double mass;              // mass
    double charge;            // coulombian charge
    string name;

public:                // Forme canonique d'une classe
    Atome();                        //constructeur
    Atome(const Atome&);            // Constructeur de copie
    Atome& operator = (const Atome&);       // Affectation par recopie
    virtual ~Atome() {};                       // Destructeur par defaut

//--------------------------------------------
// Surcharge des opérateurs +,-,*,/
// +  = addition  membre à membre.
// *  = multiplication  membre à membre (et même sous membres à sous membre. pos_x*pos_x
// On peut aussi le faire avec un réel
//--------------------------------------------


    friend Atome operator +(const Atome , const Atome);    // surcharge de l'opérateur +

    // Get components
     Vecteur3D  get_pos()  const
    {
        return pos;
    }
    Vecteur3D  get_vel()  const
    {
        return vel;
    }
    Vecteur3D  get_acc()  const
    {
        return acc;
    }
    double get_mass()   const
    {
        return mass;
    }
    double get_charge()   const
    {
        return charge;
    }
    string get_name()   const
    {
        return name;
    }


    // Set components
    void  set_pos(const Vecteur3D& new_pos)
    {
        pos = new_pos;
    }
    void  set_vel(const Vecteur3D& new_vel)
    {
        vel = new_vel;
    }
    void  set_acc(const Vecteur3D& new_acc)
    {
        acc = new_acc;
    }
    void  set_mass(const double& new_mass)
    {
        mass = new_mass;
    }
    void  set_charge(const double& new_charge)
    {
        charge = new_charge;
    }
      void  set_name(const string& new_name)
    {
        name = new_name;
    }

    // Clear components
    void  clear_pos()
    {
        pos = Vecteur3D(0.,0.,0.);
    }
    void  clear_vel()
    {
        vel = Vecteur3D(0.,0.,0.);
    }
    void  clear_acc()
    {
        acc = Vecteur3D(0.,0.,0.);
    }
    void  clear_mass()
    {
        mass = 0.;
    }
    void  clear_charge()
    {
        charge = 0.;
    }
// I did not put for name it is useless

    void  inc_pos(const Vecteur3D& d_pos)
    {
        pos += d_pos;
    }
    void  inc_vel(const Vecteur3D& d_vel)
    {
        vel += d_vel;
    }
    void  inc_acc(const Vecteur3D& d_acc)
    {
        acc += d_acc;
    }
// I did not put for mass, charge, name it is useless

    // Comparison
    bool operator == (const Atome & at) const;
    bool operator != (const Atome & at) const;


    // Affichage

    virtual void write(ostream & flux)
    {
        if (&flux == &cout)
            cout << "masse : " << mass << "\t";
        else
            flux << mass << "\t";
        if (&flux == &cout)
            cout << "charge : " << charge << "\t";
        else
            flux << charge << "\t";
        if (&flux == &cout)
            cout << "position : " << pos << "\t";
        else
            flux << pos << "\t";
        if (&flux == &cout)
            cout << "vitesse : " << vel << "\t";
        else
            flux << vel<< "\t";
        if (&flux == &cout)
            cout << "acceleration : " << acc << "\t";
        else
            flux << acc << "\t";
             if (&flux == &cout)
            cout << "name : " << name << "\t";
        else
            flux << name << "\t";
    }

    // read des données

    virtual void read(istream & flux)
    {
        if (&flux == &cin)
            cout << "enter mass : ";
        flux >> mass;
        if (&flux == &cin)
            cout << "enter charge : ";
        flux >> charge;
        if (&flux == &cin)
            cout << "enter position : ";
        flux >> pos;
        if (&flux == &cin)
            cout << "enter vitesse : ";
        flux >> vel;
        if (&flux == &cin)
            cout << "enter acceleration : ";
        flux >> acc;
           if (&flux == &cin)
            cout << "enter name : ";
        flux >> name;
    }


    // Distance carrée entre deux Atomes
    inline double dist2(const Atome & at) const;

    // Distance entre deux Atomes
    inline double dist(const Atome & at) const;

    // Distance carrée entre deux vitesses d'Atomes
    inline double dist2_vel(const Atome & at) const;

    // Distance entre deux deux vitesses d' Atomes
    inline double dist_vel(const Atome & at) const;

    inline double cosTheta (const Atome & at) const;

    inline double cos2Theta(const Atome & at) const;

    // cos and cos^2 of the angle between MM2 and axe Oz (axe Oz est theta=0).

    inline double theta(const Atome & at) const;

};


#endif
