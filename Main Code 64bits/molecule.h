/*
  Name:  classe « molecule »

  Copyright:
  Author: Daniel Comparat
  Date: 15/10/06 11:01
  SIMPLIFIE le 12/02/2012

  Description:
 Elle dérive de la classe Atome
mol un élément de la classe
contient
1) Son état externe: la classe Atome (vecteur3D pos, vel et la masse), elle ajoute
2) Son état interne: la classe Internal_state exc (niveau électronique), v,J, MJ

*/




#ifndef Molecule_SEEN
#define Molecule_SEEN

#include "Atome.h"
#include "Internal_state.h"
#include <complex>

#include <iostream>
using namespace std;

struct type_codage_react{ int n_mol; int n_laser; Vecteur3D quant_axis; complex<double> pol_vector[3]; Vecteur3D k_eff_laser; Internal_state final_internal_state;};
// In order to define a reaction we need to know the initial state (so the n_mol), the finale state (so final_internal_state)
// and to implement the recoil we need to have
// 1) For absorption or stimulated emission: the laser (n_laser) which will give the photon momentum.
// We also give k_eff_laser which is the laser (effectif) wave_vector. it can be zero in the case of lattice for instance (no momentum transfer)
// 2) for the spontaneous emission: the polarization vector (transition dipole vector) that will determined the emission probability of the emitted photon
// To help we also give the quantization axis where the dipole are calculated to then rotate to the lab axis if needed

// TODO (dc#1#): Put the FULL internal_state_finale is HUGE in term of memory. ...
//Study how to put only the adress !
//
// may be This is why the code is  slow ?!


/*** list of electronical states ***/
const int annihilation = -3; // Annihilized particles (P_bar or Ps for instance can annihilate)
const int photodetached = -2;
const int photoionized = -1; // Photoionized level
// So Negative values can be used for a "dead" level such as one in a continuum:
// usually 0 means ground electronical level,
// 1 is for an excited electronical level,
const int fond=0;
const int excite=1;
const int excite2=2;
const int excite3=3;


class Molecule : public Atome, public Internal_state
{

public:
    // vector < Internal_state *> liste_niveaux; //La liste des raies possibles est: vecteur donnant la liste des états finaux et de la force de transition

    Molecule();         // constructeurs autres que par recopie
    ~Molecule();           // destructeur

  Molecule & operator = (const Internal_state & Intern);  // opérateur d'affectation d'état interne



    friend Molecule operator +(const Molecule , const Molecule);    // surcharge de l'opérateur +

    // Affichage

    void write(ostream & flux)
    {
        this->Atome::write(flux);
        this->Internal_state::write(flux);
    }

    void write_list(ostream & flux)
    {
       this->Internal_state::write_liste_raie(flux);
    }

    // read des données

    void read(istream & flux)
    {
        this->Atome::read(flux);
        this->Internal_state::read(flux);

    }


};


#include "Field.h"

// Lecture du potentiel [en Joule] de l'état interne d'une molécule (electric, magnétique) pas dipolaire
// require that mol.Energry_cm is correct
// On peut enlèver le pot au centre
// We include the Coulomb interaction with external field (not with other particle which is done in get_coulomb_potential)
double  get_pot(const Molecule &mol, const Field &fieldB, const Field &fieldE);

// Lecture de l'énergie cinétique [en Joule]  d'une molécule
double  get_kin(const Molecule &mol);


//----------------------------------
// Surdéfinition des entrées sorties
//----------------------------------

ostream& operator << (ostream & flux, Molecule my_molecule);

istream& operator >> (istream &flux, Molecule & my_molecule);




#endif
