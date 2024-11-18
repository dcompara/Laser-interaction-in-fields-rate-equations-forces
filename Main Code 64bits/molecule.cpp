#include "molecule.h"



// ------------
// Constructors
// ------------

Molecule::Molecule()                    // Constructeur par défaut
{
};

Molecule::~Molecule()
{}
;                    // Destructeur

// opérateur d'affectation de l'état interne
// ATTENTION la copie de la liste est une copie seulement vers l'adresse du pointeur
// On en copie pas vraiment la liste. Cela implique que chaque molécule à la même liste!
Molecule & Molecule::operator = (const Internal_state & Intern)
{
    exc = Intern.exc;
    v = Intern.v;
    two_J = Intern.two_J;
    two_M = Intern.two_M;
    bound_level = Intern.bound_level;
    deg_number = Intern.deg_number;
    Energy0_cm = Intern.Energy0_cm;
    Delta_FieldB = Intern.Delta_FieldB;
    Linear_FieldB = Intern.Linear_FieldB;
    Delta_FieldE = Intern.Delta_FieldE;
    Linear_FieldE = Intern.Linear_FieldE;
    one_over_lifetime = Intern.one_over_lifetime;
    Energy_cm = Intern.Energy_cm;
    population = Intern.population;
    liste_raies = Intern.liste_raies;
    return *this;
}

//----------------------------------
// Surdéfinition des entrées sorties
//----------------------------------

ostream& operator << (ostream & flux, Molecule my_molecule)
{
    my_molecule.write(flux);
    return(flux);
}

istream& operator >> (istream &flux, Molecule & my_molecule) //my_molecule est modifiée
{
    my_molecule.read(flux);
    return(flux);
}

#include "shift_molecule.h"

// Lecture du potentiel [en Joule] de l'état interne d'une molécule (electric, magnétique) pas dipolaire
// require that mol.Energry_cm is correct
// On peut enlèver le pot au centre
// We include the Coulomb interaction with external field (not with other particle which is done in get_coulomb_potential)
double  get_pot(const Molecule &mol, const Field &fieldB, const Field &fieldE)
{
    Vecteur3D pos(0.,0.,0.);
    double Bcentre=fieldB.get_Field(pos).mag(); // Magnétique
    double Ecentre=fieldE.get_Field(pos).mag(); // Electrique
    // double E_pot_centre = HBAR*delta_field_shift_B_E(mol, Bcentre, Ecentre); // Si on veux enlever le potentiel au centre
    double E_pot_centre = 0.;
    double pot_electric = mol.get_charge()* fieldE.get_electric_potential(mol.get_pos()); // q V with FieldE= - Grad V
    double pot_gravity =   - mol.get_mass()*gravity.dot(mol.get_pos());
    double pot_internal = HBAR*(mol.Energy_cm - mol.Energy0_cm)*Conv_Ecm_delta;
    return pot_electric + pot_gravity + pot_internal - E_pot_centre;
}

// Lecture de l'énergie cinétique [en Joule]  d'une molécule
double  get_kin(const Molecule &mol)
{
    return   0.5*(mol.get_vel().mag2()) * mol.get_mass();
}
