<<<<<<< HEAD
/*
  Name: Molecule Class
  Author: Daniel Comparat
  Date: 15/10/06 (Simplified on 12/02/2012)

  Description:
  This class represents a molecule and derives from the `Atome` class.
  It encapsulates:
    1) External state: Inherited from the `Atome` class (includes position, velocity, and mass).
    2) Internal state: Inherited from the `Internal_state` class (includes electronic level, vibrational level, angular momentum J, and projection MJ).

  Features:
  - Combines functionalities of `Atome` and `Internal_state`.
  - Includes operators for assignment and addition.
  - Supports input/output operations and potential calculations.

  Additional Notes:
  - Electronic state levels are coded as integers:
      - Negative values for "dead" states (e.g., annihilation, photoionization).
      - 0 for ground states, 1 for excited states, etc.
*/

#ifndef Molecule_SEEN
#define Molecule_SEEN

#include "Atome.h"
#include "Internal_state.h"
#include "Field.h"
#include <complex>
#include <iostream>
using namespace std;

// Reaction coding structure
struct type_codage_react {
    int n_mol;                     // Molecule identifier
    int n_laser;                   // Laser identifier
    Vecteur3D quant_axis;          // Quantization axis
    complex<double> pol_vector[3]; // Polarization vector
    Vecteur3D k_eff_laser;         // Effective wave vector of the laser
    Internal_state final_internal_state; // Final internal state
};

/*
  Notes on `type_codage_react`:
  - Encodes a reaction based on the initial state (molecule index) and final state.
  - For stimulated processes, includes the effective wave vector and polarization vector.
  - Aims to model recoil effects and emission probabilities.
  - Optimization suggestion: Use references or pointers for `final_internal_state` to reduce memory usage.
*/

// Electronic state levels
const int annihilation = -3; // Annihilated particles
const int photodetached = -2;
const int photoionized = -1; // Photoionized level
const int fond = 0;          // Ground state
const int excite = 1;        // First excited state
const int excite2 = 2;       // Second excited state
const int excite3 = 3;       // Third excited state

// Molecule class definition
class Molecule : public Atome, public Internal_state {
public:
    // Constructors and Destructor
    Molecule();                         // Default constructor
    ~Molecule();                        // Destructor

    // Assignment operator for internal states
    Molecule& operator=(const Internal_state& Intern);

    // Overloaded addition operator
    friend Molecule operator+(const Molecule, const Molecule);

    // Input/Output Methods
    /**
     * Writes the molecule's state (external and internal) to the provided stream.
     * @param flux Output stream for writing data.
     */
    void write(ostream& flux) {
        this->Atome::write(flux);       // Write external state
        this->Internal_state::write(flux); // Write internal state
    }

    /**
     * Writes the list of possible transitions (internal state lines) to the stream.
     * @param flux Output stream for writing data.
     */
    void write_list(ostream& flux) {
        this->Internal_state::write_liste_raie(flux);
    }

    /**
     * Reads the molecule's state (external and internal) from the provided stream.
     * @param flux Input stream for reading data.
     */
    void read(istream& flux) {
        this->Atome::read(flux);        // Read external state
        this->Internal_state::read(flux); // Read internal state
    }
};

/**
 * Computes the potential energy of a molecule in external fields.
 * Includes contributions from electric and magnetic fields.
 * @param mol Molecule whose potential is calculated.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @return Potential energy in Joules.
 */
double get_pot(const Molecule& mol, const Field& fieldB, const Field& fieldE);

/**
 * Computes the kinetic energy of a molecule.
 * @param mol Molecule whose kinetic energy is calculated.
 * @return Kinetic energy in Joules.
 */
double get_kin(const Molecule& mol);

//----------------------------------
// Stream Overloading for Molecule
//----------------------------------

/**
 * Writes the molecule's state to the provided output stream.
 * @param flux Output stream.
 * @param my_molecule Molecule to write.
 * @return The updated output stream.
 */
ostream& operator<<(ostream& flux, Molecule my_molecule);

/**
 * Reads the molecule's state from the provided input stream.
 * @param flux Input stream.
 * @param my_molecule Molecule to populate with data.
 * @return The updated input stream.
 */
istream& operator>>(istream& flux, Molecule& my_molecule);

#endif
=======
/*
  Name: Molecule Class
  Author: Daniel Comparat
  Date: 15/10/06 (Simplified on 12/02/2012)

  Description:
  This class represents a molecule and derives from the `Atome` class.
  It encapsulates:
    1) External state: Inherited from the `Atome` class (includes position, velocity, and mass).
    2) Internal state: Inherited from the `Internal_state` class (includes electronic level, vibrational level, angular momentum J, and projection MJ).

  Features:
  - Combines functionalities of `Atome` and `Internal_state`.
  - Includes operators for assignment and addition.
  - Supports input/output operations and potential calculations.

  Additional Notes:
  - Electronic state levels are coded as integers:
      - Negative values for "dead" states (e.g., annihilation, photoionization).
      - 0 for ground states, 1 for excited states, etc.
*/

#ifndef Molecule_SEEN
#define Molecule_SEEN

#include "Atome.h"
#include "Internal_state.h"
#include "Field.h"
#include <complex>
#include <iostream>
using namespace std;

// Reaction coding structure
struct type_codage_react {
    int n_mol;                     // Molecule identifier
    int n_laser;                   // Laser identifier
    Vecteur3D quant_axis;          // Quantization axis
    complex<double> pol_vector[3]; // Polarization vector
    Vecteur3D k_eff_laser;         // Effective wave vector of the laser
    Internal_state final_internal_state; // Final internal state
};

/*
  Notes on `type_codage_react`:
  - Encodes a reaction based on the initial state (molecule index) and final state.
  - For stimulated processes, includes the effective wave vector and polarization vector.
  - Aims to model recoil effects and emission probabilities.
  - Optimization suggestion: Use references or pointers for `final_internal_state` to reduce memory usage.
*/

// Electronic state levels
const int annihilation = -3; // Annihilated particles
const int photodetached = -2;
const int photoionized = -1; // Photoionized level
const int fond = 0;          // Ground state
const int excite = 1;        // First excited state
const int excite2 = 2;       // Second excited state
const int excite3 = 3;       // Third excited state

// Molecule class definition
class Molecule : public Atome, public Internal_state {
public:
    // Constructors and Destructor
    Molecule();                         // Default constructor
    ~Molecule();                        // Destructor

    // Assignment operator for internal states
    Molecule& operator=(const Internal_state& Intern);

    // Overloaded addition operator
    friend Molecule operator+(const Molecule, const Molecule);

    // Input/Output Methods
    /**
     * Writes the molecule's state (external and internal) to the provided stream.
     * @param flux Output stream for writing data.
     */
    void write(ostream& flux) {
        this->Atome::write(flux);       // Write external state
        this->Internal_state::write(flux); // Write internal state
    }

    /**
     * Writes the list of possible transitions (internal state lines) to the stream.
     * @param flux Output stream for writing data.
     */
    void write_list(ostream& flux) {
        this->Internal_state::write_liste_raie(flux);
    }

    /**
     * Reads the molecule's state (external and internal) from the provided stream.
     * @param flux Input stream for reading data.
     */
    void read(istream& flux) {
        this->Atome::read(flux);        // Read external state
        this->Internal_state::read(flux); // Read internal state
    }
};

/**
 * Computes the potential energy of a molecule in external fields.
 * Includes contributions from electric and magnetic fields.
 * @param mol Molecule whose potential is calculated.
 * @param fieldB Magnetic field configuration.
 * @param fieldE Electric field configuration.
 * @return Potential energy in Joules.
 */
double get_pot(const Molecule& mol, const Field& fieldB, const Field& fieldE);

/**
 * Computes the kinetic energy of a molecule.
 * @param mol Molecule whose kinetic energy is calculated.
 * @return Kinetic energy in Joules.
 */
double get_kin(const Molecule& mol);

//----------------------------------
// Stream Overloading for Molecule
//----------------------------------

/**
 * Writes the molecule's state to the provided output stream.
 * @param flux Output stream.
 * @param my_molecule Molecule to write.
 * @return The updated output stream.
 */
ostream& operator<<(ostream& flux, Molecule my_molecule);

/**
 * Reads the molecule's state from the provided input stream.
 * @param flux Input stream.
 * @param my_molecule Molecule to populate with data.
 * @return The updated input stream.
 */
istream& operator>>(istream& flux, Molecule& my_molecule);

#endif
>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
