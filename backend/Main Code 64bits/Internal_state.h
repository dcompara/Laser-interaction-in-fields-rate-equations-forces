/*
  Class: Internal_state
  Author: Daniel Comparat
  Date: 15/10/06 11:01

  Description:
  This class represents the internal state of a molecule or particle, characterized by:
    - Electronic level (manifold): exc
    - Total angular momentum (J or F)
    - Projection of angular momentum (M)
    - Additional variables such as vibrational levels (v)
    - Bound state (1 for bound, 0 for unbound)
    - Degeneracy index to lift level degeneracies
    - Lifetime and energy shifts
    - Population and transition coefficients

  Key Features:
    - Initialization of all parameters to zero (except bound level and degeneracy set to 1).
    - Ability to write and read data to/from streams (e.g., files, console).
    - Handles energy shifts under external fields (electric and magnetic).
    - Supports transitions to other states and manages their properties.

  Units:
    - Dipole moments: DEBYE
    - Energy: cm^-1

  Quantum Numbers:
    - Quantum numbers (e.g., M) are typically stored as integers multiplied by 2 for precision.
      For example, M = 1/2 is stored as 1, hence the term `two_M`.

  */

#ifndef Internal_state_SEEN
#define Internal_state_SEEN

#include "Vecteur3D.h"
#include "constantes_SI.h"
#include <vector>
#include <iostream>
using namespace std;

class Internal_state {
public:
    // State properties
    int exc;                       // Electronic level (manifold)
    int bound_level;               // Bound state indicator (1 = bound, 0 = unbound)
    int two_J;                     // 2 * Total angular momentum (J or F)
    int two_M;                     // 2 * Projection of total angular momentum
    int v;                         // Additional variable (e.g., vibrational number)
    int deg_number;                // Degeneracy index
    double Energy0_cm;             // Base energy in cm^-1
    double Delta_FieldB;           // Coefficient for magnetic field energy shift
    double Linear_FieldB;          // Linear term for magnetic field shift
    double Delta_FieldE;           // Coefficient for electric field energy shift
    double Linear_FieldE;          // Linear term for electric field shift
    double one_over_lifetime;      // Inverse of lifetime (Einstein A coefficient)
    double Energy_cm;              // Actual energy considering field shifts
    double population;             // Initial population of the state

    // Transition structure: Pair of target state pointer and dipole strength
    typedef pair<Internal_state*, double> transition;
    vector<transition> liste_raies; // List of possible transitions (target states and strengths)

    // Constructors and destructors
    Internal_state();                             // Default constructor
    Internal_state(const Internal_state& Intern); // Copy constructor
    ~Internal_state();                            // Destructor
    Internal_state& operator=(const Internal_state& Intern); // Assignment operator

    // Input/Output Methods
    void write(ostream& flux);                    // Outputs state information
    void read(istream& flux);                     // Reads state information
    void write_liste_raie(ostream& flux);         // Writes the list of transitions

    // State Comparison
    bool is_equal(const Internal_state& w) const; // Compares states (based on key properties)

    //----------------------------------
    // File I/O for Levels
    //----------------------------------
    void write_Level_B(ostream& flux);           // Writes state to stream in level format
    void read_Level__B(istream& flux);           // Reads level information for magnetic fields
    void read_Level__E(istream& flux);           // Reads level information for electric fields

    //----------------------------------
    // Field Shift Calculations
    //----------------------------------
    double Energy_Shift_B_cm(const double B) const;  // Energy shift in magnetic field
    double Energy_Shift_E_cm(const double E) const;  // Energy shift in electric field
    Vecteur3D Grad_Energy_Shift_B(const double B, const Vecteur3D gradient_B2); // Gradient (SI units)
    Vecteur3D Grad_Energy_Shift_E(const double E, const Vecteur3D gradient_E2); // Gradient (SI units)
    double Field_Shift_B_cm(const double Eshift);    // Field strength for magnetic shift
    double Field_Shift_E_cm(const double Eshift);    // Field strength for electric shift

    //----------------------------------
    // Transition Management
    //----------------------------------
    void add_transition(Internal_state* v, const double strength); // Adds a transition
    double Einstein_desex_rate();                                  // Computes de-excitation rates
};

//----------------------------------
// Global Functions for Transitions
//----------------------------------
void read_Line(istream& flux, Internal_state& Upper_state, Internal_state& Lower_state,
               double& Dipole_Debye);                              // Reads transition from stream

void write_Line(ostream& flux, Internal_state& Upper_state, Internal_state& Lower_state,
                const double Dipole_Debye);                        // Writes transition to stream

//----------------------------------
// Overloaded Input/Output Operators
//----------------------------------
ostream& operator<<(ostream& flux, Internal_state Intern);         // Output operator
istream& operator>>(istream& flux, Internal_state& Intern);        // Input operator

#endif  // Internal_state_SEEN
