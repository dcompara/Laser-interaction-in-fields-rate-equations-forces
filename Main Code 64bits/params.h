/*
  Name: Parameter Class (Param) and Parameter List (FitParams)
  Description:
  This file defines:
  1. The `Param` class, which represents a named parameter with attributes like value,
     minimum/maximum bounds, step count, and temporal dependence.
  2. The `FitParams` class, which is a collection of `Param` objects, providing utilities
     for parameter management and manipulation.

  Key Features:
  - Supports parameter scanning and random parameter modifications.
  - Handles both time-dependent and fixed parameters.
  - Allows reading and initializing parameters from external files.

  Example Usage of `FitParams`:
    FitParams &params; // Vector of Param pointers.
    Param *pDayMax;
    pDayMax = params.LocateOrAddNonConstParam("DayMax");

  To set a parameter value:
    Param* param_to_set = params.LocateOrAddParam("ParamName");
    param_to_set->val = value;

  Looping through parameters:
    for (ParamIterator i = params.begin(); i != params.end(); ++i) {
        Param &p = **i;
        if (p.is_scanned) continue;
        // Perform operations on the parameter.
    }
*/

#ifndef PARAM__H
#define PARAM__H

#include <iostream>
#include <vector>
#include <string>
#include <gsl/gsl_rng.h>                   // For random number generator
#include <gsl/gsl_randist.h>               // For uniform random distribution

using namespace std;

//----------------------------------
// Parameter Class
//----------------------------------

/**
 * Represents a named parameter with attributes like value, bounds, and temporal behavior.
 */
class Param {
public:
    double val;            // Current value of the parameter.
    double minv;           // Minimum value.
    double maxv;           // Maximum value.
    double val_t0;         // Initial value at t = 0.
    double tau;            // Temporal variation rate.
    int nbstep;            // Number of steps for scanning.
    bool is_scanned;       // Indicates if the parameter is being scanned.
    bool is_time_dependent; // Indicates if the parameter changes over time.
    string name;           // Parameter name.

    // Constructors
    Param() : val(0.), minv(0.), maxv(0.), val_t0(0.), tau(1e100),
              nbstep(1), is_scanned(false), is_time_dependent(false), name("UNKNOWN") {}

    Param(const string &Name) : val(0.), minv(0.), maxv(0.), val_t0(0.), tau(1e100),
                                nbstep(1), is_scanned(false), is_time_dependent(false), name(Name) {}

    // Returns the parameter name.
    string Name() const { return name; }

    // Outputs the parameter's details.
    void dump(ostream &stream = cout, bool nicely_formatted = false) const;

    // Overloaded output operator for the parameter.
    friend ostream& operator<<(ostream &stream, const Param &p) {
        p.dump(stream);
        return stream;
    }
};

//----------------------------------
// Parameter List Class
//----------------------------------

/**
 * Represents a collection of parameters for fitting or simulation.
 * Provides utilities for parameter management.
 */
class FitParams : public vector<Param *> {
    string name; // Name of the parameter set.

public:
    // Constructors
    FitParams(const string &Name) : name(Name) {}
    FitParams() : name("NoName") {}
    FitParams(const FitParams &);

    // Sets the name of the parameter set.
    void SetName(const string &Name) { name = Name; }

    // Finds a parameter by name (const version).
    const Param* LocateParam(const string &ParamName) const;

    // Finds a parameter by name (modifiable version).
    Param* LocateParam(const string &ParamName);

    // Locates or adds a parameter to the list (const version).
    const Param* LocateOrAddParam(const string &ParamName);

    // Locates or adds a parameter to the list (modifiable version).
    Param* LocateOrAddNonConstParam(const string &ParamName);

    // Deletes a parameter by name if it exists.
    bool DeleteParamIfExists(const string &ParamName);

    // Returns the number of variable (scanned) parameters.
    size_t VariableParamCount() const;

    // Outputs the parameters and their details.
    void dump(ostream &stream = cout, const string &Message = "", bool nicely_formatted = false) const;

    // Outputs the parameter names and values in a model format.
    void dump_model(ostream &stream = cout, const string &Message = "") const;

    // Overloaded output operator for the parameter list.
    friend ostream& operator<<(ostream &stream, const FitParams &p) {
        p.dump(stream);
        return stream;
    }

    // Reads parameters from a file.
    void read(const string &FileName);

    // Initializes parameters based on external data.
    void init_value(const string &FileName);

    // Scans through the parameters.
    bool Scan_Param();

    // Randomly modifies parameters during a scan.
    bool Scan_Param_aleatoire(const gsl_rng *r);

    // Destructor
    ~FitParams() {
        for (FitParams::iterator i = begin(); i != end(); ++i) {
            delete *i;
        }
    }

    // Assignment operator for copying parameter lists.
    void operator=(const FitParams &);
};

//----------------------------------
// Helper Functions
//----------------------------------

// Iterators for looping through parameters.
typedef FitParams::iterator ParamIterator;
typedef FitParams::const_iterator ParamCIterator;

// Extracts parameter names from a file.
vector<string> FitParamNames(const string &FileName);

#endif /* PARAM__H */
