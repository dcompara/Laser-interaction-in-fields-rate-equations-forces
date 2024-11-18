/*
 Name: Laser Class
  Author: Daniel Comparat
  Date: 15/12/08

  Description:
  Represents a Gaussian laser characterized by:
    - Waist position in 3D space (`waist_pos`).
    - Direction vector (`direction`).
    - Waist dimensions (`waist`) in the x, y, and z directions.
    - Wavelength (`lambda`) in meters (vacuum wavelength).
    - Spectral width (`Gamma_Laser`) in s^-1 (FWHM).
    - Power (`Power`).
    - Polarization vector (`polarisation`), encoding circular and linear components.
    - Polarization angle (`polar_angle_degree`) relative to the x-axis.
    - Laser type (`type_laser`), such as CW, femto, Gaussian, etc.
    - Coherence with other lasers (`coherent_avec_laser_num`).

  The class includes:
    - Spectrum handling (energy attenuation and temporal intensity).
    - Calculations for irradiance (intensity) at specific points in space.
    - Euler angles for coordinate transformations.
    - Initialization and manipulation of laser parameters.

  Notes:
    - Energy and dipole values are in SI units (e.g., cm^-1 for energy, Debye for dipole).
    - Default initialization avoids division by zero for critical parameters.

  */



#ifndef my_laser_SEEN
#define my_laser_SEEN

#include <map>          // For the laser spectrum (energy and relative intensity)
#include <iostream>
#include <complex>
#include "Vecteur3D.h"
#include "constantes_SI.h" // SI constants (e.g., pi, cm, MHz)

using namespace std;

// Laser types enumeration
enum LaserType
{
    spon_emission = -1, // No laser (spontaneous emission)
    CW = 0,             // Continuous wave
    femto = 1,          // Femtosecond pulse
    multi_mode = 2,     // Multimode laser
    pulse = 3,          // Pulsed laser
    gaussien = 5,       // Gaussian profile
    lorentzien = 6,     // Lorentzian profile
    comb = 7,           // Frequency comb
    pseudo_BBR = 8,     // Pseudo black body radiation
    field_ionization = 9, // Field ionization
    pseudo_collimated_top_hat = 10 // Top-hat beam (flat-topped profile)
};


class Laser
{
protected:
    // Laser parameters
    Vecteur3D waist_pos;            // Position of the waist (X, Y, Z)
    Vecteur3D direction;            // Direction of the wave vector
    Vecteur3D waist;                // Waist dimensions in the X, Y, and Z directions
    double lambda;                  // Wavelength in meters (vacuum wavelength)
    double Gamma_Laser;             // Spectral width (FWHM)
    double Power;                   // Power in Watts
    Vecteur3D polarisation;         // Polarization vector (-1, 0, +1 for sigma-, pi, sigma+)
    double polar_angle_degree;      // Polarization angle relative to the x-axis
    int type_laser;                 // Laser type (e.g., CW, femto)
    int coherent_avec_laser_num;    // Coherent laser index (-1 if none)

public:
    // Laser spectrum and temporal attenuation
    map<double, double> spectre_Ecm_attenuation;  // Energy attenuation spectrum (energy in cm^-1, attenuation factor)
    map<double, double> Intensity_time_attenuation; // Intensity vs time (time in nanoseconds)

// Canonical class form
    Laser();                                 // Constructor
    Laser(const Laser&);                     // Copy constructor
    Laser& operator=(const Laser&);          // Assignment operator
    virtual ~Laser();                        // Destructor

// Getters
    Vecteur3D get_waist_pos() const
    {
        return waist_pos;
    }
    Vecteur3D get_direction() const
    {
        return direction.unit();    // Normalized direction vector
    }
    Vecteur3D get_waist() const
    {
        return waist;
    }
    double get_lambda() const
    {
        return lambda;
    }
    double get_omega() const
    {
        return 2 * pi * C / lambda;    // Angular frequency
    }
    double get_Gamma_Laser() const
    {
        return Gamma_Laser;
    }
    double get_Power() const
    {
        return Power;
    }
    Vecteur3D get_polarisation() const
    {
        return polarisation;
    }
    double get_polar_angle_degree() const
    {
        return polar_angle_degree;
    }
    int get_type_laser() const
    {
        return type_laser;
    }
    int get_coherent_avec_laser_num() const
    {
        return coherent_avec_laser_num;
    }

// Setters
    void set_waist_pos(const Vecteur3D& new_pos)
    {
        waist_pos = new_pos;
    }
    void set_direction(const Vecteur3D& new_dir)
    {
        direction = new_dir;
    }
    void set_waist(const Vecteur3D& new_waist)
    {
        waist = new_waist;
    }
    void set_lambda(const double& new_lambda)
    {
        lambda = new_lambda;
    }
    void set_energy_cm(const double& new_energy_cm)
    {
        lambda = 1. / (100. * new_energy_cm);
    }
    void set_Gamma_Laser(const double& new_Gamma_Laser)
    {
        Gamma_Laser = new_Gamma_Laser;
    }
    void set_Gamma_Laser_MHz(const double& new_Gamma_Laser_MHz)
    {
        Gamma_Laser = 2 * pi * new_Gamma_Laser_MHz * 1e6;
    }
    void set_Power(const double& new_Power)
    {
        Power = new_Power;
    }
    void set_polarisation(const Vecteur3D& new_pol)
    {
        polarisation = new_pol;
    }
    void set_polar_angle(const double& new_angle)
    {
        polar_angle_degree = new_angle;
    }
    void set_type_laser(const int& new_type)
    {
        type_laser = new_type;
    }
    void set_coherent_avec_laser_num(const int& num)
    {
        coherent_avec_laser_num = num;
    }


// Spectrum and intensity management

    /**
     * Reads intensity versus time from a stream.
     * @param flux Input stream containing intensity data.
     */
    void read_Intensity(istream& flux);

    /**
     * Reads intensity versus time from a file.
     * @param filename Path to the file containing intensity data.
     * @return 0 on success, non-zero on failure.
     */
    int read_Intensity(const char* filename);

    /**
     * Reads the laser spectrum (energy and attenuation) from a stream.
     * @param flux Input stream containing spectrum data.
     */
    void read_Spectrum(istream& flux);

    /**
     * Reads the laser spectrum (energy and attenuation) from a file.
     * @param filename Path to the file containing spectrum data.
     * @return 0 on success, non-zero on failure.
     */
    int read_Spectrum(const char* filename);

    /**
     * Writes the laser spectrum to a stream.
     * @param flux Output stream for spectrum data.
     */
    void write_Spectrum(ostream& flux);


// Utility methods for clearing attributes
    void clear_waist_pos()
    {
        waist_pos = Vecteur3D(0., 0., 0.);
    }
    void clear_direction()
    {
        direction = Vecteur3D(0., 0., 1.);
    }
    void clear_waist()
    {
        waist = Vecteur3D(1., 1., 0.);
    }
    void clear_lambda()
    {
        lambda = 1.;
    }
    void clear_Gamma_Laser()
    {
        Gamma_Laser = 1.;
    }
    void clear_Power()
    {
        Power = 0.;
    }
    void clear_polarisation()
    {
        polarisation = Vecteur3D(0., 0., 0.);
    }
    void clear_polar_angle()
    {
        polar_angle_degree = 0.;
    }
    void clear_type_laser()
    {
        type_laser = CW;
    }
    void clear_coherent_avec_laser_num()
    {
        coherent_avec_laser_num = -1;
    }


// Comparison operators
    bool operator==(const Laser& other) const;
    bool operator!=(const Laser& other) const;


    // Affichage

    virtual void write(ostream & flux) const
    {
        if (&flux == &cout)
            cout << "waist position : " << waist_pos << "\t";
        else
            flux << waist_pos << "\t";
        if (&flux == &cout)
            cout << "direction : " << direction << "\t";
        else
            flux << direction << "\t";
        if (&flux == &cout)
            cout << "waist : " << waist << "\t";
        else
            flux << waist<< "\t";
        if (&flux == &cout)
            cout << "lambda : " << lambda << "\t";
        else
            flux << lambda  << "\t";
        if (&flux == &cout)
            cout << "Gamma_Laser : " << Gamma_Laser << "\t";
        else
            flux << Gamma_Laser << "\t";
        if (&flux == &cout)
            cout << "Power : " << Power << "\t";
        else
            flux << Power << "\t";
        if (&flux == &cout)
            cout << "polarisation : " << polarisation << "\t";
        else
            flux << polarisation << "\t";
        if (&flux == &cout)
            cout << "polar_angle_degree : " << polar_angle_degree << "\t";
        else
            flux << polar_angle_degree << "\t";
        if (&flux == &cout)
            cout << "type_laser : " << type_laser << "\t";
        else
            flux << type_laser << "\t";
        if (&flux == &cout)
            cout << "coherent_avec_laser_num : " << coherent_avec_laser_num << "\t";
        else
            flux << coherent_avec_laser_num << "\t";

    }



    // read des données

    virtual void read(istream & flux)
    {
        if (&flux == &cin)
            cout << "entrez waist position : ";
        flux >> waist_pos;
        if (&flux == &cin)
            cout << "entrez direction : ";
        flux >> direction;
        if (&flux == &cin)
            cout << "entrez waist : ";
        flux >> waist;
        if (&flux == &cin)
            cout << "entrez lambda : ";
        flux >> lambda;
        if (&flux == &cin)
            cout << "entrez Gamma_Laser : ";
        flux >> Gamma_Laser;
        if (&flux == &cin)
            cout << "entrez Power : ";
        flux >> Power;
        if (&flux == &cin)
            cout << "entrez polarisation : ";
        flux >> polarisation;
        if (&flux == &cin)
            cout << "entrez polar_angle_degree : ";
        flux >> polar_angle_degree;
        if (&flux == &cin)
            cout << "entrez type_laser : ";
        flux >> type_laser;
        if (&flux == &cin)
            cout << "entrez coherent_avec_laser_num : ";
        flux >> coherent_avec_laser_num;

    }

// ----------------------------------
//   FONCTIONS EN LIGNES
// ----------------------------------


    // Intensity at time t. Linear Interpolated between the time given in the laser_intensity file. For time longer than the last time in the file the attenuation keep the last value
    double  intensity_t_nanosecond(const double t_nanosecond)  const;

    // transmission (1 = 100%, 0 = 0) prenant en compte le spectre à l'énergie de la transition
    double transmission_spectrum(const double energie_trans_cm)  const;


// waist au point (X,Y,Z). Used only for plot (so we did not improved it for top-hat:super gaussian)
    Vecteur3D  waist_size(const Vecteur3D& point)  const;

    /**
     * Calculates the irradiance at the laser waist.
     * @return Irradiance (intensity) at the waist.
     */
    double intensity() const;

    /**
     * Calculates the irradiance at a given point in the lab coordinate system.
     * @param point Point in the lab system where intensity is calculated.
     * @return Irradiance (intensity) at the specified point.
     */
    double intensity_lab_axis(const Vecteur3D& point) const;



// Inline functions for intensity calculations
    /**
     * Calculates the intensity at the waist in the laser coordinate system.
     * @return Intensity at the laser waist.
     */
    double intensity_repere_sur_waist(const Vecteur3D& point) const;

    /**
     * Computes the wave vector of the laser.
     * @return Wave vector as a Vecteur3D.
     */
    Vecteur3D wave_vector() const;

    /**
     * Computes the energy of the laser transition in cm^-1.
     * @return Energy in cm^-1.
     */
    double Energy_transition_laser_cm() const;




// ----------------------------------
//   FIN DES FONCTIONS EN LIGNES
// ----------------------------------

}
;



/***

Compare to the PRA 2014 we change notation now to be like Wikipedia in ZXZ convention(BE CAREFUL MATHEMATICA and Varshalovitch are in ZYZ convention)
 repère x,y,z du labo  et X,Y,Z de l'axe de quantification donné par le champ extérieur local
On utilise les angles d'Euler pour faire les rotations de repère
http://en.wikipedia.org/wiki/Euler_angles qui note (alpha,beta,gamma)


1. the first rotation is by an angle alpha (modulo 2pi, we chooose  (-pi,pi] to use atan2 ) about the z-axis,

2. the second rotation is by an angle beta in [0,pi] about the former (new) x-axis (now x')

3. the third rotation is by an angle gamma in (modulo 2pi, we chooose  (-pi,pi] to use atan2 ) about the former z-axis (now z').

***/

// Euler angles and coordinate transformations

/**
 * Calculates the Euler angles for a given direction.
 * @param direction Direction vector for the laser.
 * @return Euler angles (alpha, beta, gamma) as a Vecteur3D.
 */
Vecteur3D Euler_angles(Vecteur3D direction);

/**
 * Transforms coordinates from the lab frame to the laser frame.
 * @param point Point in lab coordinates.
 * @param alpha First Euler angle (rotation about z-axis).
 * @param beta Second Euler angle (rotation about x-axis).
 * @param gamma Third Euler angle (rotation about z'-axis, default 0).
 * @return Point in the laser frame.
 */
Vecteur3D rotation_lab_axis(const Vecteur3D& point, double alpha, double beta, double gamma = 0);

/**
 * Transforms coordinates from the laser frame to the lab frame.
 * @param point Point in laser coordinates.
 * @param alpha First Euler angle (rotation about z-axis).
 * @param beta Second Euler angle (rotation about x-axis).
 * @param gamma Third Euler angle (rotation about z'-axis, default 0).
 * @return Point in the lab frame.
 */
Vecteur3D rotation_axis_lab(const Vecteur3D& point, double alpha, double beta, double gamma = 0);




// Additional mathematical utilities


// I = ε0 c E^2 /2.
double champ_E(const double irradiance);


// (absolute value of the)  effectif dipole d.e_laser = sum_p d_p epsilon^p
// where the dipole transition vector d= sum_p d_p e^p is given in the local quantification axis
// and the polarisation vector e_laser= sum_p' epsilon^p' e_p'  is given in the laser axis
double effectif_dipole_local(const complex<double> dipole[3], const Vecteur3D& axe_quant,  const Laser& my_laser);



/**
 * Approximates the arctangent on the range [-1, 1].
 * @param z Input value.
 * @return Approximated arctangent.
 */
float ApproxAtan(float z);

/**
 * Fast approximation of atan2 for efficiency in repeated calculations.
 * @param y Y-coordinate.
 * @param x X-coordinate.
 * @return Approximated atan2 value.
 */
float atan2_approximation(float y, float x);



#endif
