/*
  Name:  classe « Laser »
  Copyright:
  Author: Daniel Comparat
  Date: 15/12/08
  Description:

 Classe Laser
 Laser supposé gaussien
 contenant
waist_pos : vecteur3D position du waist dans le repère fixe du labo,
direction :  du vecteur d'onde, non nécessairement normalisé,
waist :  Le waist est w_kx, w_ky, w_kz où kz est l'axe selon k.
lambda : longueur d'onde lambda moyenne (en SI, i.e. en metre)  DANS LE VIDE.
Gamma_Laser : largeur spectrale FWHM  en s^-1.
Power : puissance,
Polarisation est un vecteur normé qui contient epsilon^-1 (codé sur x), espilon^0 (sur y) et epsilon^+1 (sur z)
où -1 (pour sigma- =  circulaire right car l'axe de quanti est k),0 (pour pi), 1 (pour sigma+ =  circulaire left car l'axe de quanti est k)
polar_angle_degree gives the polarization angle (cf User Guide)
type_laser : le type (CW, pulsé femto, gaussien, lorentzien, comb, black body...)
coherent_avec_laser_num: pour les interférences  est le numéro du premier laser avec lequel il est cohérent.


On met aussi (surtout pour le cas façonné) un spectre en énergie spectrale (cm^-1) du laser.
C'est un tableau Energie_cm, I_attenuation (1= non atténué, 0= éteind).
L'intensité sera donc celle du laser multiplié par I_attenuation à l'énergie (immédiatement supérieure) de la transition (1 par défaut).



Il y a aussi les fonctions donnant (ATTENTION LES X,Y,Z réfèrent au repère lié au laser !)
Donc pas dans le repère du labo!

Mais il y a une fonction qui donne les angles d'Euler:  Euler_angles

intensité (irradiance) au waist: intensity()
Zone de Rayleigh: Rayleigh_range()
waist au point (X,Y,Z): waist_size(point)
intensité (irradiance) au point (X,Y,Z):  intensity(point)
wave_vector:  wave_vector()

RAPPEL: E =h\nu = \hbar \omega= h c / \lambda_{\rm vide} = h c \sigma(m-1) = \hbar c k = 100 h c sigma(cm^-1)

 * Les éléments de la classe sont initialisée à 0. Sauf le waist, direction, lambda, Gamma_Laser à (1,1,0); (0,0,1), 1, 1 pour éviter des divisions par zéros
 * On peut leur mettre des valeurs par (exemple) my_laser.set_waist_pos(new_pos)
 * On peut lire les valeurs par (exemple) my_laser.get_waist_pos()
 * == et != sont surdéfinis ATTENTION ILS COMPARENT PARFOIS DES DOUBLES (avec les erreurs d'arrondis cela doit être faux)
 * write et read permettent d'écrire et de lire dans un flux (cout (pas cerr ou clog) ou fichier)



  */




#ifndef my_laser_SEEN
#define my_laser_SEEN

#include <map>          // Pour le spectre du laser Energy, Intensité relative
#include <iostream>
using namespace std;

#include "Vecteur3D.h"
#include "constantes_SI.h"                  // SI Constantes pour pi, cm, MHz ou autre

/*
Forme canonique d'une classe T (Claude Delannoy, Programmer en C++)

class T
{
      public:
             T(...);         // constructeur de T, autres que par recopie
             T(const T &);   // constructeur de recopie de T
             ~T();           // destructeur
             T & T::operator = (const T &);  // opérateur d'affectation
             .....
};
*/


enum // Différents types de lasers
{
    spon_emission = -1, //  no_laser
    CW = 0,
    femto =1,
    multi_mode=2,
    pulse=3,
    faconne = 4, // Obsolete car il y a le spectre maintenant
    gaussien = 5,
    lorentzien = 6,
    comb = 7,
    pseudo_BBR =8,
    field_ionization =9
};



class Laser
{
protected:
    Vecteur3D waist_pos;              // (X,Y,Z) position
    Vecteur3D direction;        // (X,Y,Z) direction du vecteur d'onde k
    Vecteur3D waist;            // Le waist est w_kx, w_ky, w_kz où kz est l'axe selon k.
    double lambda;              // longueur d'onde lambda moyenne (en SI, i.e. en metre), DANS LE VIDE
    double Gamma_Laser;            // largeur spectrale (Delta lambda)
    double Power;                // puissance,
    Vecteur3D polarisation;           // Polarisation -1,0,+1: sigma-,pi,sigma+
    double  polar_angle_degree;     // angle of the polarizaiton vector versus X axis
    int type_laser; // CW, faconne , femto, pulsé , ...
    int coherent_avec_laser_num; // numéro du premier laser avec lequel il est cohérent -1 ou lui même si pas cohérent

public:
    map < double, double > spectre_Ecm_attenuation ; //La liste du spectre laser atténué
// C'est un tableau Energie_cm, I_attenuation (1= non atténué, 0= éteind).




public:                // Forme canonique d'une classe
    Laser();                        //constructeur
    Laser(const Laser&);            // Constructeur de copie
    Laser& operator = (const Laser&);       // Affectation par recopie
    virtual ~Laser();                       // Destructeur par defaut

    // Get components
    Vecteur3D  get_waist_pos()  const
    {
        return waist_pos;
    }
    Vecteur3D  get_direction()  const
    {
        return direction.unit(); // vecteur directeur normalisé
    }
    Vecteur3D  get_waist()  const
    {
        return waist;
    }
    double get_lambda()  const
    {
        return lambda;
    }
    double get_omega()  const // Pulsation
    {
        return 2*pi*C/lambda;
    }
    double get_Gamma_Laser()  const
    {
        return Gamma_Laser;
    }
    double get_Power()  const
    {
        return Power;
    }
    Vecteur3D get_polarisation()  const
    {
        return polarisation;
    }
    double get_polar_angle_degree()  const
    {
        return polar_angle_degree;
    }
    int get_type_laser()  const
    {
        return type_laser;
    }
    int get_coherent_avec_laser_num()  const
    {
        return coherent_avec_laser_num;
    }



    // Set components
    void  set_waist_pos(const Vecteur3D& new_pos)
    {
        waist_pos = new_pos;
    }
    void  set_direction(const Vecteur3D& new_dir)
    {
        direction = new_dir;
    }
    void  set_waist(const Vecteur3D& new_waist)
    {
        waist = new_waist;
    }
    void  set_lambda(const double& new_lambda)
    {
        lambda = new_lambda;
    }
    void  set_energy_cm(const double& new_energy_cm)
    {
        lambda = 1./(100.*new_energy_cm);
    }
    void  set_Gamma_Laser(const double& new_Gamma_Laser)
    {
        Gamma_Laser = new_Gamma_Laser;
    }
    void  set_Gamma_Laser_MHz(const double& new_Gamma_Laser_MHz)
    {
        Gamma_Laser = 2*pi*new_Gamma_Laser_MHz*1e6;
    }
    void  set_Power(const double& new_Power)
    {
        Power = new_Power;
    }
    void  set_polarisation(const Vecteur3D& new_pol)
    {
        polarisation = new_pol;
    }

    void  set_polar_angle(const double& new_polar_angle_degree)
    {
        polar_angle_degree = new_polar_angle_degree;
    }

    void set_type_laser(const int& new_type)
    {
        type_laser = new_type;
    }
    void set_coherent_avec_laser_num(const int& new_type)
    {
        coherent_avec_laser_num = new_type;
    }


// Lit les fichiers Energy_cm atténuation
// Si le fichier est inconnu cela ne rajoute rien à la liste
    void read_Spectrum(istream & flux);

    int read_Spectrum(const char *nom_file);

// Ecrit le spectre
    void write_Spectrum(ostream & flux);

    // Clear components
    void  clear_waist_pos()
    {
        waist_pos = Vecteur3D(0.,0.,0.);
    }
    void  clear_direction()
    {
        direction = Vecteur3D(0.,0.,1.);
    }
    void  clear_waist()
    {
        waist = Vecteur3D(1.,1.,0.);
    }
    void  clear_lambda()
    {
        lambda = 1.;
    }
    void  clear_Gamma_Laser()
    {
        Gamma_Laser = 1.;
    }
    void  clear_Power()
    {
        Power = 0.;
    }
    void  clear_polarisation()
    {
        polarisation =  Vecteur3D(0.,0.,0.);
    }
    void  clear_polar_angle()
    {
        polar_angle_degree = 0.;
    }

    void  clear_type_laser()
    {
        type_laser = CW;
    }
    void  clear_coherent_avec_laser_num()
    {
        coherent_avec_laser_num = -1;
    }



    // Comparison
    bool operator == (const Laser & Laser) const;
    bool operator != (const Laser & Laser) const;


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

// intensité au waist
    double  intensity()  const;

    // transmission (1 = 100%, 0 = 0) prenant en compte le spectre à l'énergie de la transition
    double transmission_spectrum(const double energie_trans_cm)  const;

// Zone de Rayleigh Vecteur3D
// ZRZ est la zone de Railieh moyenne
    Vecteur3D  Rayleigh_range()  const;

// waist au point (X,Y,Z)
    Vecteur3D  waist_size(const Vecteur3D& point)  const;

// intensité au point (X,Y,Z) dans le repère laser centré sur le waist
    double  intensity_repere_sur_waist(const Vecteur3D& point)  const;


// intensité au point (x,y,z)
// I.E. lié au repère du labo
    double  intensity_lab_axis(const Vecteur3D& point)  const;


// wave_vector k
// h c / \lambda = h c \sigma(m-1) = \hbar c k
    Vecteur3D wave_vector()  const;


    // Energie de la transition laser en cm^-1
    double Energy_transition_laser_cm()  const;


// ----------------------------------
//   FIN DES FONCTIONS EN LIGNES
// ----------------------------------

}
;



//----------------------------------
// AUTRES fonctions
//--------------------------------

// I = ε0 c E^2 /2.
double champ_E(const double irradiance);


// (absolute value of the)  effectif dipole d.e_laser = sum_p d_p epsilon^p
// where the dipole transition vector d= sum_p d_p e^p is given in the local quantification axis
// and the polarisation vector e_laser= sum_p' epsilon^p' e_p'  is given in the laser axis
double effectif_dipole_local(const Vecteur3D& dipole, const Vecteur3D& axe_quant,  const Laser& my_laser);



/***

Compare to the PRA 2014 we change notation now to be like Wikipedia in ZXZ convention(BE CAREFUL MATHEMATICA and Varshalovitch are in ZYZ convention)
 repère x,y,z du labo  et X,Y,Z de l'axe de quantification donné par le champ extérieur local
On utilise les angles d'Euler pour faire les rotations de repère
http://en.wikipedia.org/wiki/Euler_angles qui note (alpha,beta,gamma)


1. the first rotation is by an angle alpha (modulo 2pi, we chooose  (-pi,pi] to use atan2 ) about the z-axis,

2. the second rotation is by an angle beta in [0,pi] about the former (new) x-axis (now x')

3. the third rotation is by an angle gamma in (modulo 2pi, we chooose  (-pi,pi] to use atan2 ) about the former z-axis (now z').

***/

// Calcul des angles d'EULER. Pour un repère donné uniquement par son vecteur OZ=direction
// Il reste donc un arbitraire pour  choisir l'angle de rotation autour de cet axe pour son repère
// Nous choissons les angles tel que le repère soit le repère polaire dont OZ est la direction et OX selon le méridien
Vecteur3D  Euler_angles( Vecteur3D direction);

// Passage des coordonnées point(x,y,z) (labo) à (X,Y,Z): donné par les angles d'Euler alpha beta et gamma (qui font passer de e_x,e_y,e_z à e_X e_Y e_Z)
// cf http://mathworld.wolfram.com/EulerAngles.html
//Laser coordinates where point=(x,y,z) is lab coordinate
Vecteur3D  rotation_lab_axis(const Vecteur3D& point, double alpha, double beta, double gamma=0.);

// Passage des coordonnées point (X,Y,Z) donné dans le REPERE à labo (x,y,z) (repère labo).
// Le repère XYZ est donnée donné par les angles d'Euler alpha beta et gamma par rapport à xyz
Vecteur3D  rotation_axis_lab(const Vecteur3D& point, double alpha, double beta, double gamma=0.);

#endif
