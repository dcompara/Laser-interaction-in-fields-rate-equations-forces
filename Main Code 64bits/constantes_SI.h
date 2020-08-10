/*
  Name: Constantes physiques du système d'unité international
  Copyright:
  Author: Daniel COMPARAT
  Date: 23/10/06 09:10
  Description:
      Basé sur  CLHEP - a Class Library for High Energy Physics.
      http://wwwinfo.cern.ch/asd/lhc++/clhep/index.html

// HEP coherent Physical Constants
// Below is a non exhaustive list of Physical CONSTANTS,
//
// The basic units are :
//      meter  (not millimeter as in CLHEP)
//   second (not nanosecond )
//   electron Volt  (not Mev)
//   degree Kelvin
//              amount of substance (mole)
//              luminous intensity (candela)
//   radian
//              steradian

*/


#ifndef constantes_SI_SEEN
#define constantes_SI_SEEN


 #include <math.h>
 #include "Vecteur3D.h"

using namespace std;


const double  VERY_LARGE_NUMBER = 1.e90;
const double VERY_SMALL_NUMBER = 1./VERY_LARGE_NUMBER;
const double SMALL_NUMBER = 1.e-10; // Utile pour les différences finies (pour éviter les erreurs d'arrondis)
const double  LARGE_NUMBER = 1./SMALL_NUMBER;
const double SMALL_NUMBER_RATE = 1.; // Rate smaller than this (1 second) will not be considered
const double SMALL_DIPOLE_DEBYE = 1.e-8; // Rate smaller than this (1 second) will not be considered


const int aucune = -1; // numéro pour indiquer aucune molécule
const int all = -1; // numéro pour indiquer toute les molécule

static const double     pi  = 3.14159265358979323846;
static const double  twopi  = 2*pi;
static const double halfpi  = pi/2;
static const double     pi2 = pi*pi;
static const double SQRT2 = sqrt(2.);

static const double g_grav = 9.80665; // Gravité standard
static const Vecteur3D gravity(0.,0.,-g_grav);

static const double EAU = 219474.63137032; // Energy atomic units -> cm^(-1)
static const double hartree = EAU; // Energy atomic units -> cm^(-1)
static const double hartreeJ = 4.35974417*1e-18; // Hartree en J
static const double BAU = 2.35051809*1e5;
static const double ME = 9.1093697*1e-31;
static const double FAU = 5.1422082*1e11;
static const double ASO = 2*554.0406*hartree/3;
static const double MAU = 1.6605402*1e-27;
static const double Mproton = 1.67262178*1e-27;
static const double mau = 1.6605402*1e-27;
static const double MCs = 132.905442*MAU;
static const double MRb87 = 86.90918052*MAU;
static const double MRb85 = 84.911789732*MAU;
static const double MNa = 22.98976928*MAU;
static const double MH = 1.00794*MAU;
static const double MCs2 = 2.*132.905442*MAU;
static const double MCO = 28.*MAU;
static const double MNH = 15.0146*MAU;
static const double MBaF = 156.325*MAU;
static const double MLi7Cs = (7+133)*MAU;
static const double MLi6Cs = (6+133)*MAU;
static const double MRb85Cs = (85+133)*MAU;
static const double MRb87Cs = (87+133)*MAU;
static const double MPs = 2.*ME;
static const double MRb = 87.*mau;
static const double MC2moins = 24.*MAU;
static const double kCoulomb = 8.9875517873681764*1e9;  // 1/(4 pi epsilon_0)

static const double H = 6.6260755*1e-34;
static const double hPlanck = 6.6260755*1e-34;
static const double HBAR = H/(2*pi);
static const double C = 299792458;
static const double QE = -1.60217733*1e-19;
static const double QION = -QE;
static const double A0 = 0.529177249*1e-10;
static const double a0 = 0.529177249*1e-10;
static const double MU0 = 4*pi*1e-7;
static const double EPSILON0 = 1/(MU0*C*C);
static const double e2 = QE*QE*C*C*1e-7;  // e2 = q_e^2/(4 Pi Epsilon0)
static const double ALPHA = e2/(HBAR*C);
static const double MUBOHR = -QE*HBAR/(2*ME); // ATTENTION POSITIF!!
static const double MHz = 1e6 ;
static const double MHzcm = C/10000 ;//  cm^(-1) -> MHz
static const double VparcmMHz = -QE*ALPHA*C/hartreeJ/1e6; // champ électrique V/cm -> MHz
static const double CenCm =QE*A0*A0/(ALPHA*C); // MHz -> C*m
static const double kB = 1.380658*1e-23;
static const double MUKMHZ = 1e-6*kB/(1e6*H);
static const double mW = 1e-3;
static const double nW = 1e-9;
static const double MICRON = 1e-6;
static const double Mus = MICRON; // microseconde
// static const double conversionUAEnergieMHz=pow(QE*A0,2)/(H*pow(MICRON,3)*MHz);// utilisé dans les calculs de potentiels
static const double Debye = 1.e-21/C;
static const double Spol_Debye_A_s = (8e6*pi*pi*C*C*C*Debye*Debye)/(3.*EPSILON0*C*C*C*HBAR); // 3.13618932*10^-7 = Conversion HonlLondon (en Debye) en A Einstein (s^-1) si energie en cm^-1
static const double Conv_Ecm_delta = 2.*pi*MHz*MHzcm ;//  cm^(-1) -> s^-1;
static const double  sigmaSB = 5.670367e-8 ; // Stefan–Boltzmann



// Conversion entre l'unité des énergie dans les fichiers et celle des calcul de taux

//
// Length [L]
//
static const double millimeter  = 0.001;
static const double millimeter2 = millimeter*millimeter;
static const double millimeter3 = millimeter*millimeter*millimeter;

static const double centimeter  = 10.*millimeter;
static const double centimeter2 = centimeter*centimeter;
static const double centimeter3 = centimeter*centimeter*centimeter;

static const double meter  = 1000.*millimeter;
static const double meter2 = meter*meter;
static const double meter3 = meter*meter*meter;

static const double kilometer = 1000.*meter;
static const double kilometer2 = kilometer*kilometer;
static const double kilometer3 = kilometer*kilometer*kilometer;

static const double parsec = 3.0856775807e+16*meter;

static const double micrometer = 1.e-6 *meter;
static const double  nanometer = 1.e-9 *meter;
static const double  angstrom  = 1.e-10*meter;
static const double  fermi     = 1.e-15*meter;

static const double      barn = 1.e-28*meter2;
static const double millibarn = 1.e-3 *barn;
static const double microbarn = 1.e-6 *barn;
static const double  nanobarn = 1.e-9 *barn;
static const double  picobarn = 1.e-12*barn;

// symbols
static const double nm  = nanometer;
static const double um  = micrometer;

static const double mm  = millimeter;
static const double mm2 = millimeter2;
static const double mm3 = millimeter3;

static const double cm  = centimeter;
static const double cm2 = centimeter2;
static const double cm3 = centimeter3;

static const double m  = meter;
static const double m2 = meter2;
static const double m3 = meter3;

static const double km  = kilometer;
static const double km2 = kilometer2;
static const double km3 = kilometer3;

static const double pc = parsec;

//
// Angle
//
static const double radian      = 1.;
static const double milliradian = 1.e-3*radian;
static const double degree = (3.14159265358979323846/180.0)*radian;

static const double   steradian = 1.;

// symbols
static const double rad  = radian;
static const double mrad = milliradian;
static const double sr   = steradian;
static const double deg  = degree;

//
// Time [T]
//
static const double nanosecond  = 1.e-9;
static const double second      = 1.e+9 *nanosecond;
static const double millisecond = 1.e-3 *second;
static const double microsecond = 1.e-6 *second;
static const double  picosecond = 1.e-12*second;

static const double hertz = 1./second;
static const double kilohertz = 1.e+3*hertz;
static const double megahertz = 1.e+6*hertz;

// symbols
static const double ns = nanosecond;
static const double  s = second;
static const double ms = millisecond;

//
// Electric charge [Q]
//
// static const double eplus = 1. ;// positron charge
static const double eplus = 1.60217733e-19 ;// positron charge
static const double e_SI  = 1.60217733e-19;// positron charge in coulomb
static const double coulomb = eplus/e_SI;// coulomb = 6.24150 e+18 * eplus



//
// Energy [E]
//
//static const double megaelectronvolt = 1.e6 ; // CERN value
static const double megaelectronvolt = 1.e6*e_SI ; // SI_value
static const double     electronvolt = 1.e-6*megaelectronvolt;
static const double kiloelectronvolt = 1.e-3*megaelectronvolt;
static const double gigaelectronvolt = 1.e+3*megaelectronvolt;
static const double teraelectronvolt = 1.e+6*megaelectronvolt;
//static const double petaelectronvolt = 1.e+9*megaelectronvolt;
//

// Joule ici 1 !
static const double joule = electronvolt/e_SI;// joule = 6.24150 e+12 * MeV




// symbols
static const double MeV = megaelectronvolt;
static const double  eV = electronvolt;
static const double keV = kiloelectronvolt;
static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;

//
// Mass [E][T^2][L^-2]
//
static const double  kilogram = joule*second*second/(meter*meter);
static const double      gram = 1.e-3*kilogram;
static const double milligram = 1.e-3*gram;

// symbols
static const double  kg = kilogram;
static const double   g = gram;
static const double  mg = milligram;

//
// Power [E][T^-1]
//
static const double watt = joule/second;// watt = 6.24150 e+3 * MeV/ns

//
// Force [E][L^-1]
//
static const double newton = joule/meter;// newton = 6.24150 e+9 * MeV/mm

//
// Pressure [E][L^-3]
//

static const double hep_pascal = newton/m2;   // pascal = 6.24150 e+3 * MeV/mm3
static const double bar        = 100000*hep_pascal; // bar    = 6.24150 e+8 * MeV/mm3
static const double atmosphere = 101325*hep_pascal; // atm    = 6.32420 e+8 * MeV/mm3

//
// Electric current [Q][T^-1]
//
static const double      ampere = coulomb/second; // ampere = 6.24150 e+9 * eplus/ns
static const double milliampere = 1.e-3*ampere;
static const double microampere = 1.e-6*ampere;
static const double  nanoampere = 1.e-9*ampere;

//
// Electric potential [E][Q^-1]
//
//static const double megavolt = megaelectronvolt/eplus;
static const double megavolt = 1.e6;
static const double kilovolt = 1.e-3*megavolt;
static const double     volt = 1.e-6*megavolt;

//
// Electric resistance [E][T][Q^-2]
//
static const double ohm = volt/ampere;// ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

//
// Electric capacitance [Q^2][E^-1]
//
static const double farad = coulomb/volt;// farad = 6.24150e+24 * eplus/Megavolt
static const double millifarad = 1.e-3*farad;
static const double microfarad = 1.e-6*farad;
static const double  nanofarad = 1.e-9*farad;
static const double  picofarad = 1.e-12*farad;

//
// Magnetic Flux [T][E][Q^-1]
//
static const double weber = volt*second;// weber = 1000*megavolt*ns

//
// Magnetic Field [T][E][Q^-1][L^-2]
//
static const double tesla     = volt*second/meter2;// tesla =0.001*megavolt*ns/mm2
static const double milli_tesla     = 0.001;// tesla =0.001*megavolt*ns/mm2

static const double gauss     = 1.e-4*tesla;
static const double kilogauss = 1.e-1*tesla;

//
// Inductance [T^2][E][Q^-2]
//
static const double henry = weber/ampere;// henry = 1.60217e-7*MeV*(ns/eplus)**2

//
// Temperature
//
static const double kelvin = 1.;
static const double  mK = 1e-3*kelvin;
static const double  muK = 1e-6*kelvin;

//
// Amount of substance
//
static const double mole = 1.;

//
// Activity [T^-1]
//
static const double becquerel = 1./second ;
static const double curie = 3.7e+10 * becquerel;

//
// Absorbed dose [L^2][T^-2]
//
static const double      gray = joule/kilogram ;
static const double  kilogray = 1.e+3*gray;
static const double milligray = 1.e-3*gray;
static const double microgray = 1.e-6*gray;

//
// Luminous intensity [I]
//
static const double candela = 1.;

//
// Luminous flux [I]
//
static const double lumen = candela*steradian;

//
// Illuminance [I][L^-2]
//
static const double lux = lumen/meter2;

//
// Miscellaneous
//
static const double perCent     = 0.01 ;
static const double perThousand = 0.001;
static const double perMillion  = 0.000001;



//
//
//
static const double Avogadro = 6.0221367e+23/mole;

//
// c   = 299.792458 mm/ns
// c^2 = 898.7404 (mm/ns)^2
//
static const double c_light   = 2.99792458e+8 * m/s;
static const double c_squared = c_light * c_light;

//
// h     = 4.13566e-12 MeV*ns
// hbar  = 6.58212e-13 MeV*ns
// hbarc = 197.32705e-12 MeV*mm
//
static const double h_Planck      = 6.6260755e-34 * joule*s;
static const double hbar_Planck   = h_Planck/twopi;
static const double hbarc         = hbar_Planck * c_light;
static const double hbarc_squared = hbarc * hbarc;

//
//
//
static const double electron_charge = - eplus; // see SystemOfUnits.h
static const double e_squared = eplus * eplus;

//
// amu_c2 - atomic equivalent mass unit
// amu    - atomic mass unit
//
static const double electron_mass_c2 = 0.51099906 * MeV;
static const double   proton_mass_c2 = 938.27231 * MeV;
static const double  neutron_mass_c2 = 939.56563 * MeV;
static const double           amu_c2 = 931.49432 * MeV;
static const double              amu = amu_c2/c_squared;

//
// permeability of free space mu0    = 2.01334e-16 Mev*(ns*eplus)^2/mm
// permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV*mm)
//
static const double mu0      = 4*pi*1.e-7 * henry/m;
static const double epsilon0 = 1./(c_squared*mu0);

//
// electromagnetic coupling = 1.43996e-12 MeV*mm/(eplus^2)
//
static const double elm_coupling           = e_squared/(4*pi*epsilon0);
static const double fine_structure_const   = elm_coupling/hbarc;
static const double classic_electr_radius  = elm_coupling/electron_mass_c2;
static const double electron_Compton_length = hbarc/electron_mass_c2;
static const double Bohr_radius = electron_Compton_length/fine_structure_const;


//
// static const double k_Boltzmann = 8.617385e-11 * MeV/kelvin;
static const double   k_Boltzmann = 1.3806504e-23;

#endif
