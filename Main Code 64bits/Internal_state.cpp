#include "Internal_state.h"

// ------------
// Constructors
// ------------

Internal_state::Internal_state()                    // Constructeur par défaut
{
    exc = 0;        // 0 est l'état fondamental
    two_Spin = two_Lambda = two_Omega = 0;
    v =0;
    two_J = two_N = two_M = 0;
    Sym = 1;                  // Parity
    deg_number = 0;           // nombre de dégénérescence
    Energy0_cm = Delta_FieldB = Linear_FieldB = Delta_FieldE = Linear_FieldE = 0.;
    one_over_lifetime=0.;
    Energy_cm = population = 0.;
    liste_raies = vector < transition > ();
};

Internal_state::~Internal_state()
{}
;                    // Destructeur


Internal_state::Internal_state(const Internal_state & Intern)  // Constructeur de (re)copie
{
    exc = Intern.exc;
    two_Spin = Intern.two_Spin;
    two_Lambda = Intern.two_Lambda;
    two_Omega = Intern.two_Omega;
    v = Intern.v;
    two_J = Intern.two_J;
    two_N = Intern.two_N;
    two_M = Intern.two_M;
    Sym = Intern.Sym;
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

};


Internal_state & Internal_state::operator = (const Internal_state & Intern)  // opérateur d'affectation
{
    if (this != &Intern) // On vérifie que les objets ne sont pas les mêmes !
    {
        exc = Intern.exc;
        two_Spin = Intern.two_Spin;
        two_Lambda = Intern.two_Lambda;
        two_Omega = Intern.two_Omega;
        v = Intern.v;
        two_J = Intern.two_J;
        two_N = Intern.two_N;
        two_M = Intern.two_M;
        Sym = Intern.Sym;
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
    }

    return *this;
}

//  void  set_int(const Internal_state & Intern)  // opérateur d'affectation
//  {
//    exc = Intern.exc;
//    v = Intern.v;
//    J = Intern.J;
//    MJ = Intern.MJ;
//  }



// Affichage (on ne met pas la liste des raies!)
void Internal_state::write(ostream & flux)
{
    if (&flux == &cout)
        cout << "excitation : " << exc << "\t";
    else
        flux << exc << "\t";

    if (&flux == &cout)
        cout << "Spin  : " << two_Spin/2. << "\t";
    else
        flux << two_Spin/2. << "\t";

    if (&flux == &cout)
        cout << "Lambda : " << two_Lambda/2. << "\t";
    else
        flux << two_Lambda/2. << "\t";

    if (&flux == &cout)
        cout << "Omega : " << two_Omega/2. << "\t";
    else
        flux << two_Omega/2. << "\t";

    if (&flux == &cout)
        cout << "vibration : " << v << "\t";
    else
        flux << v << "\t";

    if (&flux == &cout)
        cout << "rotation J  : " << two_J/2. << "\t";
    else
        flux << two_J/2. << "\t";

    if (&flux == &cout)
        cout << "rotation N  : " << two_N/2. << "\t";
    else
        flux << two_N/2. << "\t";

    if (&flux == &cout)
        cout << "projection du moment total M  : " << two_M/2. << "\t";
    else
        flux << two_M/2. << "\t";

    if (&flux == &cout)
        cout << "Parité (+/- 1) : " << Sym << "\t";
    else
        flux << Sym << "\t";

    if (&flux == &cout)
        cout << "Numéro pour lever la dégénérescence " << deg_number << "\t";
    else
        flux << deg_number << "\t";

    if (&flux == &cout)
        cout << "Energy (cm^-1) " << Energy0_cm << "\t";
    else
        flux << Energy0_cm << "\t";

    if (&flux == &cout)
        cout << "Delta Field B " << Delta_FieldB << "\t";
    else
        flux << Delta_FieldB << "\t";

    if (&flux == &cout)
        cout << "Linear Field B " << Linear_FieldB << "\t";
    else
        flux << Linear_FieldB << "\t";

    if (&flux == &cout)
        cout << "Delta Field E " << Delta_FieldE << "\t";
    else
        flux << Delta_FieldE << "\t";

    if (&flux == &cout)
        cout << "Linear Field E " << Linear_FieldE << "\t";
    else
        flux << Linear_FieldE << "\t";

    if (&flux == &cout)
        cout << "1/Lifetime " << one_over_lifetime << "\t";
    else
        flux << one_over_lifetime << "\t";

    if (&flux == &cout)
        cout << "Energie shiftée (cm-1) " << Energy_cm << "\t";
    else
        flux << Energy_cm << "\t";

    if (&flux == &cout)
        cout << "Population " << population << "\t";
    else
        flux << population << "\t";

}


// Ecrit la liste des raies
void Internal_state::write_liste_raie(ostream & flux)
{
    if ((int) liste_raies.size() ==0)
    {
        cout << " VIDE " << endl;
        flux << " VIDE " << endl;
        return;
    }

    if (&flux == &cout)
        for( int i = 0; i < (int) liste_raies.size(); i++ )
            cout << "raie " << *(liste_raies[i].first) << " Force " << liste_raies[i].second << endl ;
    else
        for( int i = 0; i < (int) liste_raies.size(); i++ )
            flux << *(liste_raies[i].first) <<  liste_raies[i].second << endl;
}




// read des données (sans la liste des raies)
void Internal_state::read(istream & flux)
{
    if (&flux == &cin)
        cout << "entrez excitation : ";
    flux >> exc;

    if (&flux == &cin)
        cout << "Spin (2*spin) : ";
    flux >> two_Spin ;

    if (&flux == &cin)
        cout << "Lambda (2*) : "  ;
    flux >> two_Lambda ;

    if (&flux == &cin)
        cout << "Omega (2*) : " ;
    flux >> two_Omega ;


    if (&flux == &cin)
        cout << "vibration : "  ;
    flux >> v ;

    if (&flux == &cin)
        cout << "rotation J (2*) : "  ;
    flux >> two_J ;

    if (&flux == &cin)
        cout << "rotation N (2*) : "  ;
    flux >> two_N ;

    if (&flux == &cin)
        cout << "projection du moment total M (2*) : "  ;
    flux >> two_M ;

    if (&flux == &cin)
        cout << "Parité (+/- 1) : "  ;
    flux >> Sym ;

    if (&flux == &cin)
        cout << "Numéro pour lever la dégénérescence "  ;
    flux >> deg_number ;

    if (&flux == &cin)
        cout << "Energy0_cm (cm^-1) "  ;
    flux >> Energy0_cm ;

    if (&flux == &cin)
        cout << "Delta Field B "  ;
    flux >> Delta_FieldB ;

    if (&flux == &cin)
        cout << "Linear Field B "  ;
    flux >> Linear_FieldB ;

    if (&flux == &cin)
        cout << "Delta Field E "  ;
    flux >> Delta_FieldE ;

    if (&flux == &cin)
        cout << "Linear Field E "  ;
    flux >> Linear_FieldE ;

    if (&flux == &cin)
        cout << "1/Lifetime "  ;
    flux >> one_over_lifetime ;

    if (&flux == &cin)
        cout << "Energie shiftée "  ;
    flux >> Energy_cm;

    if (&flux == &cin)
        cout << "Population "  ;
    flux >>  population;

}


// True si l'état est = à w, false sinon
// On ne teste que l'état électronique, M, parité et numéro.
bool Internal_state::is_equal  (const Internal_state & w) const
{
    if( (exc == w.exc) && (two_M == w.two_M) && (Sym == w.Sym) && (deg_number == w.deg_number) )
    {
        return true;
    }

    return false;
}



//----------------------------------
// Lecture des fichiers en Energie
//----------------------------------

// Read Energy List
// THE FILE CONTAINS:
// Manifold (1 for upper or 0 for lower typically)
// M (in Field it is real M)
// Sym(Parity +/-)
// #(number to discriminate the levels)
// Population
// v
// J
// N
// Omega
// Energy0(in 0 field)
// Delta
// C


// Ecriture dans le flux du niveau avec le format souhaité
void Internal_state::write_Level_PgopherB(ostream & flux)
{
    flux << exc << "    ";
    flux << two_M << "    ";
    flux << Sym  << "    ";
    flux << deg_number  << "    ";
    flux << population  << "    ";
    flux << v  << "    ";
    flux << two_J << "    ";
    flux << two_N << "    ";
    flux << two_Omega << "    ";
    flux << Energy0_cm  << "    ";
    flux << Delta_FieldB  << "    ";
    flux << Linear_FieldB  << endl;
}


void Internal_state::read_Level_Pgopher_B(istream & flux)
{
    flux >> exc;
    flux >> two_M;
    flux >> Sym ;
    flux >> deg_number ;
    flux >> population ;
    flux >> v ;
    flux >> two_J;
    flux >> two_N;
    flux >> two_Omega;
    flux >> Energy0_cm ;
    Energy_cm = Energy0_cm;
    flux >> Delta_FieldB ;
    flux >> Linear_FieldB ;
}

void Internal_state::read_Level_Pgopher_E(istream & flux)
{
    flux >> exc;
    flux >> two_M;
    flux >> Sym ;
    flux >> deg_number ;
    flux >> population ;
    flux >> v ;
    flux >> two_J;
    flux >> two_N;
    flux >> two_Omega;
    flux >> Energy0_cm ;
    Energy_cm = Energy0_cm;
    flux >> Delta_FieldE ;
    flux >> Linear_FieldE ;
}


// Signe d'une variable -1 si négatif et +1 si positif et 0 si nul
int signe(const double x)
{
    if (x < 0.) return -1;
    if (x > 0.) return 1;
    return 0 ;
}


// Calcul du shift en champ
//via E0 + signe(C)*(-Delta/2+sqrt((Delta/2)^2+(C F)^2))
double  Internal_state::Energy_Shift_B_cm(const double B) const
{
    return signe(Linear_FieldB)*(-Delta_FieldB/2.+sqrt((Delta_FieldB/2.)*(Delta_FieldB/2.)+(Linear_FieldB * B)*(Linear_FieldB * B)));
}

double   Internal_state::Energy_Shift_E_cm(const double E) const
{
    return signe(Linear_FieldE)*(-Delta_FieldE/2.+sqrt((Delta_FieldE/2.)*(Delta_FieldE/2.)+(Linear_FieldE * E)*(Linear_FieldE * E)));
}



// Calcul du gradient du shift (en unité SI pas cm^-1) en champ
// Grad(F.F) *  C*C*signe(C)*/[2*sqrt((Delta/2)^2+(C F)^2))]
Vecteur3D  Internal_state::Grad_Energy_Shift_B(const double B, const Vecteur3D gradient_B2)
{
    if (abs(Linear_FieldB) < VERY_SMALL_NUMBER) return Vecteur3D(0.,0.,0.); // Pour éviter des divisions par zéro
    return HBAR*Conv_Ecm_delta*0.5*gradient_B2*signe(Linear_FieldB)*Linear_FieldB*Linear_FieldB/sqrt((Delta_FieldB/2.)*(Delta_FieldB/2.)+(Linear_FieldB * B)*(Linear_FieldB * B));
}

Vecteur3D Internal_state::Grad_Energy_Shift_E(const double E, const Vecteur3D gradient_E2)
{
    if (abs(Linear_FieldE) < VERY_SMALL_NUMBER) return Vecteur3D(0.,0.,0.); // Pour éviter des divisions par zéro
    return HBAR*Conv_Ecm_delta*0.5*gradient_E2*signe(Linear_FieldE)*Linear_FieldE*Linear_FieldE/sqrt((Delta_FieldE/2.)*(Delta_FieldE/2.)+(Linear_FieldE * E)*(Linear_FieldE * E));
}



// Calcul du champ F correspondant au shift E_cm: via E0 + signe(C)*(-Delta/2+sqrt((Delta/2)^2+(C F)^2)) = E_cm = Eshift
double   Internal_state::Field_Shift_B_cm(const double Eshift)
{
    if (abs(Linear_FieldB) < VERY_SMALL_NUMBER) return 0.; // Pour éviter des divisions par zéro
    return sqrt(Eshift*(Eshift+ Delta_FieldB*signe(Linear_FieldB)))/(signe(Linear_FieldB)*Linear_FieldB);
}



double   Internal_state::Field_Shift_E_cm(const double Eshift)
{
    if (abs(Linear_FieldE) < VERY_SMALL_NUMBER) return 0.; // Pour éviter des divisions par zéro
    return sqrt(Eshift*(Eshift+ Delta_FieldE*signe(Linear_FieldE)))/(signe(Linear_FieldE)*Linear_FieldE);
}



// Somme les taux de desexcitation
double Internal_state::Einstein_desex_rate()
{
    double Gamma = VERY_SMALL_NUMBER; // Not 0. To avoid division by zero

    for( int i = 0; i < (int) liste_raies.size(); i++ )
    {
        double dip_Debye = sqrt(3.*liste_raies[i].second); // dipole en Debye de la transiton
        double energie_cm =  Energy_cm - (liste_raies[i].first)->Energy_cm ; // Energie en cm^-1 de la transition
        if (energie_cm > 0.) // Transition  en émission
            Gamma += dip_Debye*dip_Debye*Spol_Debye_A_s*energie_cm*energie_cm*energie_cm;
        //  Spol_Debye_A_s = (8e6*pi*pi*C*C*C*Debye*Debye)/(3.*EPSILON0*C*C*C*HBAR); // 3.13618932*10^-7 = Conversion HonlLondon (en Debye) en A Einstein (s^-1) si energie en cm^-1
    }
    return Gamma;
}


//----------------------------------
// Lecture des fichiers  en Raies
//----------------------------------


void Internal_state::add_transition (Internal_state *v,const double strenght)
{
    transition raie;
    raie = make_pair (v,strenght);
    liste_raies.push_back (raie);
    return;
}


/*
 PGOPHER GIVES:

Molecule    responsible for this transition.
M'          Manifold for the upper state of the transition
J'          Upper state total angular momentum or M'. (M' in the presence of static field, otherwise J' in the absence of hyperfine structure, otherwise F').
S'          Symmetry of the upper state.
#'          Eigenvalue number for upper state. This gives the index (starting from 1) of the upper energy level with respect to other levels of the same total angular momentum and symmetry.
M''         for the lower state of the transition
J''         Lower state total angular momentum or M". (M" in the presence of static field, otherwise J' in the absence of hyperfine structure, otherwise F').
S''         Symmetry of the lower state.
#''         Eigenvalue number for lower state. This gives the index (starting from 1) of the upper energy level with respect to other levels of the same total angular momentum and symmetry.
Position    The position of the transition; the units for this are the units used in the main simulation window. For line position fitting, alter this number to the position of the observed peak, either manually or (more usually) with the assignment process described under line position fits.
Std Dev     The (relative) standard deviation of this transition, used to calculate the weight of the observation in a line position fit. If this is negative, zero or blank this transition will not be included in the fit. Entries will initially have this blank or negative, and the automatic assignment process will set this to 1. Less certain measurements can be given larger values. Note that the parameters produced by a fit only depend on the relative values of the weights. If UNREGISTERED EVALUATION VERSION is set at the top level, then this column is set to the negative of the estimated uncertaInternal_statey of the transition, based on the errors in the most recent fit, if available.
strengths    The Internal_stateensity of the rotational transition.
Width       The width of the transition. (Note that this does not include the Gaussian or Lorentzian width set on the plot window.)
deltaJ(J'') The transition strength in branch notation, such as P(1). The display will depend on the molecule type and the electronic and nuclear angular momenta included in the calculation.
Name        Details of the transition, including other quantum numbers where appropriate. The assignment process will also add the source of the measurement (typically a filename) where possible.

Name contains:
    The manifold and state name
    J: The J quantum number; not shown if ShowJ is false at the Molecule level
    N: The N quantum number; not shown if ShowN is false at the Molecule level or all states are singlet states
    Ω: The Ω quantum number; not shown if ShowOmega is false at the Molecule level (the default) or all states are singlet states
    Fn: The component of the multiplet numbered from 1 in order of increasing energy; not shown if ShowFNumber is false at the Molecule level or all states are singlet states.
                This contains the same information as the Ω quantum number, so it does not usually make sense to show both.
    e/f: The parity; not shown if Showef is false at the Molecule level.
    Hyperfine: quantum numbers are added at the end as required.

 BUT HERE WE USE ONLY:
 UpperManifold	M'	Sym'	#'	LowerManifold	M"	Sym"	#"	Position	Intensity	Eupper	Elower	Spol

REMEMBER That we use Pgopher with the option 2J thus tjhe quantum numbers are *2: two_J
*/


void read_Line_Pgopher(istream & flux, Internal_state & Upper_state,Internal_state & Lower_state, double & Spol, double & position, double  & Intensity)
{

    flux >> Upper_state.exc;
    flux >> Upper_state.two_M;
    flux >> Upper_state.Sym ;
    flux >> Upper_state.deg_number ;

    flux >> Lower_state.exc;
    flux >> Lower_state.two_M;
    flux >> Lower_state.Sym ;
    flux >> Lower_state.deg_number ;

    flux >> position;
    flux >> Intensity;

    flux >> Upper_state.Energy0_cm;
    flux >>	Lower_state.Energy0_cm;
    flux >> Spol;

    return;
}



// Ecriture dans le flux la transition (ligne) avec le format souhaité
void write_Line_Pgopher(ostream & flux, Internal_state & Upper_state, Internal_state & Lower_state, const double  Spol, const double position, const double Intensity)
{
    flux << Upper_state.exc << "    ";
    flux << Upper_state.two_M << "    ";
    flux << Upper_state.Sym  << "    ";
    flux << Upper_state.deg_number  << "    ";

    flux << Lower_state.exc << "    ";
    flux << Lower_state.two_M << "    ";
    flux << Lower_state.Sym  << "    ";
    flux << Lower_state.deg_number  << "    ";

    flux << position << "    ";
    flux << Intensity << "    ";

    flux << Upper_state.Energy0_cm << "    ";
    flux <<	Lower_state.Energy0_cm << "    ";
    flux << Spol << endl;

    return;
}


//----------------------------------
// Surdéfinition des entrées sorties  (sans la liste des raies)
//----------------------------------

ostream & operator << (ostream & flux, Internal_state Intern)
{
    Intern.write(flux);
    return(flux);
}

istream& operator >> (istream &flux, Internal_state & Intern) //Intern est modifié!
{
    Intern.read(flux);
    return(flux);
}


