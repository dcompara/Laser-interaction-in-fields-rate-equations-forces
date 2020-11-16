#include "Internal_state.h"

// ------------
// Constructors
// ------------

Internal_state::Internal_state()                    // Constructeur par défaut
{
<<<<<<< HEAD
    exc = 0;         // manifold
=======
    exc = 0;        // 0 est l'état fondamental
>>>>>>> f1d67ca6be17196db0b5ab5615163dc80d1182e6
    v =0;
    two_M = 0;
    two_J =0;
    bound_level  = 1;                  // Bound level 0=continuum, 1 = bound level
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
    v = Intern.v;
    two_M = Intern.two_M;
    two_J = Intern.two_J;
    bound_level  = Intern.bound_level ;
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
        v = Intern.v;
        two_M = Intern.two_M;
        two_J = Intern.two_J;
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
    }

    return *this;
}



// Affichage (on ne met pas la liste des raies!)
void Internal_state::write(ostream & flux)
{
    if (&flux == &cout)
        cout << "excitation : " << exc << "\t";
    else
        flux << exc << "\t";

    if (&flux == &cout)
        cout << "vibration (or variable) : " << v << "\t";
    else
        flux << v << "\t";

    if (&flux == &cout)
        cout << "J : " << two_J/2. << "\t";
    else
        flux << v << "\t";

    if (&flux == &cout)
        cout << "projection du moment total M  : " << two_M/2. << "\t";
    else
        flux << two_M/2. << "\t";

    if (&flux == &cout)
        cout << "bound level (0=no, 1=yes) : " << bound_level << "\t";
    else
        flux << bound_level << "\t";

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
            cout << "raie " << *(liste_raies[i].first) << " dipole_debye " << liste_raies[i].second << endl ;
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
        cout << "vibration (or variable): "  ;
    flux >> v ;

    if (&flux == &cin)
        cout << "projection du moment total M (2*) : "  ;
    flux >> two_M ;

     if (&flux == &cin)
        cout << "moment total J (2*) : "  ;
    flux >> two_J ;

    if (&flux == &cin)
        cout << "bound level (0=no, 1=yes) : "  ;
    flux >> bound_level ;

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
// On ne teste que l'état électronique, M, bound_level et numéro.
bool Internal_state::is_equal  (const Internal_state & w) const
{
    if( (exc == w.exc) && (two_M == w.two_M) && (bound_level  == w.bound_level) && (deg_number == w.deg_number) )
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
// bound_level
// #(number to discriminate the levels)
// Population
// v
// Energy0(in 0 field)
// Delta
// C


// Ecriture dans le flux du niveau avec le format souhaité
void Internal_state::write_Level_B(ostream & flux)
{
    flux << exc << "    ";
    flux << two_M << "    ";
    flux << bound_level  << "    ";
    flux << deg_number  << "    ";
    flux << population  << "    ";
    flux << v  << "    ";
    flux << two_J  << "    ";
    flux << Energy0_cm  << "    ";
    flux << Delta_FieldB  << "    ";
    flux << Linear_FieldB  << endl;
}


void Internal_state::read_Level__B(istream & flux)
{
    flux >> exc;
    flux >> two_M;
    flux >> bound_level ;
    flux >> deg_number ;
    flux >> population ;
    flux >> v ;
    flux >> two_J;
    flux >> Energy0_cm ;
    Energy_cm = Energy0_cm;
    flux >> Delta_FieldB ;
    flux >> Linear_FieldB ;
}

void Internal_state::read_Level__E(istream & flux)
{
    flux >> exc;
    flux >> two_M;
    flux >> bound_level ;
    flux >> deg_number ;
    flux >> population ;
    flux >> v ;
    flux >> two_J;
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
        double dip_Debye = liste_raies[i].second; // dipole en Debye de la transiton
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



// Lecture des fichierq Lignes qui contiennent
// UpperManifold	M'	bound_level'	#'	LowerManifold	M"	bound_level"	#"	Eupper	Elower	Dipole_Debye
void read_Line(istream & flux, Internal_state & Upper_state,Internal_state & Lower_state, double & Dipole_Debye)
{

    flux >> Upper_state.exc;
    flux >> Upper_state.two_M;
    flux >> Upper_state.bound_level ;
    flux >> Upper_state.deg_number ;

    flux >> Lower_state.exc;
    flux >> Lower_state.two_M;
    flux >> Lower_state.bound_level ;
    flux >> Lower_state.deg_number ;

    flux >> Upper_state.Energy0_cm;
    flux >>	Lower_state.Energy0_cm;
    flux >> Dipole_Debye;

    return;
}



// Ecriture dans le flux la transition (ligne) avec le format souhaité
void write_Line(ostream & flux, Internal_state & Upper_state, Internal_state & Lower_state, const double  Dipole_Debye)
{
    flux << Upper_state.exc << "    ";
    flux << Upper_state.two_M << "    ";
    flux << Upper_state.bound_level  << "    ";
    flux << Upper_state.deg_number  << "    ";

    flux << Lower_state.exc << "    ";
    flux << Lower_state.two_M << "    ";
    flux << Lower_state.bound_level  << "    ";
    flux << Lower_state.deg_number  << "    ";

    flux << Upper_state.Energy0_cm << "    ";
    flux <<	Lower_state.Energy0_cm << "    ";
    flux << Dipole_Debye << endl;

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


