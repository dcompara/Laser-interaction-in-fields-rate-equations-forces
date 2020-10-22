/*
  Name:  classe ��internal_state �
  Copyright:
  Author: Daniel Comparat
  Date: 15/10/06 11:01



Contient:
Level electronique: exc (excitation de la molecule, par exemple niveau �lectronique)
moment total J (ou F)
projection du moment total J (ou F) sur un axe suppos� connu: M
other variable (such as vibration, etc...): v
bound level:  +1  (or 0 for unbound states)
# = nombre pour lever la d�g�n�rescence
1/dur�e de vie (A Einstein du niveau)
shift en �nergie > 0 si le niveau initial gagne de l'�nergie (se rapproche des �tats excit� vers la ionisation)
population: La population (non normalis�e)

Coefficient pour un fit en champ F via formule: Energy (en cm^-1) + signe(Linear)*(-Delta/2+sqrt((Delta/2)^2+(Linear F)^2)). Ici l'unit� n'est pas sp�cifi�e mais souvent nous utiliserons comme unit� d'�nergie le cm^-1.
Liste des raies = la liste des transitions vers d'autres �tats

 * Tout initialis�e � 0. (sauf bound level et d�g�n�rescence +1)

 * write et read permettent d'�crire et de lire dans un flux (cout (pas cerr ou clog) ou fichier)
 * write_list_raie �crit la liste des raies


 unit�s DEBYE pour dipole et CM^-1 pour �nergie
* Levels
* Liste des raies = la liste des transitions vers d'autres �tats


La comparaison des niveaux se fait uniquement sur
Level electronique (Manifold), M, bound_level  et nombre pour lever la d�g�n�rescence

  */

/***
SOUVENT LES NOMBRES QUANTIQUES SONT EN FAIT *2
Ainsi M=1/2 est stock� en 1 d'ou le nom de  two_M.
***/


#ifndef Internal_state_SEEN
#define Internal_state_SEEN

#include "Vecteur3D.h"
#include "constantes_SI.h"
#include <vector>

#include <iostream>
using namespace std;



class Internal_state
{
public:
    int exc;                 // �tat excit� ou non (0 est l'�tat fondamental)
    int bound_level;                  // Parity
    int two_J;  // 2J
    int two_M;  // 2M
    int v;  // Other variable (such a vibrational number)
    int deg_number;           // nombre de d�g�n�rescence
    double Energy0_cm, Delta_FieldB, Linear_FieldB, Delta_FieldE, Linear_FieldE; // Energie en cm^-1. Coefficient pour un fit de l'�nergie en champ B ou E via formule: E0 + signe(C)*(-Delta/2+sqrt((Delta/2)^2+(C F)^2))
    double one_over_lifetime; // 1/(dur�e de vie en seconde)
    double Energy_cm, population; //  Energie r�elle due au shift en �nergie et population initiale de l'�tat

    typedef pair<Internal_state *,double> transition; // Une transition est vers un (pointeur vers) un �tat final (INTERNAL STATE) avec un dipole de transition (un double)

    vector < transition > liste_raies; //La liste des raies possibles est: vecteur donnant la liste des �tats finaux et de la force de transition

public:
    Internal_state();         // constructeurs autres que par recopie
    Internal_state(const Internal_state & Intern);   // constructeur de recopie
    ~Internal_state();           // destructeur
    Internal_state & operator = (const Internal_state & Intern);  // op�rateur d'affectation

// Affichage (on ne met pas la liste des raies!)
    void write(ostream & flux);

    // read des donn�es (sans la liste des raies)
    void read(istream & flux);

    // Ecrit la liste des raies
    void write_liste_raie(ostream & flux);


// True si l'�tat est = � w, false sinon
// On ne teste que l'�tat �lectronique, v, M, bound_level  et num�ro.
    bool is_equal (const Internal_state & w) const;


//----------------------------------
// Lecture des fichiers  de niveaux
//----------------------------------

// Les fichiers contiennent
// Manifold,  M, bound_level, #, J, N, Omega, Energy0, Delta, C


// Ecriture dans le flux du niveau avec le format souhait�
    void write_Level_B(ostream & flux);

// Si variation en B ou E
    void read_Level__B(istream & flux);
    void read_Level__E(istream & flux);

// Calcul du shift en champ: via E0 + signe(C)*(-Delta/2+sqrt((Delta/2)^2+(C F)^2))
    double   Energy_Shift_B_cm(const double B) const;
    double   Energy_Shift_E_cm(const double E) const;

// Calcul du gradient du shift (en unit� SI pas cm^-1) en champ
// Grad(F.F) *  signe(C)*/[2*sqrt((Delta/2)^2+(C F)^2))]
    Vecteur3D  Grad_Energy_Shift_B(const double B, const Vecteur3D gradient_B2);
    Vecteur3D  Grad_Energy_Shift_E(const double E, const Vecteur3D gradient_E2);

// Calcul du champ F correspondant au shift E_cm: via E0 + signe(C)*(-Delta/2+sqrt((Delta/2)^2+(C F)^2)) = E_cm
    double   Field_Shift_B_cm(const double Eshift);
    double   Field_Shift_E_cm(const double Eshift);

//----------------------------------
// Etudes des raies
//----------------------------------

// ajoute la transition -->v et la force de raie � la liste
    void add_transition (Internal_state  *v,const double strenght);

// Somme les taux de desexcitation
    double Einstein_desex_rate();


};



// Lecture des fichier  qui contiennent
// UpperManifold	M'	bound_level'	#'	LowerManifold	M"	bound_level"	#"	Eupper	Elower	Dipole_Debye
void read_Line(istream & flux, Internal_state & Upper_state,Internal_state & Lower_state,
                       double & Dipole_Debye);

// Ecriture dans le flux la transition (ligne) avec le format souhait�
void write_Line(ostream & flux, Internal_state & Upper_state, Internal_state & Lower_state, const double  Dipole_Debye);


//----------------------------------
// Surd�finition des entr�es sorties (sans la liste des raies)
//----------------------------------

ostream & operator << (ostream & flux, Internal_state Intern);

istream& operator >> (istream &flux, Internal_state & Intern);


#endif
