/*

Lecture de
Facteurs de Franck-Condon
Energies des niveaux de vibration (la rotation est ensuite calculée par le Bv)
Calculs de facteurs angulaires liés à la rotation


Il y a FC_abs pour l'absorption et FC_em pour l'emission

On utilise de façon générique le nom
X pour l'état fond
A pour l'état excité
donc vX, vA pour les niveaux vibrationnels et
JX, JA pour la rotation




Il y a donc des tableaux
FC_abs[vA][vX]  // Qui contient un Franck Condon pour l'absorption (Normalisés à 1 pour avoir le bon taux de transition)
proba_abs[JA][MJA][JX][MJX] et proba_em[[JA][MJA][JX][MJX] // Qui contient les facteurs genre Hönl London (Normalisés à 1 pour avoir le bon taux de transition)
// En fait pour eviter des indices <0 on note [JA][MJA+JA][JX][MJX+JX]

FC_em[vA][vX] // Qui contient un Franck Condon pondéré par un terme en (omega/omega_moyen)^3 pour tenir compte de la durée de vie

EvX[vX][JX] // Energie des états vX, JX, Pour l'instant très simple avec Bv seul
EvA[vA][JA] // Energie des états vA, JA


Les bordures sont absorbantes
i.e. les états  de vmax ou Jmax ont des Franck Condon de zéro et restent ainsi piégés


LES UNITES ICI SONT CELLES DES FICHIERS (cm^-1 typiquement)

ATTENTION
La notation  [NA][NX] indique qu'il y a
NA lignes et NX colonnes dans les tableaux
Il faut donc bien vérifier que dans les fichiers
vA indice les lignes et vX indice les colonnes

En effet les boucles sont vA puis vX


*/

#include "Transitions_initialisation.h"


/*** Program to read Line or Energy list from Pgopher ***/

// Author: Daniel Comparat
//  Date: 12/02/2012



/** ENERGY LIST **/
// PGOPHER Gives
// Molecule Manifold   M  Sym  #  g Population Label  Energy Linear   Dipole  Err  Quadratic       Err  Two_Level Energy Delta     C              Dipole2    Err
// WE DO NOT USE IT LIKE THAT BUT USING ORIGIN WE MODIFY AND USE ONLY
// Manifold(1 for upper or 0 for lower typically)  M(in Field it is real M)  Sym(Parity +/-)  #(number to discriminate the levels) Population J N  Energy0(in 0 field)   Delta   C
// Where for fit in Energy in Fields  Energy +/- Sqrt(Delta^2/4 + C^2*F^2) where F is field
// C= Linear coefficient is in cm-1/T
// Energy and Delta are in cm-1

/*** Program to read Line or Energy list from Pgopher ***/
// PGOPHER Gives
// Molecule Manifold   M  Sym  #  g Population Label  Energy Linear   Dipole  Err  Quadratic       Err  Two_Level Energy Delta     C              Dipole2    Err
// WE DO NOT USE IT LIKE THAT BUT USING ORIGIN WE MODIFY AND USE ONLY
// Manifold(1 for upper or 0 for lower typically)  M(in Field it is real M)  Sym(Parity +/-)  #(number to discriminate the levels) J N  Energy0(in 0 field)   Delta   C
// Where for fit in Energy in Fields  Energy +/- Sqrt(Delta^2/4 + C^2*F^2) where F is field
// C= Linear coefficient is in cm-1/T
// Energy and Delta are in cm-1

// Retourne le nb de niveaux (1 si contient seulement Level[0])
int Pgopher_Level_List(const char *nom_file, vector <Internal_state> &Level, FitParams &params)
{
    Level.clear();
    ifstream file(nom_file);
    if ( !file )
    {
        cerr << "Erreur d'ouverture fichier Level"  << nom_file << endl;
        return 0;
    }

    int i=0; // Compteur du nombre de niveaux

    int type_of_field_for_internal_state_shift = params.LocateParam("type_of_field_for_internal_state_shift")->val;

    while (!file.eof())
    {
        Internal_state Current_Level;
        if (type_of_field_for_internal_state_shift == 0) // Zeeman effect
            Current_Level.read_Level_Pgopher_B(file); // read the new level
        if (type_of_field_for_internal_state_shift == 10) // Stark effect
            Current_Level.read_Level_Pgopher_E(file); // read the new level
        Level.push_back(Current_Level);
        i++;
        // cout << i << "  "  << Level[i].Energy0_cm << "  " << Level[i].Delta_FieldB << endl;
    }
    file.close();
    return i;
}

/*** LINE LIST *****/
// Read the file containing the position, assignment and intensity of lines in the simulated spectrum.
// Retourne le nb de transitions
//  C'est une liste qui depuis l'état de départ pointe vers tous les états d'arrivés possible
int Pgopher_Line_List(const char *nom_file, vector <Internal_state> &Level, FitParams &params)
{
    int nb_levels = Level.size();

    ifstream file(nom_file);

    if ( !file )
    {
        cerr << "Erreur d'ouverture fichier Lignes" << nom_file << endl;
        return 0;
    }

    Internal_state Up_state, Low_state;
    double Spol,position,intensity; // Intensité de la raie en Debye^2

    int i=0; // Compteur du nombre de transitions

    while (!file.eof())
    {
        read_Line_Pgopher(file, Up_state, Low_state, Spol,position,intensity); // Transition Up_state --> Low_state avec force Spol
// On ajoute ensuite à la liste des raies le vrai niveau car il contient plus d'informations (en particulier les variations en champ E et B)
// que celui donné dans le fichier de liste des raies
        int n_up = -1;
        int n_low = -1;
        for (int j = 0; j < nb_levels; j++) // Boucle pour voir dans la liste des niveaux où sont Up_state et Low_state
        {
            if (Level[j].is_equal (Up_state))
                n_up=j;
            if (Level[j].is_equal (Low_state))
                n_low=j;

        }

        // cout << " Transition entre niveau num. " << n_up << " et " << n_low << endl;
        if ((n_up != -1) && (n_low != -1)) // if we have find the corresponding levels then add the transition
        {
            Level[n_up].add_transition ( &(Level[n_low]), Spol); // Ajoute la transition  Up_state --> Low_state avec force Spol
            Level[n_low].add_transition ( &(Level[n_up]), Spol);
            i++;
        }
    }
    file.close();
    return i;

}



/************************************************************************/
/************************** INITIALISATION FC, Energie, HL **************/
/************************************************************************/

// Initialisation des durées de vies
void initialisation_Gamma( vector <Internal_state> &Level, FitParams &params)
{
    int nb_Levels = Level.size();
    for (int i = 0; i < nb_Levels; i++) // Boucle pour voir dans la liste des niveaux où sont Up_state et Low_state
        Level[i].one_over_lifetime =  Level[i].Einstein_desex_rate();
}



// Initialisation des probabilités de peuplement de départ
void initialisation_proba(const gsl_rng *r, vector <Molecule> &Mol, const int Nb_molecules, const vector <Internal_state> &Level)
{

    int nb_Levels = Level.size();

    double *proba = new double[nb_Levels]; // We can not use   vector <double> proba because of gsl_ran_discrete_preproc

    for (int i=0; i<nb_Levels; i++)
                proba[i]=Level[i].population;


    gsl_ran_discrete_t * g = gsl_ran_discrete_preproc(nb_Levels, proba); // crée un générateur de nombre aléatoire selon les poids donné par le tableau proba[]

    for (int nb = 0; nb < Nb_molecules; nb ++)
    {
        int numero = gsl_ran_discrete(r, g); // Numéro du niveau tiré au Hazard
        Mol[nb] = Level[numero]; // The vector Mol has been created before so no push_pack() it is a modification of the content !. c'est l'état interne de la molecule

//  DEBUG          Mol[nb] = Level[nb]; // The vector Mol has been created before so no push_pack() it is a modification of the content !. c'est l'état interne de la molecule
//    cerr << "ATTENTION A ENLEVER dans initialisation_proba test pour Ps " << endl;
    }
    gsl_ran_discrete_free (g);

    delete[] proba;
    return;
}

