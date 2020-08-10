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

    if (( (int) params.LocateParam("if_fixed_lifetime")->val) == ((int) true))
    {
        for (int i = 0; i < nb_Levels; i++)
        {
            if ( Level[i].exc > 0) // Not ground state (0) neither continuum ones (<0)
                Level[i].one_over_lifetime = 1./params.LocateParam("duree_vie")->val;
        }
    }
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
            Mol[nb] = Level[nb]; // The vector Mol has been created before so no push_pack() it is a modification of the content !. c'est l'état interne de la molecule
    cerr << "ATTENTION A ENLEVER dans initialisation_proba test pour Ps " << endl;
    }
    gsl_ran_discrete_free (g);

    delete[] proba;
    return;
}




/******************************************************************************************/
/** Modification des fichiers de niveaux et de transitions: FC et des Energies           **/
/******************************************************************************************/


void Modification_FC_trans_mol(const char *nom_file_Levels, const char *nom_file_Lines , vector <Internal_state> &Level, FitParams &params, DataCards data)
{

    /*******************************************************/
    /******************** Molecular data (cm-1) ************/
    /*******************************************************/


    /** Nom files **/

    const char *nom_file_FC_vAvX = data.SParam("nom_file_FC_vAvX").c_str(); // File containing  FC factors between  vA (colonne) and vX (lignes)
    const char *nom_file_E_vA = data.SParam("nom_file_E_vA").c_str(); // File containingles energies of the A state (vA, E_vA, Bv)
    const char *nom_file_E_vX = data.SParam("nom_file_E_vX").c_str(); // File containingles energies of the A state (vX, E_vX, Bv)

    /*********** Declaration of datas

    BE CAREFUL !!!!!!!!!!!!!!!
    ORDER FOR THE NOTATIONS
     [A][X] and [v][J]
    so:

       EvX[vX][JX] // Energy of vX, JX states
       EvA[vA][JA] // Energy of vA, JA states

       *************/


    const int NXmax=params.LocateParam("NXmax")->val;       // Nb of vibrational levels in X state; 0,1,...,NXmax-1
    const int NAmax=params.LocateParam("NAmax")->val;       // Nb of vibrational levels in A state; 0,1,...,NAmax-1

    const int NX_out=params.LocateParam("NX_out")->val;       // Nb of vibrational levels in X state (dans le fichier final)
    const int NA_out=params.LocateParam("NA_out")->val;       // Nb of vibrational levels in A state (dans le fichier final)



    const int N_Two_JX_out_max = params.LocateParam("N_Two_JX_out_max")->val;       // Nb of vibrational levels used to calculate X state; 0,1,...,NX-1
    const int N_Two_JA_out_max = params.LocateParam("N_Two_JA_out_max")->val;       // Nb of vibrational levels used to calculate A; 0,1,...,NA-1


    double **EvX=new double*[NXmax];
    for (int vX=0; vX<NXmax; vX++)
    {
        EvX[vX]=new double[N_Two_JX_out_max];
    }

    double **EvA=new double*[NAmax];
    for (int vA=0; vA<NAmax; vA++)
    {
        EvA[vA]=new double[N_Two_JA_out_max];
    }

    double **FC=new double*[NAmax];
    for (int vA=0; vA<NAmax; vA++)
    {
        FC[vA]=new double[NXmax];
    }



    /*****************************************************************/
    /******* Coeur du programme                              *********/
    /*****************************************************************/


    initialisation_energie(EvX, nom_file_E_vX, NXmax, N_Two_JX_out_max); // Initialisation avec Energie v=0 J=0 à zéro
    initialisation_energie(EvA, nom_file_E_vA, NAmax, N_Two_JA_out_max); // Initialisation avec Energie v=0 J=0 à zéro
    initialisation_FC(FC, nom_file_FC_vAvX, NAmax, NXmax);

    New_Pgopher_Level_Line(nom_file_Levels, nom_file_Lines, Level, FC, EvX, EvA, NX_out, NA_out, N_Two_JA_out_max, N_Two_JX_out_max, params);

    /*****************************************************************/
    /******* Destruction of variables and dynamic arrays *********/
    /*****************************************************************/


    for (int vX=0; vX<NXmax; vX++)
    {
        delete[] EvX[vX];
    }
    delete[] EvX;

    for (int vA=0; vA<NAmax; vA++)
    {
        delete[] EvA[vA];
    }
    delete[] EvA;

    for (int vA=0; vA<NAmax; vA++)
    {
        delete[] FC[vA];
    }
    delete[] FC;

}



// Initialization Of FC factors for absorption and emission
void initialisation_FC(double **FCondon, const char *nom_file_FC, const int NvAmax, const int NvXmax)
{
    ifstream file_FC(nom_file_FC); // The file format is FC[vA][vX], vA lines, vX columns

    for (int vA = 0; vA < NvAmax; vA++)
        for (int vX = 0; vX < NvXmax; vX++)
        {
            file_FC >> FCondon[vA][vX] ; // Read the file FC_vA_vX
            // cout  << "vX = " << vX << " vA = " << vA <<  " FC_abs[vA][vX] = " << FCondon[vA][vX] << endl;
        }

    file_FC.close();

    return;
}


// Energy of vibrational levels
// The file format is v Ev Bv
void initialisation_energie(double **Ev, const char *nom_file_E_v, const int Nvmax, const int N_Two_J_max, const double E_v0_J0)
{

    ifstream file_EV(nom_file_E_v); // The file format is v Ev Bv

    for (int vi = 0; vi < Nvmax; vi++)
    {
        int v;
        double Bv;
        file_EV >> v;
        file_EV >> Ev[v][0];  // We calculate the  hypothetic energy of J=0
        file_EV >> Bv;


        for (int two_J = 0; two_J < N_Two_J_max; two_J++)
        {
            Ev[v][two_J] = Ev[v][0]+Bv*(two_J/2.)*(two_J/2.+1.); // Rotationel levels Bv J (J+1)
        }
    }

    double decalage = E_v0_J0 - Ev[0][0];
    for (int v = 0; v < Nvmax; v++)
        for (int two_J = 0; two_J < N_Two_J_max; two_J++)
            Ev[v][two_J] += decalage; // Dacay of energy in order to get the right origin of energy

    file_EV.close();

    return;
}

#include <sstream>  // pour ostringstream
#include <iomanip>  // pour steprecision
void New_Pgopher_Level_Line(const char *nom_file_initial_Level, const char *nom_file_initial_Line, vector <Internal_state> &Level, double **FC, double **EvX, double **EvA, int NXout, int NAout, int N_Two_JA_out_max, int N_Two_JX_out_max, FitParams &params)
{
    int nb_Levels = Pgopher_Level_List(nom_file_initial_Level, Level, params);

    std::ostringstream oss;
    oss << "_vX_" << NXout << "_vA_" << NAout << "_2JX_" << N_Two_JX_out_max-1 << "_2JA_" << N_Two_JA_out_max-1;
    std::string num = oss.str(); // num = "[i]"

    string nom_file_short = (string(nom_file_initial_Level)).substr(0,string(nom_file_initial_Level).length() -4); // enlève le .dat
    string nom_file_long = nom_file_short + num + ".dat";
    const char *nom_file_new_Level_Out = (nom_file_long).c_str(); // Nom_file[numero_laser].dat
    ofstream file_Level_out(nom_file_new_Level_Out);   // Fichier contenant les données
    file_Level_out << setprecision(8); // POur avoir bien les tous les chiffres significatifs


    nom_file_short = (string(nom_file_initial_Line)).substr(0,string(nom_file_initial_Line).length() -4); // enlève le .dat
    nom_file_long = nom_file_short + num + ".dat";
    const char *nom_file_new_Line_Out = (nom_file_long).c_str(); // Nom_file[numero_laser].dat
    ofstream file_Line_out(nom_file_new_Line_Out);   // Fichier contenant les données
    file_Line_out << setprecision(8); // POur avoir bien les tous les chiffres significatifs

    /*********************  LEVELS ***************************/

    double Bv0X= (EvX[0][2]-EvX[0][0])/2.;   // entre J=1 (two_J=2) et J=0 on a Bv J(J+1) = 2 Bv
    double Bv0A= (EvA[0][2]-EvA[0][0])/2.;   // entre J=1 (two_J=2) et J=0 on a Bv J(J+1) = 2 Bv

    int Max_deg_Number = 10000; // Sécurité pour le max de dégénérescence
    for (int i=0; i<nb_Levels; i++) // ON lit les Level et on ajoute tous les vX si le niveaux est X et les vA si A.
    {
        Internal_state New_Level =  Level[i];

        if (Level[i].exc ==0)  // etat X
            for (int vX = 0; vX < NXout; vX++)
            {
                New_Level.v=vX;
                New_Level.deg_number = Level[i].deg_number + Max_deg_Number*vX; // Nouvelle dégénérescence (on met une sécurité avec LARGE_NUMBER qui est le max possible de la dégénéréscence)
                New_Level.Energy0_cm = Level[i].Energy0_cm - Bv0X* (Level[i].two_J/2.)*( Level[i].two_J/2.+1.) + EvX[vX][Level[i].two_J] ; // Rotationel levels Bv J (J+1). Permet de conserver la vrai progression pour v=0. Ensuite c'est le Bv théorique
                if (vX != 0)  New_Level.population = 0.; //On met les populations des v>0 à 0
                if (Level[i].two_J < N_Two_JX_out_max) New_Level.write_Level_PgopherB(file_Level_out);
            }

        if (Level[i].exc == 1)  // etat A
            for (int vA = 0; vA < NAout; vA++)
            {
                New_Level.v=vA;
                New_Level.deg_number = Level[i].deg_number + Max_deg_Number*vA; // Nouvelle dégénérescence (on met une sécurité avec Nb_niveau qui est le max possible de la dégénéréscence)
                New_Level.Energy0_cm = Level[i].Energy0_cm + EvA[vA][Level[i].two_J] - Bv0A* (Level[i].two_J/2.)*( Level[i].two_J/2.+1.);
                if (vA != 0)  New_Level.population = 0.;
                if (Level[i].two_J < N_Two_JA_out_max) New_Level.write_Level_PgopherB(file_Level_out);
            }
    }

    file_Level_out.close();

    /*********************  LIGNES ***************************/

    Internal_state Up_state, Low_state, New_Up, New_Low;
    double Spol, New_Spol;
    double position, New_position, Intensity, New_Intensity;

    ifstream file(nom_file_initial_Line);

    while (true)
    {
        read_Line_Pgopher(file, Up_state, Low_state, Spol, position, Intensity);

// On cherche le vrai niveau en parcourant la liste des niveaux car il contient plus d'informations (en particulier J)
// que celui donné dans le fichier de liste des raies
        int n_up=0;
        int n_low=0;
        for (int j = 0; j < nb_Levels; j++) // Boucle pour voir dans la liste des niveaux où sont Up_state et Low_state
        {
            if (Level[j].is_equal (Up_state))
                n_up=j;
            if (Level[j].is_equal (Low_state))
                n_low=j;
        }

        Up_state = Level[n_up]; // Copie toutes les infos du niveau
        Low_state = Level[n_low];

        for (int vX = 0; vX < NXout; vX++)
            for (int vA = 0; vA < NAout; vA++)
            {
                New_Low = Low_state;
                New_Low.v=vX;
                New_Low.deg_number = Low_state.deg_number + Max_deg_Number*vX; // Nouvelle dégénérescence (on met une sécurité avec LARGE_NUMBER qui est le max possible de la dégénéréscence)
                New_Low.Energy0_cm = Low_state.Energy0_cm + EvX[vX][Low_state.two_J] - Bv0X* (Low_state.two_J/2.)*( Low_state.two_J/2.+1.); // Rotationel levels Bv J (J+1). Permet de conserver la vrai progression pour v=0. Ensuite c'est le Bv théorique
                if (Low_state.two_J >= N_Two_JX_out_max) continue; // Si pas le bon niveau on sort

                New_Up = Up_state;
                New_Up.v=vA;
                New_Up.deg_number = Up_state.deg_number + Max_deg_Number*vA; // Nouvelle dégénérescence (on met une sécurité avec Nb_niveau qui est le max possible de la dégénéréscence)
                New_Up.Energy0_cm = Up_state.Energy0_cm + EvA[vA][Up_state.two_J] - Bv0A* (Up_state.two_J/2.)*( Up_state.two_J/2.+1.);
                if (Up_state.two_J >= N_Two_JA_out_max) continue;

                New_Spol = Spol*FC[vA][vX]/FC[0][0];
                New_position =  New_Up.Energy0_cm- New_Low.Energy0_cm;
                New_Intensity = Intensity * FC[vA][vX]/FC[0][0];

                if (New_Spol < sqrt(SMALL_NUMBER)) continue; // On ne prend que les transition importantes dipole ~ sqrt(Spol)

                write_Line_Pgopher(file_Line_out, New_Up, New_Low,  New_Spol,  New_position, New_Intensity);
            }
        if ( file.eof() )
        {
            break;      // test de fin de fichier
        }
    }
    file_Line_out.close();
    file.close();
}

