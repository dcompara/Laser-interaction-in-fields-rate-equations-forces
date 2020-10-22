/*
  Name: Laser cooling of translational degree of freedom of a molecule.


  Author: Daniel Comparat
  Date: 17/12/08 (mise � jour 25/4/2012)

  Compilateur : Code::Blocks
  Librairy Math: GSL (GNU)
  Affichage: OPEN GL + GLUT


les mol�cules sont d�crites par une classe Molecule
leur param�tres (nb, vitesse, etc.. sont dans le fichier Liste_Param.h) (ce n'est pas un fichier C++ le .h est l� juste pour avoir un bon �diteur))
Les transtions mol�culaires niveaux, transitions (Franck-COndon, H�nl London,...) sont lus
dans Initialisation_programme qui appelle les fichiers dans Transitions_initialisation

Les taux de transitions proviennent d'une liste qui souvent vient du programme
Ainsi les unit�s sont
DEBYE^2 pour force de raie (dipole^2) et
CM^-1 pour �nergie

Elles sont excit� par laser (classe Laser)
qui excite le nuage gaussien  (taille sigma_SAMPLEx,y,z) via plusieurs
laser eux aussi suppos�s gaussiens (waist_x,y,z). Voir d�tail dans chaque .h des classes)

Il y a un mouvement dans un champ (�lectrique ou magn�tique) d�crit par une classe Field

Le programme Transitions_rate_calcul calcul les taux de transition (pas d'�quations de Bloch seulement de taux).
Il est bas� sur l'interaction avec un laser de largeur fini.
La saturation est (un peu) prise en compte.


l'�volution temporelle (de l'�tat interne) est bas�e sur une �quation de taux bien adapt�
� l'�volution cin�tique Monte Carlo faite dans le programme Kinetic_Monte_Carlo

Plus pr�cisement on calcul un temps dt_KMC pour l'�volution de l'�tat interne
il faut �videmment qu'il soit petit devant le temps dt_dyn caract�ristique d'une petite �volution des taux
(par exemple le temps de transit dans le laser qui modifie les taux d'excitation)
Dans le cas contraire on fait �voluer les positions des particules d'une fraction de dt_dyn
sous l'effets des forces Zeeman, Stark et dipolaires selon un algorithme "V�locity Verlet" dans "one_body"
et on recommence.

La mise � jours des �tats des mol�cules et de leur �nergie ce fait dans "shift_molecule".

La sortie est graphique  dans "Affichage" et dans des fichiers de donn�es dans "sortie_donnees"

Des statistiques peuvent �tre faites sur les vitesses et positions (temp�ratures, �nergie) dans "Stat_Molecule"

Il y a des
TODO:  l� o� il faut faire des choses mieux
DEBUG l� o� je met des fonctions utilis�es pour le d�buggage actuel

*/

#include  <iostream>                       // to include cout, cin
#include  <cmath>                          // to include sqrt(), etc.
#include  <cstdlib>                        // for atoi() and atof()
#include  <stdlib.h>
#include   <stdio.h>
#include  <unistd.h>                       // for getopt()
#include  <fstream>                        // to read the data and put them in files
#include <algorithm>                       // Pour le binary search
#include <numeric>                         // Pour accumulate
#include <iomanip>
#include <list>                            // for list


#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include  <gl/glut.h>                       // GLUT bibliotheque
#include <GL/glu.h>                         // OPENGL (GLU)
#include <GL/gl.h>                          // OPENGL
#include <GL/glext.h>

//#include "GltZpr/zpr.h"               //Pour zoomer, translater et tourner avec la souris
// Ne marche pas actuelement
//#include "GltZpr/zpr.c"               //Pour zoomer, translater et tourner avec la souris


#include <gsl/gsl_rng.h>                   // for random generator
#include <gsl/gsl_randist.h>               // for gaussian random generator


#include "algorithmes.h"                       // Pour le binary search et le wait
#include "Internal_state.h"                       // Classe �tat interne
#include "constantes_SI.h"                  // SI Constantes
#include "Kinetic_Monte_Carlo.h"            // Algorithme MC cin�tique
#include "Laser.h"                       // Classe Laser
#include "Molecule.h"                       // Classe molecule
#include "Field.h"                       // Classe molecule
#include "datacards.h"                      // Pour lire le fichier de param�tres
#include "params.h"                         // Pour faire varier des param�tres

#include "Transitions_initialisation.h"       // Initialisation du programme
#include "Initialisation_programme.h"       // Initialisation du programme
#include  "shift_molecule.h"              // Pour le shift en �nergie delta
#include  "diagonalization.h"              // Pour le shift en �nergie delta
#include "Transition_rate_calcul.h"                     // Sauvegarde des donn�es (+ sortie diverses)
#include "Stat_Molecule.h"                       // Classe molecule
#include "Affichage.h"                     // Affichage �cran
#include "Affichage_Mol.h"                     // Affichage des mol�cules
#include "sortie_donnees.h"                     // Sauvegarde des donn�es (+ sortie diverses)
#include "one_body.h"                   // Algorithme pour calculer les d�placements
#include "scaling_velocities.h"

using namespace std;


// fonction appel�e en cas de d�passement de m�moire
void deborde ()
{
    cerr << "m�moire insuffisante - arr�t de l'ex�cution " << endl;
    exit (1);
}


/************************************************************************/
/************************** Programme Principal *************************/
/********  main appelle RePaint () dans OPENGL + GLUT *******************/
/************************************************************************/


// fonction g�n�rale: initialisation + boucle KMC et N-corps
void RePaint ()
{
    /*****************************************************************/
    /** Cr�ation des variables utiles (param�tres � scanner ou pas) **/
    /*****************************************************************/

    string nomdat = "Data/Liste_Param.h" ;  // nom du fichier de param�tres
    DataCards data(nomdat.c_str()); // Lit le fichier et cr�e les datacards.

    bool Graphics = (bool) data.IParam("Graphics"); // affichage graphique ou non.
    const double SIZE_affichage = data.DParam("SIZE_affichage"); // Taille de la zone d'affichage
    MC_algorithmes Algorithme_MC = (MC_algorithmes) data.IParam("Choix_algorithme_Monte_Carlo"); // Choix de l'algorithme Mont� Carlo
    N_Body_algorithmes Algorithme_N_body = (N_Body_algorithmes) data.IParam("Choix_algorithme_N_corps"); // Choix de l'alogithme N coprs

    /***  NOM DES FICHIERS ***/
    string nom_sortie_temp_string, nom_sortie_scal_string, nom_file_Levels_string, nom_file_Lines_string, nom_sortie_donnees_string, nom_sortie_rate_string, nom_fichier_random_gen_string,
           nom_file_Laser_Spectrum_string, nom_file_Laser_Intensity_string, nom_file_Magn_Field_3D_string, nom_file_Elec_Field_3D_string;

    const char *nom_sortie_temp, *nom_sortie_scal, *nom_file_Levels, *nom_file_Lines, *nom_sortie_donnees, *nom_sortie_rate, *nom_fichier_random_gen,
          *nom_file_Laser_Spectrum, *nom_file_Laser_Intensity, *nom_file_Magn_Field_3D, *nom_file_Elec_Field_3D;

    nom_file_Levels_string =  data.SParam("nom_file_Levels");      // Fichier contenant les Levels (etat, energie, ...)
    nom_file_Lines_string = data.SParam("nom_file_Lines");         // Fichier contenant les transitions
    nom_sortie_donnees_string = data.SParam("nom_sortie_donnees");
    nom_sortie_rate_string = data.SParam("nom_sortie_rate");
    nom_fichier_random_gen_string = data.SParam("nom_fichier_random_gen");
    nom_file_Laser_Spectrum_string = data.SParam("nom_file_Laser_Spectrum");         // Fichier contenant les transitions
    nom_file_Laser_Intensity_string = data.SParam("nom_file_Laser_Intensity");         // Fichier contenant les transitions
    nom_file_Magn_Field_3D_string =  data.SParam("nom_file_Magn_Field_3D");
    nom_file_Elec_Field_3D_string =  data.SParam("nom_file_Elec_Field_3D");

    nom_file_Levels =  nom_file_Levels_string.c_str();      // Fichier contenant les Levels (etat, energie, ...)
    nom_file_Lines =  nom_file_Lines_string.c_str();         // Fichier contenant les transitions
    nom_sortie_donnees = nom_sortie_donnees_string.c_str();
    nom_sortie_rate = nom_sortie_rate_string.c_str();
    nom_fichier_random_gen = nom_fichier_random_gen_string.c_str();
    nom_file_Laser_Spectrum = nom_file_Laser_Spectrum_string.c_str();         // Fichier contenant les transitions
    nom_file_Laser_Intensity = nom_file_Laser_Intensity_string.c_str();         // Fichier contenant les transitions
    nom_file_Magn_Field_3D =  nom_file_Magn_Field_3D_string.c_str();
    nom_file_Elec_Field_3D =  nom_file_Elec_Field_3D_string.c_str();

    FitParams params; // this is a vector<Param *>
    params.read(nomdat); // Lit le ficher des param�tres entre "BEGIN_OF_FITPARAMS" et "END_OF_FITPARAMS".
    params.init_value(nomdat); // Initialise les param�tres selon la valeur non scann�e de la dataCard (cf datacards.h) ou � la valeur min si scann�e


    /******************************************************************/
    /** cr�ation abstraite des objets (champs, mol�cules, laser ...) **/
    /******************************************************************/

    vector <Internal_state> Level; // Level is a list of Internal_state
    vector <Laser> lasers; //Lasers utilis�s pour le refroidissement

    Field champB, champE; //champs magn�tique et �lectrique

    const int Nb_type_of_Mol =  data.IParam("Nb_type_of_Mol");
    int N_Mol[1];
    N_Mol[0] = (int)  params.LocateParam("N_Mol[0]")->val; // The number of laser cooled molecules
    vector <Molecule> Mol;

    vector <double> rate; // Taux (int�gr� en temps) des transitions.
    vector <type_codage_react> reaction_list; // liste des r�actions selon le format donn� par type_codage_react souvent (n� du syst�me, �tat final du syst�me)

    ofstream file_out(nom_sortie_donnees); // ajouter ", ios::app" si on veux �crire � la fin du fichier. Sinon efface le fichier � la r�ecriture
    file_out<< setprecision(8);             // 8 d�cimales dans les fichiers de sortie
    ofstream file_rate(nom_sortie_rate);

    clock_t t_start,t_end; // Pour mesurer le temps de d�roulement du programme
    t_start = clock();


    /*****************************************************************/
    /********************** Initialisation du GSL ********************/
    /*****************************************************************/

    gsl_rng * r; // G�n�rateur de nb al�atoire
    // set_random_generator(r, (int) params.LocateParam("Seed_Init_Random_Number_Generator")->val, nom_fichier_random_gen);  // Initialize the random_generator (Mersenne Twister here)
// TODO (Daniel#8#): DOes not work (file fo rrandom generator). But averything is almost done, so should be easy to fix

    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;  // "Mersenne Twister" generator by default. Plus rapide que RANLUX. cf. chapitre 17.12 de gsl_ref.
    r = gsl_rng_alloc (T);

    initialisation_trans_mol(nom_file_Levels, nom_file_Lines, Level, params); // Lecture des fichiers de niveaux et de transitions

    bool  test_fin_boucle_param; // Param�tre de fin de boucle sur les param�tres � scanner

    if (data.IParam("is_DataCard_out") == 1)
    {
        file_out << " Nb card " << data.NbCards() << endl << "DATA " << endl << data << endl;
    }
    // AU DEBUT SORT TOUTES LES DONNEES du fichier de donn�e initiale
    if (data.IParam("is_DataCard_out") == 2)
    {
        ofstream file_out_datacard(data.SParam("nom_sortie_donnees_Data").c_str()); // ajouter ", ios::app" si on veux �crire � la fin du fichier. Sinon efface le fichier � la r�ecriture
        file_out_datacard << " Nb card " << data.NbCards() << endl << "DATA " << endl << data << endl; // AU DEBUT SORT TOUTES LES DONNEES du fichier de donn�e initiale
        file_out_datacard.close();
    }

    do // boucle sur les param�tres � scann�s (voir fonction Modif_Param)
    {
        double t=0.; //temps de d�but de la simulation en �s
        double t_mise_a_jour=0.; //temps de faire la mise a jour (potentiels, acc�l�ration, N corps ..) de toutes les mol�cules
        /** On  ne met � jour les potentiels que de la mol�cule dont l'�tat interne � �t� modifi� ou lorsque le pas temporel du N-corps est atteind on les modifie toutes **/
        double t_dia=0.; //temps de sortie fichier
        double t_out=0.; //temps d'affichage
        double t_fin =  params.LocateParam("t_fin")->val; // temps final
        double dt_KMC; // Temps pour qu'une r�action apparaisse
        double dt_dyn = 0; // temps dynamique
        double dt_dia = params.LocateParam("dt_dia")->val; // time interval between diagnostics (in cout) output
        double dt_out = params.LocateParam("dt_out")->val; // time interval between output of snapshots (draw particles)

        /*****************************************************************/
        /********************** Initialisation du GSL ********************/
        /*****************************************************************/

        gsl_rng_set(r,(int) params.LocateParam("Seed_Init_Random_Number_Generator")->val);
        save_random_generator(r, nom_fichier_random_gen);

        /*****************************************************************/
        /********************** Initialisation des mol�cules *************/
        /*****************************************************************/


        Init_Field(champB, champE, params, nom_file_Magn_Field_3D, nom_file_Elec_Field_3D);  // Initialise les champs.
        // Create the list of molecules (create the vector Mol here) its size will not be modified then
        Init_Molecule(r, Mol, champB, champE, Nb_type_of_Mol, params, data); // Position, vitesse
        initialisation_proba(r, Mol, N_Mol[0], Level); // Etat des populations de mol�cules au d�part en fonction de la population voulue

        const int Nb_laser = data.IParam("Nb_laser");  // number of used laser (we could have more in the Liste_Param file, but this will be the number used in this run)

        Init_Laser(lasers,Nb_laser, params, nom_file_Laser_Spectrum, nom_file_Laser_Intensity); // Initialise les lasers
        // Sortie_laser_spectrum(file_out, lasers, params,0); // Debug
        // Sortie_laser_intensity(file_out, lasers, params,0);

        /**
        on calcul un temps dt_KMC pour l'�volution de l'�tat interne.
        Si dt_KMC < temps dt_evol_ext_typ caract�ristique de l'�volution des taux (par exemple le temps de transit dans le laser qui modifie les taux d'excitation) on effectue la transition et on bouge les particules
        Dans le cas contraire on fait �voluer les particules durant dt_dyn qui est d'une fraction de dt_evol_ext_typ et on recommence.
        **/

        int number_mol = aucune;  // num�ro de la mol�cule affect�e par une modification (aucune au d�but d'o� la valeur -1)
        int number_photons = 0;  // nombre de photons en jeu (absorption, emission spontan�e ou stimul� ...

        //Create_dipole_Lines_from_Matrices("matrice_dipole.dat");

        while(true) // Infinite loop untill t_end is reached
        {
            Init_Field(champB, champE, params);  // Re0-Initialise les champs. Important si scan des param�tres car ils ont chang�s. We do not add the files because we do not modify them
            Init_Laser(lasers, Nb_laser, params, nom_file_Laser_Spectrum, nom_file_Laser_Intensity); // Initialise les lasers. // On pourait ne pas remetre � jour le fichier des niveaux
            dt_dyn = (params.LocateParam("dt_dyn_epsilon_param")->val);

            calcul_rates_molecules(Level, Algorithme_MC, reaction_list, rate, Mol, champB, champE, lasers, t, number_mol, N_Mol[0], params); // Calcul les taux de transition de toutes les mol�cules si numero_mol = aucune. Sinon on ne recalcule que celui de la mol�cule numero_mol
            if (t >= t_dia)
            {
                /*** Here, or just before if you want to have output all the time, put you output files ***/
                // I suggest that you look at the  Sortie_rate_example and Sortie_donnee_example in the sortie_donnees.cpp to inspire you for the Sortie_rate and Sortie_donnee or Sortie_donnee_pop_v file or whatever you want to have such as Sortie_laser_spectrum ***/

                // Sortie_rate(file_rate, rate, Level, reaction_list, Mol, champB, champE, lasers, N_Mol[0], t, params);
                // Sortie_donnee(file_out, Mol, Level, champB, champE, lasers, t, (int) Mol.size(),params,  data, number_photons);  // sortie de toutes les donn�es mol�culaires
                t_dia += dt_dia;
            }

            int n_reac = find_reaction(Algorithme_MC, r, rate,  dt_KMC); // Trouve la r�action KMC

            if (dt_KMC < dt_dyn || Algorithme_MC == Fast_Rough_Method) // L'�volution se fera durant un temps dt_KMC or if Fast_Rough_Method in order to be sure to evolve the internal states (because it is not random)
            {
                evolve_step(Algorithme_N_body, Mol, champB, champE, lasers, Mol.size(), t, dt_KMC, params); // Evolution N corps (avant le changement de vitesse). Cela change aussi le temps;
                number_mol = do_reaction(Algorithme_MC, r, rate, reaction_list, Mol, n_reac, lasers, dt_KMC, file_rate, true, number_photons, params); // On fait l'�volution de l'�tat interne
            }
            else // L'�volution se fera sans r�action durant un temps dt_dyn
            {
                evolve_step(Algorithme_N_body, Mol, champB, champE, lasers, Mol.size(), t, dt_dyn, params); // Evolution N corps (avant le changement de vitesse)
            }

// On peut peut �tre am�liorer. Ne faut'il pas mettre � jour toujours la mol�cule modif�e surtout dans le cas dipolaire ?
            if (t >= t_mise_a_jour || Algorithme_MC == Fast_Rough_Method)
            {
                number_mol = aucune; // Il faudra faire une mise a jour totale (cf calcul rate) donc que toutes les mol�cules ont �t� affect�es
                // set_pot_all_mol(Mol, champB, champE, lasers, t, N_Mol[0], params); //Met � jour (pour recalculer les bonnes transitions) de tous les potentiels (gravit�, dipolaire, magn�tique, electrique, ...) avec la nouvelle position
                t_mise_a_jour += dt_dyn;
            }
// On pourrait penser acc�lerer, en faisant �voluer le N_corps jusqu'a dt_KMC sans recalculer dt_KMC � chaque fois. Mais cela est risqu� comme on le voit avec un mol�cule loin d'un waist --> tKMC immense mais qui va diminuer vite et si on ne recalcule pas on va faire une erreur

            if (Graphics && t >= t_out)
            {
                Draw(Mol, Level,  champB, champE, lasers, SIZE_affichage, t, Mol.size(), Nb_type_of_Mol, params);       // Affichage des points
                t_out += dt_out;
                wait(params.LocateParam("t_wait_affichage")->val);   // Permet de ne pas avoir un affichage trop rapide
            }

            if (t > t_fin) // FIN DE LA BOUCLE
            {
                for (ParamIterator i=params.begin(); i != params.end(); ++i)
                {
                    Param &p = **i;
                    if( p.is_scanned == true )
                        cout << " " << p.val_t0  ;
                }
                cout << endl;
                break;
            }

            Modif_Param(t, params); // On modifie en dynamique la longueur d'onde le waist etc... des lasers en fonction du temps
        }


        if (data.SParam("is_Scan_Random") == "true")
            test_fin_boucle_param = params.Scan_Param_aleatoire(r); // Modification des param�tres
        else
            test_fin_boucle_param = params.Scan_Param();
    }

    while (!test_fin_boucle_param); // test de fin des boucles sur les param�tres

    file_out.close();
    file_rate.close();
    t_end = clock();

    cout << "dur�e du programme (s) = " << (t_end - t_start)/double(CLOCKS_PER_SEC) << endl;


// FIN du programme


    exit(1);
    system("PAUSE");
    return;

}




// fonction main appellant RePaint ().
// C'est ainsi que fonctionne OPENGL + GLUT
int main(int argc, char** argv)
{


    string nomdat = "Data/Liste_Param.h" ;  // nom du fichier de param�tres
    DataCards data(nomdat.c_str()); // Lit le fichier et cr�e les datacards.
    bool Graphics = (bool) data.IParam("Graphics"); // affichage graphique ou non.

    if (Graphics)
    {
        int size_screen =600;
        glutInit(&argc, argv);                            //Initialisation de la GLUT
        glutInitDisplayMode(GLUT_DEPTH|                   //On active la profondeur
                            GLUT_DOUBLE|                  //Un seul tampon d'affichage
                            GLUT_RGBA);                   //Couleurs au format RGBA
        glutInitWindowPosition(700, 100);                 // coordonn�es de la fen�tre
        glutInitWindowSize(size_screen, size_screen);                     //taille Fen�tre
        glutCreateWindow("Cooling");   //Cr�ation de la fen�tre
        // zprInit();                          //Pour zoomer, translater et tourner avec la souris. NE marche pas

        glutDisplayFunc(RePaint);                      //On lui dit d'appeler la fonction renderFunc pour afficher


        glutReshapeFunc(reshape);                      // N�cessaire pour recadrer la fonction
        glutMainLoop();                                   //Boucle infinie du programme
        exit(1);
        return 0;
    }
    else
        RePaint ();
    return 0;
    exit(1);
}
