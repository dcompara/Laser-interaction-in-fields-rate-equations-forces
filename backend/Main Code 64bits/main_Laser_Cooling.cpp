<<<<<<< HEAD
/*
  Name: Laser cooling of the translational degree of freedom of a molecule.

  Author: Daniel Comparat
  Date: 17/12/08 (updated on 25/4/2012)

  Compiler: Code::Blocks
  Math Library: GSL (GNU)
  Visualization: OPEN GL + GLUT

  Molecules are described using a class called Molecule.
  Their parameters (number, velocity, etc.) are defined in the file Liste_Param.h
  (note: this is not a C++ file; the .h extension is used for compatibility with text editors).

  Molecular transitions (levels, Franck-Condon, Hönl-London, etc.) are initialized
  in Initialisation_programme, which reads from files defined in Transitions_initialisation.

  Transition rates come from a list, often derived from a specific program.
  Units used are:
  - DEBYE^2 for line strength (dipole^2)
  - CM^-1 for energy.

  Molecules are excited by lasers (Laser class), which interact with the Gaussian-shaped
  particle cloud (defined by sigma_SAMPLEx, y, z) via several lasers, also assumed Gaussian
  (waist_x, y, z). For more details, see the .h files of each class.

  Molecular motion and interactions are calculated...


Molecules move in a field (electric or magnetic) described by the Field class.

The Transitions_rate_calcul program calculates transition rates (not Bloch equations, only rates).
It is based on interactions with a laser of finite linewidth.
Saturation is (partially) taken into account.

The temporal evolution of the internal state is based on a rate equation, well-suited
for the Kinetic Monte Carlo algorithm implemented in the Kinetic_Monte_Carlo program.

More specifically:
- A time step `dt_KMC` is calculated for the internal state evolution.
  This step must be much smaller than `dt_dyn`, the characteristic time for a small change in rates
  (e.g., transit time in the laser, which modifies excitation rates).
- If this condition is not met, the particle positions are evolved by a fraction of `dt_dyn`
  under the effects of Zeeman, Stark, and dipolar forces using a "Velocity Verlet" algorithm
  implemented in the "one_body" module, and the process repeats.

The states and energies of the molecules are updated ...

/*
Updates to molecular states and energies are performed in "shift_molecule."

Output:
- Graphical output is handled in "Affichage."
- Data files are generated in "sortie_donnees."

Statistics on velocities and positions (e.g., temperatures, energy) are calculated in "Stat_Molecule."

Tags in the code:
- TODO: Indicates areas requiring improvement.
- DEBUG: Marks functions used for current debugging.
*/



#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <GL/glut.h>   // GLUT library
#include <GL/glu.h>    // OpenGL (GLU)
#include <GL/gl.h>     // OpenGL
#include <GL/glext.h>

// #include "GltZpr/zpr.h"   // For zooming, translating, and rotating with the mouse
// Currently not working
// #include "GltZpr/zpr.c"   // For zooming, translating, and rotating with the mouse

#include <gsl/gsl_rng.h>       // For random number generator
#include <gsl/gsl_randist.h>   // For Gaussian random generator

#include "algorithmes.h"       // For binary search and waiting algorithms
#include "Internal_state.h"    // Internal state class
#include "constantes_SI.h"     // SI constants
#include "Kinetic_Monte_Carlo.h" // Kinetic Monte Carlo algorithm
#include "Laser.h"                   // Laser class
#include "Molecule.h"                // Molecule class
#include "Field.h"                   // Field class
#include "datacards.h"               // For reading the parameter file
#include "params.h"                  // For varying parameters

#include "Transitions_initialisation.h" // Initialization of transitions
#include "Initialisation_programme.h"   // General program initialization
#include "shift_molecule.h"          // For energy shifts (delta)
#include "diagonalization.h"         // For energy shifts (delta)
#include "Transition_rate_calcul.h"  // Calculates transition rates
#include "Stat_Molecule.h"           // For molecule statistics
#include "Affichage.h"               // Screen display
#include "Affichage_Mol.h"           // Molecule display
#include "sortie_donnees.h"          // Data output


using namespace std;


// Function called in case of memory overflow
void deborde()
{
    std::cerr << "Insufficient memory - execution stopped" << std::endl;
    exit(1);
}


/************************************************************************/
/************************** Main Program *******************************/
/******** main calls RePaint() in OPENGL + GLUT ***********************/
/************************************************************************/

// General function: initialization + KMC and N-body loop
void RePaint()
{
    /*****************************************************************/
    /** Creation of useful variables (parameters to be scanned or not) **/
    /*****************************************************************/

    std::string nomdat = "Data/Liste_Param.h";  // Name of the parameter file
    DataCards data(nomdat.c_str()); // Reads the file and creates the datacards.

    bool Graphics = (bool)data.IParam("Graphics"); // Enable or disable graphical output.
    const double SIZE_affichage = data.DParam("SIZE_affichage"); // Size of the display area
    MC_algorithmes Algorithme_MC = (MC_algorithmes)data.IParam("Choix_algorithme_Monte_Carlo"); // Choice of Monte Carlo algorithm
    N_Body_algorithmes Algorithme_N_body = (N_Body_algorithmes)data.IParam("Choix_algorithme_N_corps"); // Choice of N-body algorithm

    /*** FILE NAMES ***/
    std::string nom_sortie_temp_string, nom_sortie_scal_string, nom_file_Levels_string, nom_file_Lines_string, nom_sortie_donnees_string, nom_sortie_rate_string, nom_fichier_random_gen_string,
        nom_file_Laser_Spectrum_string, nom_file_Laser_Intensity_string, nom_file_Magn_Field_3D_string, nom_file_Elec_Field_3D_string;

    const char *nom_sortie_temp, *nom_sortie_scal, *nom_file_Levels, *nom_file_Lines, *nom_sortie_donnees, *nom_sortie_rate, *nom_fichier_random_gen,
          *nom_file_Laser_Spectrum, *nom_file_Laser_Intensity, *nom_file_Magn_Field_3D, *nom_file_Elec_Field_3D;

    nom_file_Levels_string = data.SParam("nom_file_Levels");      // File containing levels (state, energy, ...)
    nom_file_Lines_string = data.SParam("nom_file_Lines");        // File containing transitions
    nom_sortie_donnees_string = data.SParam("nom_sortie_donnees");
    nom_sortie_rate_string = data.SParam("nom_sortie_rate");
    nom_fichier_random_gen_string = data.SParam("nom_fichier_random_gen");
    nom_file_Laser_Spectrum_string = data.SParam("nom_file_Laser_Spectrum");   // File containing transitions
    nom_file_Laser_Intensity_string = data.SParam("nom_file_Laser_Intensity"); // File containing transitions
    nom_file_Magn_Field_3D_string = data.SParam("nom_file_Magn_Field_3D");
    nom_file_Elec_Field_3D_string = data.SParam("nom_file_Elec_Field_3D");

    nom_file_Levels = nom_file_Levels_string.c_str();      // File containing levels (state, energy, ...)
    nom_file_Lines = nom_file_Lines_string.c_str();        // File containing transitions
    nom_sortie_donnees = nom_sortie_donnees_string.c_str();
    nom_sortie_rate = nom_sortie_rate_string.c_str();
    nom_fichier_random_gen = nom_fichier_random_gen_string.c_str();
    nom_file_Laser_Spectrum = nom_file_Laser_Spectrum_string.c_str();   // File containing transitions
    nom_file_Laser_Intensity = nom_file_Laser_Intensity_string.c_str(); // File containing transitions
    nom_file_Magn_Field_3D = nom_file_Magn_Field_3D_string.c_str();
    nom_file_Elec_Field_3D = nom_file_Elec_Field_3D_string.c_str();

    FitParams params; // This is a vector<Param *>
    params.read(nomdat); // Reads the parameter file between "BEGIN_OF_FITPARAMS" and "END_OF_FITPARAMS".
    params.init_value(nomdat); // Initializes parameters based on the non-scanned value of the dataCard (cf datacards.h) or to the minimum value if scanned.

    /******************************************************************/
    /** Abstract creation of objects (fields, molecules, lasers, ...) **/
    /******************************************************************/

    std::vector<Internal_state> Level; // Level is a list of Internal_state
    std::vector<Laser> lasers;         // Lasers used for cooling

    Field champB, champE;              // Magnetic and electric fields

    const int Nb_type_of_Mol = data.IParam("Nb_type_of_Mol");
    int N_Mol[1];
    N_Mol[0] = (int)params.LocateParam("N_Mol[0]")->val; // The number of laser-cooled molecules
    std::vector<Molecule> Mol;

    std::vector<double> rate;                   // Rates (time-integrated) of transitions.
    std::vector<type_codage_react> reaction_list; // List of reactions according to the format given by type_codage_react, often (system number, final system state).

    std::ofstream file_out(nom_sortie_donnees); // Add ", ios::app" to write to the end of the file; otherwise, it erases the file on rewrite.
    file_out << std::setprecision(8);          // 8 decimals in output files
    std::ofstream file_rate(nom_sortie_rate);

    clock_t t_start, t_end; // For measuring program runtime
    t_start = clock();

    /******************************************************************/
    /********************** GSL Initialization ************************/
    /******************************************************************/

    gsl_rng *r; // Random number generator
    // set_random_generator(r, (int)params.LocateParam("Seed_Init_Random_Number_Generator")->val, nom_fichier_random_gen);  // Initialize the random generator (Mersenne Twister here)
    // TODO (Daniel#8#): Does not work (file for random generator). But everything is almost done, so it should be easy to fix.

    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;  // "Mersenne Twister" generator by default. Faster than RANLUX. See chapter 17.12 of gsl_ref.
    r = gsl_rng_alloc(T);

    initialisation_trans_mol(nom_file_Levels, nom_file_Lines, Level, params); // Reads level and transition files

    bool test_fin_boucle_param; // End-of-loop parameter for scanned parameters

    if (data.IParam("is_DataCard_out") == 1)
    {
        file_out << " Nb card " << data.NbCards() << std::endl << "DATA " << std::endl << data << std::endl;
    }
    // At the beginning, output all data from the initial data file
    if (data.IParam("is_DataCard_out") == 2)
    {
        std::ofstream file_out_datacard(data.SParam("nom_sortie_donnees_Data").c_str()); // Add ", ios::app" to write to the end of the file; otherwise, it erases the file on rewrite.
        file_out_datacard << " Nb card " << data.NbCards() << std::endl << "DATA " << std::endl << data << std::endl; // At the beginning, output all data from the initial data file
        file_out_datacard.close();
    }

    /*** Hamiltonian + dipole matrix element ***/
    MatrixXcd E0_cm;
    MatrixXcd Zeeman_cm_B;
    MatrixXcd H;        // Hamiltonian Matrix. It is a Hermitian matrix, so I use complex, not MatrixXd
    MatrixXd d0[3];
    MatrixXcd d[3];

    Read_Energy_Zeeman_dipole_for_Diagonalization(E0_cm, Zeeman_cm_B, d0); // Initialize energy, Zeeman, and dipole matrices

    /*****************************************************************/
    /** Global loop (over scanned parameters) ***********************/
    /*****************************************************************/

    do
    {
        double t = 0.0;                 // Start time of the simulation in µs
        double t_mise_a_jour = 0.0;     // Time to update potentials, acceleration, N-body, etc., for all molecules
        /** Potentials are only updated for the molecule whose internal state was modified, or when the N-body time step is reached, all are updated. **/
        double t_dia = 0.0;             // Time for file output
        double t_out = 0.0;             // Time for graphical output
        double t_fin = params.LocateParam("t_fin")->val; // Final time
        double dt_KMC;                  // Time for a reaction to occur
        double dt_dyn = 0.0;            // Dynamic time
        double dt_dia = params.LocateParam("dt_dia")->val; // Interval between diagnostics output (in `cout`)
        double dt_out = params.LocateParam("dt_out")->val; // Interval between graphical output snapshots (particle drawing)

        /*****************************************************************/
        /********************** Initialization of GSL ********************/
        /*****************************************************************/
        gsl_rng_set(r, (int)params.LocateParam("Seed_Init_Random_Number_Generator")->val);
        save_random_generator(r, nom_fichier_random_gen);

        /*****************************************************************/
        /******************** Initialization of molecules ***************/
        /*****************************************************************/

        Init_Field(champB, champE, params, nom_file_Magn_Field_3D, nom_file_Elec_Field_3D);  // Initialize fields
        // Create the list of molecules (vector `Mol`), its size will not be modified later
        Init_Molecule(r, Mol, champB, champE, Nb_type_of_Mol, params, data); // Position, velocity
        initialisation_proba(r, Mol, N_Mol[0], Level); // Initialize the populations of molecules at the start based on the desired population

        const int Nb_laser = data.IParam("Nb_laser");  // Number of lasers used (could have more in `Liste_Param`, but this is the number for this run)

        Init_Laser(lasers, Nb_laser, params, nom_file_Laser_Spectrum, nom_file_Laser_Intensity); // Initialize lasers
        // Sortie_laser_spectrum(file_out, lasers, params, 0); // Debug
        // Sortie_laser_intensity(file_out, lasers, params, 0);

        int number_mol = aucune; // Index of the molecule affected by a modification (none at the start, hence value -1)
        int number_photons = 0;  // Number of photons involved (absorption, spontaneous, or stimulated emission)

        /*****************************************************************/
        /** Global time loop *********************************************/
        /*****************************************************************/

        while (true) // Infinite loop until `t_end` is reached
        {
            Init_Field(champB, champE, params); // Re-initialize fields; important if scanning parameters change them
            Init_Laser(lasers, Nb_laser, params, nom_file_Laser_Spectrum, nom_file_Laser_Intensity); // Initialize lasers; could skip re-reading level files

            dt_dyn = params.LocateParam("dt_dyn_epsilon_param")->val;

            calcul_rates_molecules(Level, Algorithme_MC, reaction_list, rate, Mol, champB, champE, lasers, t, number_mol, N_Mol[0], params, H, E0_cm, Zeeman_cm_B, d0, d); // Compute transition rates for all molecules if `numero_mol == aucune`. Otherwise, recompute only for the molecule `numero_mol`

            if (t >= t_dia)
            {
                /*** Insert output logic here or just before for continuous file output ***/
                // Suggested examples: `Sortie_rate_example` and `Sortie_donnee_example` in `sortie_donnees.cpp`
                Sortie_donnee(file_out, Mol, Level, champB, champE, lasers, t, (int)Mol.size(), params, data, number_photons); // Output all molecular data
                t_dia += dt_dia;
            }

            int n_reac = find_reaction(Algorithme_MC, r, rate, dt_KMC); // Find the KMC reaction

            if (dt_KMC < dt_dyn || Algorithme_MC == Fast_Rough_Method)
            {
                evolve_step(Algorithme_N_body, Mol, champB, champE, lasers, Mol.size(), t, dt_KMC, params); // N-body evolution (before velocity update). Also modifies time.
                number_mol = do_reaction(Algorithme_MC, r, rate, reaction_list, Mol, n_reac, lasers, dt_KMC, file_rate, true, number_photons, params); // Perform the internal state evolution
            }
            else // Evolution will proceed without reaction for `dt_dyn`
            {
                evolve_step(Algorithme_N_body, Mol, champB, champE, lasers, Mol.size(), t, dt_dyn, params); // N-body evolution (before velocity update)
            }

// TODO (Daniel#6#): Maybe separate dt_dyn and t_mise_a_jour

            if (t >= t_mise_a_jour || Algorithme_MC == Fast_Rough_Method)
            {
                number_mol = aucune; // A complete update is needed (see rate calculation), meaning all molecules will be affected.
                // set_pot_all_mol(Mol, champB, champE, lasers, t, N_Mol[0], params);
                // Update all potentials (gravity, dipolar, magnetic, electric, etc.) with the new position to recalculate the correct transitions.
                t_mise_a_jour += dt_dyn;
            }

// One could consider accelerating the process by evolving the N-body system up to dt_KMC without recalculating dt_KMC each time.
// However, this is risky, as seen with a molecule far from the waist -> tKMC becomes immense but decreases quickly, and if we don't recalculate, we might make an error.


            if (Graphics && t >= t_out)
            {
                Draw(Mol, Level, champB, champE, lasers, SIZE_affichage, t, Mol.size(), Nb_type_of_Mol, params); // Display points
                t_out += dt_out;
                wait(params.LocateParam("t_wait_affichage")->val); // Prevent excessively fast rendering
            }

            if (t > t_fin) // End of loop
            {
                for (ParamIterator i = params.begin(); i != params.end(); ++i)
                {
                    Param &p = **i;
                    if (p.is_scanned == true)
                        std::cout << " " << p.val_t0;
                }
                std::cout << std::endl;
                break;
            }

            Modif_Param(t, params); // Dynamically modify laser wavelength, waist, etc., based on time
        }

        if (data.SParam("is_Scan_Random") == "true")
            test_fin_boucle_param = params.Scan_Param_aleatoire(r); // Modify parameters randomly
        else
            test_fin_boucle_param = params.Scan_Param();
    }
    while (!test_fin_boucle_param);   // End-of-loop check for parameters

    file_out.close();
    file_rate.close();
    t_end = clock();

    std::cout << "Execution time (s): " << (t_end - t_start) / double(CLOCKS_PER_SEC) << std::endl;
    std::cout << "Have you checked the to-do list, SMALL_NUMBER_RATE, and other parameters?" << std::endl;

    exit(1);
    system("PAUSE");
    return;
}



// fonction main call RePaint ().
// This is how it works in OPENGL + GLUT
int main(int argc, char** argv)
{


    string nomdat = "Data/Liste_Param.h" ;  // nom du fichier de paramètres
    DataCards data(nomdat.c_str()); // Lit le fichier et crée les datacards.
    bool Graphics = (bool) data.IParam("Graphics"); // affichage graphique ou non.

    if (Graphics)
    {
        const int size_screen = 600;
        glutInit(&argc, argv);                            // Initialisation de la GLUT
        glutInitDisplayMode(GLUT_DEPTH |                 // Active la profondeur
                            GLUT_DOUBLE |                // Double buffering
                            GLUT_RGBA);                  // Couleurs au format RGBA
        glutInitWindowPosition(0, 0);                    // Coordonnées de la fenêtre
        glutInitWindowSize(size_screen, size_screen);    // Taille de la fenêtre

        int windowID = glutCreateWindow("Cooling");      // Création de la fenêtre
        if (windowID == 0)
        {
            std::cerr << "Error: Unable to create window!" << std::endl;
            return -1;
        }

        // zprInit();                          //Pour zoomer, translater et tourner avec la souris. NE marche pas

        glutDisplayFunc(RePaint);                      //On lui dit d'appeler la fonction renderFunc pour afficher


        glutReshapeFunc(reshape);                      // Nécessaire pour recadrer la fonction
        glutMainLoop();                                   //Boucle infinie du programme
        exit(1);
        return 0;
    }
    else
        RePaint ();
    return 0;
    exit(1);
}
=======
/*
  Name: Laser cooling of the translational degree of freedom of a molecule.

  Author: Daniel Comparat
  Date: 17/12/08 (updated on 25/4/2012)

  Compiler: Code::Blocks
  Math Library: GSL (GNU)
  Visualization: OPEN GL + GLUT

  Molecules are described using a class called Molecule.
  Their parameters (number, velocity, etc.) are defined in the file Liste_Param.h
  (note: this is not a C++ file; the .h extension is used for compatibility with text editors).

  Molecular transitions (levels, Franck-Condon, Hönl-London, etc.) are initialized
  in Initialisation_programme, which reads from files defined in Transitions_initialisation.

  Transition rates come from a list, often derived from a specific program.
  Units used are:
  - DEBYE^2 for line strength (dipole^2)
  - CM^-1 for energy.

  Molecules are excited by lasers (Laser class), which interact with the Gaussian-shaped
  particle cloud (defined by sigma_SAMPLEx, y, z) via several lasers, also assumed Gaussian
  (waist_x, y, z). For more details, see the .h files of each class.

  Molecular motion and interactions are calculated...


Molecules move in a field (electric or magnetic) described by the Field class.

The Transitions_rate_calcul program calculates transition rates (not Bloch equations, only rates).
It is based on interactions with a laser of finite linewidth.
Saturation is (partially) taken into account.

The temporal evolution of the internal state is based on a rate equation, well-suited
for the Kinetic Monte Carlo algorithm implemented in the Kinetic_Monte_Carlo program.

More specifically:
- A time step `dt_KMC` is calculated for the internal state evolution.
  This step must be much smaller than `dt_dyn`, the characteristic time for a small change in rates
  (e.g., transit time in the laser, which modifies excitation rates).
- If this condition is not met, the particle positions are evolved by a fraction of `dt_dyn`
  under the effects of Zeeman, Stark, and dipolar forces using a "Velocity Verlet" algorithm
  implemented in the "one_body" module, and the process repeats.

The states and energies of the molecules are updated ...

/*
Updates to molecular states and energies are performed in "shift_molecule."

Output:
- Graphical output is handled in "Affichage."
- Data files are generated in "sortie_donnees."

Statistics on velocities and positions (e.g., temperatures, energy) are calculated in "Stat_Molecule."

Tags in the code:
- TODO: Indicates areas requiring improvement.
- DEBUG: Marks functions used for current debugging.
*/



#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <GL/glut.h>   // GLUT library
#include <GL/glu.h>    // OpenGL (GLU)
#include <GL/gl.h>     // OpenGL
#include <GL/glext.h>

// #include "GltZpr/zpr.h"   // For zooming, translating, and rotating with the mouse
// Currently not working
// #include "GltZpr/zpr.c"   // For zooming, translating, and rotating with the mouse

#include <gsl/gsl_rng.h>       // For random number generator
#include <gsl/gsl_randist.h>   // For Gaussian random generator

#include "algorithmes.h"       // For binary search and waiting algorithms
#include "Internal_state.h"    // Internal state class
#include "constantes_SI.h"     // SI constants
#include "Kinetic_Monte_Carlo.h" // Kinetic Monte Carlo algorithm
#include "Laser.h"                   // Laser class
#include "Molecule.h"                // Molecule class
#include "Field.h"                   // Field class
#include "datacards.h"               // For reading the parameter file
#include "params.h"                  // For varying parameters

#include "Transitions_initialisation.h" // Initialization of transitions
#include "Initialisation_programme.h"   // General program initialization
#include "shift_molecule.h"          // For energy shifts (delta)
#include "diagonalization.h"         // For energy shifts (delta)
#include "Transition_rate_calcul.h"  // Calculates transition rates
#include "Stat_Molecule.h"           // For molecule statistics
#include "Affichage.h"               // Screen display
#include "Affichage_Mol.h"           // Molecule display
#include "sortie_donnees.h"          // Data output


using namespace std;


// Function called in case of memory overflow
void deborde()
{
    std::cerr << "Insufficient memory - execution stopped" << std::endl;
    exit(1);
}


/************************************************************************/
/************************** Main Program *******************************/
/******** main calls RePaint() in OPENGL + GLUT ***********************/
/************************************************************************/

// General function: initialization + KMC and N-body loop
void RePaint()
{
    /*****************************************************************/
    /** Creation of useful variables (parameters to be scanned or not) **/
    /*****************************************************************/

    std::string nomdat = "Data/Liste_Param.h";  // Name of the parameter file
    DataCards data(nomdat.c_str()); // Reads the file and creates the datacards.

    bool Graphics = (bool)data.IParam("Graphics"); // Enable or disable graphical output.
    const double SIZE_affichage = data.DParam("SIZE_affichage"); // Size of the display area
    MC_algorithmes Algorithme_MC = (MC_algorithmes)data.IParam("Choix_algorithme_Monte_Carlo"); // Choice of Monte Carlo algorithm
    N_Body_algorithmes Algorithme_N_body = (N_Body_algorithmes)data.IParam("Choix_algorithme_N_corps"); // Choice of N-body algorithm

    /*** FILE NAMES ***/
    std::string nom_sortie_temp_string, nom_sortie_scal_string, nom_file_Levels_string, nom_file_Lines_string, nom_sortie_donnees_string, nom_sortie_rate_string, nom_fichier_random_gen_string,
        nom_file_Laser_Spectrum_string, nom_file_Laser_Intensity_string, nom_file_Magn_Field_3D_string, nom_file_Elec_Field_3D_string;

    const char *nom_sortie_temp, *nom_sortie_scal, *nom_file_Levels, *nom_file_Lines, *nom_sortie_donnees, *nom_sortie_rate, *nom_fichier_random_gen,
          *nom_file_Laser_Spectrum, *nom_file_Laser_Intensity, *nom_file_Magn_Field_3D, *nom_file_Elec_Field_3D;

    nom_file_Levels_string = data.SParam("nom_file_Levels");      // File containing levels (state, energy, ...)
    nom_file_Lines_string = data.SParam("nom_file_Lines");        // File containing transitions
    nom_sortie_donnees_string = data.SParam("nom_sortie_donnees");
    nom_sortie_rate_string = data.SParam("nom_sortie_rate");
    nom_fichier_random_gen_string = data.SParam("nom_fichier_random_gen");
    nom_file_Laser_Spectrum_string = data.SParam("nom_file_Laser_Spectrum");   // File containing transitions
    nom_file_Laser_Intensity_string = data.SParam("nom_file_Laser_Intensity"); // File containing transitions
    nom_file_Magn_Field_3D_string = data.SParam("nom_file_Magn_Field_3D");
    nom_file_Elec_Field_3D_string = data.SParam("nom_file_Elec_Field_3D");

    nom_file_Levels = nom_file_Levels_string.c_str();      // File containing levels (state, energy, ...)
    nom_file_Lines = nom_file_Lines_string.c_str();        // File containing transitions
    nom_sortie_donnees = nom_sortie_donnees_string.c_str();
    nom_sortie_rate = nom_sortie_rate_string.c_str();
    nom_fichier_random_gen = nom_fichier_random_gen_string.c_str();
    nom_file_Laser_Spectrum = nom_file_Laser_Spectrum_string.c_str();   // File containing transitions
    nom_file_Laser_Intensity = nom_file_Laser_Intensity_string.c_str(); // File containing transitions
    nom_file_Magn_Field_3D = nom_file_Magn_Field_3D_string.c_str();
    nom_file_Elec_Field_3D = nom_file_Elec_Field_3D_string.c_str();

    FitParams params; // This is a vector<Param *>
    params.read(nomdat); // Reads the parameter file between "BEGIN_OF_FITPARAMS" and "END_OF_FITPARAMS".
    params.init_value(nomdat); // Initializes parameters based on the non-scanned value of the dataCard (cf datacards.h) or to the minimum value if scanned.

    /******************************************************************/
    /** Abstract creation of objects (fields, molecules, lasers, ...) **/
    /******************************************************************/

    std::vector<Internal_state> Level; // Level is a list of Internal_state
    std::vector<Laser> lasers;         // Lasers used for cooling

    Field champB, champE;              // Magnetic and electric fields

    const int Nb_type_of_Mol = data.IParam("Nb_type_of_Mol");
    int N_Mol[1];
    N_Mol[0] = (int)params.LocateParam("N_Mol[0]")->val; // The number of laser-cooled molecules
    std::vector<Molecule> Mol;

    std::vector<double> rate;                   // Rates (time-integrated) of transitions.
    std::vector<type_codage_react> reaction_list; // List of reactions according to the format given by type_codage_react, often (system number, final system state).

    std::ofstream file_out(nom_sortie_donnees); // Add ", ios::app" to write to the end of the file; otherwise, it erases the file on rewrite.
    file_out << std::setprecision(8);          // 8 decimals in output files
    std::ofstream file_rate(nom_sortie_rate);

    clock_t t_start, t_end; // For measuring program runtime
    t_start = clock();

    /******************************************************************/
    /********************** GSL Initialization ************************/
    /******************************************************************/

    gsl_rng *r; // Random number generator
    // set_random_generator(r, (int)params.LocateParam("Seed_Init_Random_Number_Generator")->val, nom_fichier_random_gen);  // Initialize the random generator (Mersenne Twister here)
    // TODO (Daniel#8#): Does not work (file for random generator). But everything is almost done, so it should be easy to fix.

    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;  // "Mersenne Twister" generator by default. Faster than RANLUX. See chapter 17.12 of gsl_ref.
    r = gsl_rng_alloc(T);

    initialisation_trans_mol(nom_file_Levels, nom_file_Lines, Level, params); // Reads level and transition files

    bool test_fin_boucle_param; // End-of-loop parameter for scanned parameters

    if (data.IParam("is_DataCard_out") == 1)
    {
        file_out << " Nb card " << data.NbCards() << std::endl << "DATA " << std::endl << data << std::endl;
    }
    // At the beginning, output all data from the initial data file
    if (data.IParam("is_DataCard_out") == 2)
    {
        std::ofstream file_out_datacard(data.SParam("nom_sortie_donnees_Data").c_str()); // Add ", ios::app" to write to the end of the file; otherwise, it erases the file on rewrite.
        file_out_datacard << " Nb card " << data.NbCards() << std::endl << "DATA " << std::endl << data << std::endl; // At the beginning, output all data from the initial data file
        file_out_datacard.close();
    }

    /*** Hamiltonian + dipole matrix element ***/
    MatrixXcd E0_cm;
    MatrixXcd Zeeman_cm_B;
    MatrixXcd H;        // Hamiltonian Matrix. It is a Hermitian matrix, so I use complex, not MatrixXd
    MatrixXd d0[3];
    MatrixXcd d[3];

    Read_Energy_Zeeman_dipole_for_Diagonalization(E0_cm, Zeeman_cm_B, d0); // Initialize energy, Zeeman, and dipole matrices

    /*****************************************************************/
    /** Global loop (over scanned parameters) ***********************/
    /*****************************************************************/

    do
    {
        double t = 0.0;                 // Start time of the simulation in µs
        double t_mise_a_jour = 0.0;     // Time to update potentials, acceleration, N-body, etc., for all molecules
        /** Potentials are only updated for the molecule whose internal state was modified, or when the N-body time step is reached, all are updated. **/
        double t_dia = 0.0;             // Time for file output
        double t_out = 0.0;             // Time for graphical output
        double t_fin = params.LocateParam("t_fin")->val; // Final time
        double dt_KMC;                  // Time for a reaction to occur
        double dt_dyn = 0.0;            // Dynamic time
        double dt_dia = params.LocateParam("dt_dia")->val; // Interval between diagnostics output (in `cout`)
        double dt_out = params.LocateParam("dt_out")->val; // Interval between graphical output snapshots (particle drawing)

        /*****************************************************************/
        /********************** Initialization of GSL ********************/
        /*****************************************************************/
        gsl_rng_set(r, (int)params.LocateParam("Seed_Init_Random_Number_Generator")->val);
        save_random_generator(r, nom_fichier_random_gen);

        /*****************************************************************/
        /******************** Initialization of molecules ***************/
        /*****************************************************************/

        Init_Field(champB, champE, params, nom_file_Magn_Field_3D, nom_file_Elec_Field_3D);  // Initialize fields
        // Create the list of molecules (vector `Mol`), its size will not be modified later
        Init_Molecule(r, Mol, champB, champE, Nb_type_of_Mol, params, data); // Position, velocity
        initialisation_proba(r, Mol, N_Mol[0], Level); // Initialize the populations of molecules at the start based on the desired population

        const int Nb_laser = data.IParam("Nb_laser");  // Number of lasers used (could have more in `Liste_Param`, but this is the number for this run)

        Init_Laser(lasers, Nb_laser, params, nom_file_Laser_Spectrum, nom_file_Laser_Intensity); // Initialize lasers
        // Sortie_laser_spectrum(file_out, lasers, params, 0); // Debug
        // Sortie_laser_intensity(file_out, lasers, params, 0);

        int number_mol = aucune; // Index of the molecule affected by a modification (none at the start, hence value -1)
        int number_photons = 0;  // Number of photons involved (absorption, spontaneous, or stimulated emission)

        /*****************************************************************/
        /** Global time loop *********************************************/
        /*****************************************************************/

        while (true) // Infinite loop until `t_end` is reached
        {
            Init_Field(champB, champE, params); // Re-initialize fields; important if scanning parameters change them
            Init_Laser(lasers, Nb_laser, params, nom_file_Laser_Spectrum, nom_file_Laser_Intensity); // Initialize lasers; could skip re-reading level files

            dt_dyn = params.LocateParam("dt_dyn_epsilon_param")->val;

            calcul_rates_molecules(Level, Algorithme_MC, reaction_list, rate, Mol, champB, champE, lasers, t, number_mol, N_Mol[0], params, H, E0_cm, Zeeman_cm_B, d0, d); // Compute transition rates for all molecules if `numero_mol == aucune`. Otherwise, recompute only for the molecule `numero_mol`

            if (t >= t_dia)
            {
                /*** Insert output logic here or just before for continuous file output ***/
                // Suggested examples: `Sortie_rate_example` and `Sortie_donnee_example` in `sortie_donnees.cpp`
                Sortie_donnee(file_out, Mol, Level, champB, champE, lasers, t, (int)Mol.size(), params, data, number_photons); // Output all molecular data
                t_dia += dt_dia;
            }

            int n_reac = find_reaction(Algorithme_MC, r, rate, dt_KMC); // Find the KMC reaction

            if (dt_KMC < dt_dyn || Algorithme_MC == Fast_Rough_Method)
            {
                evolve_step(Algorithme_N_body, Mol, champB, champE, lasers, Mol.size(), t, dt_KMC, params); // N-body evolution (before velocity update). Also modifies time.
                number_mol = do_reaction(Algorithme_MC, r, rate, reaction_list, Mol, n_reac, lasers, dt_KMC, file_rate, true, number_photons, params); // Perform the internal state evolution
            }
            else // Evolution will proceed without reaction for `dt_dyn`
            {
                evolve_step(Algorithme_N_body, Mol, champB, champE, lasers, Mol.size(), t, dt_dyn, params); // N-body evolution (before velocity update)
            }

// TODO (Daniel#6#): Maybe separate dt_dyn and t_mise_a_jour

            if (t >= t_mise_a_jour || Algorithme_MC == Fast_Rough_Method)
            {
                number_mol = aucune; // A complete update is needed (see rate calculation), meaning all molecules will be affected.
                // set_pot_all_mol(Mol, champB, champE, lasers, t, N_Mol[0], params);
                // Update all potentials (gravity, dipolar, magnetic, electric, etc.) with the new position to recalculate the correct transitions.
                t_mise_a_jour += dt_dyn;
            }

// One could consider accelerating the process by evolving the N-body system up to dt_KMC without recalculating dt_KMC each time.
// However, this is risky, as seen with a molecule far from the waist -> tKMC becomes immense but decreases quickly, and if we don't recalculate, we might make an error.


            if (Graphics && t >= t_out)
            {
                Draw(Mol, Level, champB, champE, lasers, SIZE_affichage, t, Mol.size(), Nb_type_of_Mol, params); // Display points
                t_out += dt_out;
                wait(params.LocateParam("t_wait_affichage")->val); // Prevent excessively fast rendering
            }

            if (t > t_fin) // End of loop
            {
                for (ParamIterator i = params.begin(); i != params.end(); ++i)
                {
                    Param &p = **i;
                    if (p.is_scanned == true)
                        std::cout << " " << p.val_t0;
                }
                std::cout << std::endl;
                break;
            }

            Modif_Param(t, params); // Dynamically modify laser wavelength, waist, etc., based on time
        }

        if (data.SParam("is_Scan_Random") == "true")
            test_fin_boucle_param = params.Scan_Param_aleatoire(r); // Modify parameters randomly
        else
            test_fin_boucle_param = params.Scan_Param();
    }
    while (!test_fin_boucle_param);   // End-of-loop check for parameters

    file_out.close();
    file_rate.close();
    t_end = clock();

    std::cout << "Execution time (s): " << (t_end - t_start) / double(CLOCKS_PER_SEC) << std::endl;
    std::cout << "Have you checked the to-do list, SMALL_NUMBER_RATE, and other parameters?" << std::endl;

    exit(1);
    system("PAUSE");
    return;
}



// fonction main call RePaint ().
// This is how it works in OPENGL + GLUT
int main(int argc, char** argv)
{


    string nomdat = "Data/Liste_Param.h" ;  // nom du fichier de paramètres
    DataCards data(nomdat.c_str()); // Lit le fichier et crée les datacards.
    bool Graphics = (bool) data.IParam("Graphics"); // affichage graphique ou non.

    if (Graphics)
    {
        const int size_screen = 600;
        glutInit(&argc, argv);                            // Initialisation de la GLUT
        glutInitDisplayMode(GLUT_DEPTH |                 // Active la profondeur
                            GLUT_DOUBLE |                // Double buffering
                            GLUT_RGBA);                  // Couleurs au format RGBA
        glutInitWindowPosition(0, 0);                    // Coordonnées de la fenêtre
        glutInitWindowSize(size_screen, size_screen);    // Taille de la fenêtre

        int windowID = glutCreateWindow("Cooling");      // Création de la fenêtre
        if (windowID == 0)
        {
            std::cerr << "Error: Unable to create window!" << std::endl;
            return -1;
        }

        // zprInit();                          //Pour zoomer, translater et tourner avec la souris. NE marche pas

        glutDisplayFunc(RePaint);                      //On lui dit d'appeler la fonction renderFunc pour afficher


        glutReshapeFunc(reshape);                      // Nécessaire pour recadrer la fonction
        glutMainLoop();                                   //Boucle infinie du programme
        exit(1);
        return 0;
    }
    else
        RePaint ();
    return 0;
    exit(1);
}
>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
