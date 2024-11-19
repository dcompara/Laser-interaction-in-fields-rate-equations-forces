<<<<<<< HEAD
/*
  Name:
  Author:
  Date: 23/10/06 10:01
  Description: Parameters for displaying molecules using OpenGL (see also Affichage.h).

  The main logic is in the function `affichage_one_bod`, which handles the display
  of a single particle.
*/

#ifndef AFFICHAGE_MOL_SEEN
#define AFFICHAGE_MOL_SEEN

#include <iostream>
using namespace std;

#include "constantes_SI.h"
#include "Vecteur3D.h"
#include "Molecule.h"
#include "Field.h"
#include "Affichage.h"

/************************************************************************/
/************************** Program Functions ***************************/
/************************************************************************/

/**
 * TODO: Add more flexible color customization.
 * Initializes the colors used for displaying particles.
 *
 * @param couleur Array of Vecteur3D to store colors.
 * @param Nmax Maximum number of states/colors to initialize.
 */
void init_couleurs(Vecteur3D *couleur, const int Nmax);

/**
 * TODO: Make sizes adjustable based on user-defined parameters.
 * Initializes the sizes used for displaying particles on screen.
 *
 * @param taille Array of doubles to store sizes.
 * @param Nmax Maximum number of sizes to initialize.
 */
void init_tailles(double *taille, const int Nmax);

/**
 * Displays a single particle.
 * The size is proportional to its vibration, and the rotation angles (J, MJ)
 * define its orientation along Ox and Oy axes.
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Particle to display.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_one_bod(const Body my_bod, FitParams &params);

/**
 * Displays the legend for particle states.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_state Number of distinct particle states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_legend(const vector <Body> &bod, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const double t, const int Nb_body, const int Nb_state, FitParams &params);

/**
 * Displays the color legend for particle states.
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Single particle.
 * @param Nb_state Number of states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_stat_color(Body my_bod, const int Nb_state, FitParams &params);

/**
 * Displays the particles within a defined zone.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Nb_body Total number of particles.
 * @param SIZE_real Scaling factor for size.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_zone(const vector<Body> &bod, const int Nb_body, const double SIZE_real, FitParams &params);

/**
 * Displays the entire system, including particles, states, and fields.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param fieldB Magnetic field.
 * @param fieldE Electric field.
 * @param lasers Vector of laser beams.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_of_different_type_body Number of distinct types of particles.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void Draw(const vector<Body> &bod, const vector<Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector<Laser> &lasers, const double t, const int Nb_body, const int Nb_of_different_type_body, FitParams &params);

/**************************************
Note: Templates must be implemented in the header file.
If a template is used in multiple files, it cannot be separated into .h and .cpp.
For more details: http://www.commentcamarche.net/faq/11194-les-templates-en-c
***************************************/




/************************************************************************/
/********  PARAMETERS TO DRAW MOLECULES (colors, sizes, etc.) ***********/
/************************************************************************/

/**
 * Basic colors in RGB format, normalized for OpenGL (values divided by 255).
 */
const Vecteur3D Black = Vecteur3D(0, 0, 0) / 255;
const Vecteur3D White = Vecteur3D(255, 255, 255) / 255;
const Vecteur3D Red = Vecteur3D(255, 0, 0) / 255;
const Vecteur3D Lime = Vecteur3D(0, 255, 0) / 255;
const Vecteur3D Blue = Vecteur3D(0, 0, 255) / 255;
const Vecteur3D Yellow = Vecteur3D(255, 255, 0) / 255;
const Vecteur3D Cyan_Aqua = Vecteur3D(0, 255, 255) / 255;
const Vecteur3D Magenta_Fuchsia = Vecteur3D(255, 0, 255) / 255;
const Vecteur3D Silver = Vecteur3D(192, 192, 192) / 255;
const Vecteur3D Gray = Vecteur3D(128, 128, 128) / 255;
const Vecteur3D Maroon = Vecteur3D(128, 0, 0) / 255;
const Vecteur3D Olive = Vecteur3D(128, 128, 0) / 255;
const Vecteur3D Green = Vecteur3D(0, 128, 0) / 255;
const Vecteur3D Purple = Vecteur3D(128, 0, 128) / 255;
const Vecteur3D Teal = Vecteur3D(0, 128, 128) / 255;
const Vecteur3D Navy = Vecteur3D(0, 0, 128) / 255;

/**
 * Class Param_Draw_Mol: Handles parameters for drawing molecules.
 * TODO: Implement this class if advanced molecule rendering is required.
 */
/*
class Param_Draw_Mol {
protected:
    string name;
    double mass;
    map<int, Vecteur3D> color_electronical_state;
    map<int, Vecteur3D> size_electronical_state;
    double scale_vibrational_state; // v
    double angle_rotational_state;  // J
    double angle_angular_state;     // M

public:
    Param_Draw_Mol();              // Constructor
    virtual ~Param_Draw_Mol();     // Destructor

    virtual void read(istream &flux) {
        flux >> name;
        flux >> mass;
        flux >> scale_vibrational_state;
        flux >> angle_rotational_state;
        flux >> angle_angular_state;
    }
};

istream& operator>>(istream &flux, Param_Draw_Mol &at) {
    at.read(flux);
    return flux;
}

// Default Constructor
Param_Draw_Mol::Param_Draw_Mol() {
    name = "Cs";
    mass = MCs;
};
*/

/**
 * Initializes colors for particle states.
 *
 * - Green: Ground state.
 * - Yellow: Excited state.
 * - Blue: Second excited state.
 * - Olive: Third excited state.
 * - Others: A gradient based on state index.
 */
void init_couleurs(Vecteur3D *couleur, const int Nmax) {
    couleur[fond] = Green;
    couleur[excite] = Yellow;
    couleur[excite2] = Blue;
    couleur[excite3] = Olive;
    for (int i = excite3; i < Nmax; i++) {
        couleur[i] = Vecteur3D(50 + 2 * i, 50 + 2 * i, 0) / 255;
    }
}

/**
 * Initializes sizes for particle states.
 *
 * - Size is proportional to the excitation state.
 */
void init_tailles(double *taille, const int Nmax) {
    taille[fond] = 0.01;
    taille[excite] = 0.015;
    taille[excite2] = 0.02;
    taille[excite3] = 0.03;
    for (int i = excite3; i < Nmax; i++) {
        taille[i] = 0.005 * i;
    }
}

/**
 * Displays a single particle.
 *
 * - Size is proportional to vibration.
 * - Orientation is defined by rotation angles J (around Ox) and MJ (around Oy).
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Particle to display.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_one_bod(const Body my_bod, FitParams &params) {
    int etat = my_bod.exc; // Excitation state
    if (etat < 0) etat = excite3; // Handle lost/ionized particles
    if (my_bod.get_name() == "P_bar") etat = excite3; // Special handling for "P_bar"

    int v = my_bod.v;
    int two_J = my_bod.two_J;
    int two_M = my_bod.two_M;

    const int Nb_max = 100; // Conservative limit
    if (etat >= Nb_max) etat = excite2; // Avoid array out-of-bounds

    Vecteur3D couleur_bod[Nb_max];
    double tailles_bod[Nb_max];
    init_couleurs(couleur_bod, Nb_max);
    init_tailles(tailles_bod, Nb_max);

    int two_J_max = 2 * params.LocateParam("two_J_maX")->val; // Maximum J for display
    int v_max = params.LocateParam("v_maX")->val; // Maximum vibration state

    // Rotate based on particle state
    glRotated(90. * two_J / (two_J_max + 1.1), 0., 0., 1.);
    glRotated(90. * two_M / (two_J + 1.1), 0., 1., 0.);

    double DIAM = (1.0 / 50.0);
    double length = 0.2 * (v + 1.2) / (v_max + 1.1);

    glColor3d(couleur_bod[etat].x(), couleur_bod[etat].y(), couleur_bod[etat].z());
    glutSolidSphere(tailles_bod[etat], 10, 10); // Draw sphere

    glPushMatrix();
    GLUquadricObj *quad = gluNewQuadric();
    glRotatef(90., 0.0, 1.0, 0.0); // Necessary adjustment

    gluCylinder(quad, length * DIAM, length * DIAM, 0.8 * length, 10, 10);
    glPopMatrix();
    gluDeleteQuadric(quad);

    glTranslatef(length, 0., 0.); // Position of second atom
    glColor3d(1, 1, 1);
    glutSolidSphere(tailles_bod[etat], 10, 10);
}

/**
 * Displays the legend for particle states.
 *
 * - Includes details like state population, velocity, and temperature.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_state Number of distinct particle states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_legend(const vector <Body> &bod, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const double t, const int Nb_body, const int Nb_state, FitParams &params)
{

    int Nb_max=100;
    Vecteur3D couleur_bod [Nb_max]; //  For safety because the last color is excite3 3;
    double tailles_bod [Nb_max];

    init_couleurs(couleur_bod,Nb_max);
    init_tailles(tailles_bod,Nb_max);


    int size_screen =600;     // Cela correspond à la taille dans main
    glViewport(0,0,size_screen,size_screen);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,size_screen,0,size_screen);
    glScalef(1,-1,1);
    glTranslatef(0,-size_screen,0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    Body bod_avg, bod_var;
    // vector <Stat> stat_Mol;
    // stat_molecule(bod, stat_Mol, Level,  fieldB, fieldE, Nb_body, Nb_state, params);


    // Affiche le nombre d'molecules excités dans chaque état exc
    int taille_police = 18;
    int position_y =70;
    int step_y = 1;


    position_y -= taille_police+2;



    int   i=1;
    glColor3d ( couleur_bod[i].x(),couleur_bod[i].y(),couleur_bod[i].z());
    bitmap_output(25,position_y,const_cast<char *>(" excité "),GLUT_BITMAP_TIMES_ROMAN_24); // Le const_cast<char *> est là pour éviter le warning: deprecated conversion from string constant to 'char*'
    position_y -= taille_police+step_y;
    i=0;
    glColor3d (couleur_bod[i].x(),couleur_bod[i].y(),couleur_bod[i].z());
    bitmap_output(25,position_y,const_cast<char *>(" fond "),GLUT_BITMAP_TIMES_ROMAN_24);

    Stat stat_Mol;
    stat_molecule_un_niveau(bod, stat_Mol, Level,  fieldB, fieldE, all_levels , Nb_body, params);

    int nb_fond = stat_Mol.population;

    glColor3d(1, 1, 1);    // Blanc
    position_y += taille_police+360;

    texte_output(t/nanosecond," t (ns) = ",5, position_y);
    position_y += taille_police+step_y;

    texte_output(stat_Mol.mean_v.x()/cm," vx (cm/s) = ",5, position_y);
    texte_output(stat_Mol.Temp_1D_50.x()/muK,"                                     Tx (mu K) = ",5, position_y); // Calcul la température kB T=m v^2 sur chaque axe

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_v.y()/cm," vy (cm/s) = ",5, position_y);
    texte_output(stat_Mol.Temp_1D_50.y()/muK,"                                     Ty (mu K) = ",5, position_y); // On suppose que toutes les particules ont la même masse

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_v.z()/cm," vz (cm/s) = ",5, position_y);
    texte_output(stat_Mol.Temp_1D_50.z()/muK,"                                      Tz (mu K) = ",5, position_y);

    position_y += taille_police+step_y;
    texte_output((stat_Mol.E_pot/kB)/mK/1.5/nb_fond," E_pot(mK) ",5, position_y);
    texte_output((stat_Mol.E_cin/kB)/mK/1.5/nb_fond,"                              E_cin(mK)  ",5, position_y);

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_pos.x()/MICRON," x (mu m) = ",5, position_y);
    texte_output(stat_Mol.sigma_pos.x()/MICRON,"                                         Dx (mu m) = ",5, position_y);

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_pos.y()/MICRON," y (mu m) = ",5, position_y);
    texte_output(stat_Mol.sigma_pos.y()/MICRON,"                                         Dy (mu m) = ",5, position_y);

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_pos.z()/MICRON," z (mu m) = ",5, position_y);
    texte_output(stat_Mol.sigma_pos.z()/MICRON,"                                         Dz (mu m) = ",5, position_y);

    /******  Stats of all molecules, of N_Mol[0] and others ******/
    position_y += taille_police+step_y;
    texte_output(((stat_Mol.E_pot+stat_Mol.E_cin)/kB)/mK/3./nb_fond," E_tot(mK) ",5, position_y);
    texte_output(stat_Mol.Temp_3D_50/mK,"                               T_all (mK)= ",3, position_y);
// Attention relative à la température E_pot = 3/2 k T, E_cin =3/2 kT; E tot=3 kT

    int N_Mol_0 = params.LocateParam("N_Mol[0]")->val ;
    vector <Molecule> liste_Adresses_Mol_dans_niveau; //List of Molecule of first type
    for (int i=0; i!= N_Mol_0; i++)
        liste_Adresses_Mol_dans_niveau.push_back(bod[i]);
    stat_molecule_form_list(bod, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  Nb_body, params);
    texte_output(stat_Mol.Temp_3D_50/mK,"                                                                T[0]= ",3, position_y);

    liste_Adresses_Mol_dans_niveau.clear();
    for (int i=N_Mol_0; i!= Nb_body; i++)
        liste_Adresses_Mol_dans_niveau.push_back(bod[i]);
    stat_molecule_form_list(bod, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  Nb_body, params);
    texte_output(stat_Mol.Temp_3D_50/mK,"                                                                                    T[1]= ",3, position_y);



}


/**
 * Displays color-coded statistics for particles.
 *
 * - Assigns positions to particles based on state index.
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Single particle.
 * @param Nb_state Number of states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_stat_color(Body my_bod, const int Nb_state, FitParams &params) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW); // Select transformation stack for "point of view"
    glLoadIdentity(); // Reset transformation stack

    glTranslatef(-0.95, 0.95, 0);

    for (int i = 0; i != Nb_state; i++) {
        glPushMatrix(); // Apply transformations only to subsequent elements
        my_bod.set_pos(Vecteur3D(0, -(i + 2) / 10., 0));

        affichage_one_bod(my_bod, params);
        glPopMatrix(); // Restore previous transformations
    }
}

/**
 * Displays particles in a defined region.
 *
 * - Adjusts positions to the given scaling factor.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Nb_body Total number of particles.
 * @param SIZE_real Scaling factor for size.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_zone(const vector<Body> &bod, const int Nb_body, const double SIZE_real, FitParams &params) {
    for (int i = 0; i < Nb_body; i++) {
        Vecteur3D position(bod[i].get_pos() / SIZE_real); // Scale position

        glPushMatrix();
        glTranslatef(position.x(), position.y(), -position.z());
        affichage_one_bod(bod[i], params);
        glPopMatrix();
    }
}

/**
 * Main function to draw the entire system.
 *
 * - Includes particles, states, fields, and the coordinate system.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param fieldB Magnetic field.
 * @param fieldE Electric field.
 * @param lasers Vector of laser beams.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_of_different_type_body Number of distinct types of particles.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void Draw(const vector<Body> &bod, const vector<Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector<Laser> &lasers, const double SIZE_real, const double t, const int Nb_body, const int Nb_of_different_type_body, FitParams &params) {
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);

    int Nb_state = Level.size();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-1.0, 1.0, -1.0, 1.0, 1, 3); // Define viewing volume

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, 0, -Z0);

    double angle = 90.0; // Default angle for gravity alignment
    double rot_axe_x = params.LocateParam("rot_axe_x")->val;
    double rot_axe_y = params.LocateParam("rot_axe_y")->val;
    double rot_axe_z = params.LocateParam("rot_axe_z")->val;
    glRotated(angle, rot_axe_x, rot_axe_y, rot_axe_z);

    affichage_zone(bod, Nb_body, SIZE_real, params); // Draw particles
    traceRepere(0.5); // Draw coordinate system
    affichage_legend(bod, Level,  fieldB, fieldE, t, Nb_body,Nb_state, params) ;     // Affiche les statistiques sur les molecules

    glutSwapBuffers(); // Avoid flickering during rendering
}

#endif

=======
/*
  Name:
  Author:
  Date: 23/10/06 10:01
  Description: Parameters for displaying molecules using OpenGL (see also Affichage.h).

  The main logic is in the function `affichage_one_bod`, which handles the display
  of a single particle.
*/

#ifndef AFFICHAGE_MOL_SEEN
#define AFFICHAGE_MOL_SEEN

#include <iostream>
using namespace std;

#include "constantes_SI.h"
#include "Vecteur3D.h"
#include "Molecule.h"
#include "Field.h"
#include "Affichage.h"

/************************************************************************/
/************************** Program Functions ***************************/
/************************************************************************/

/**
 * TODO: Add more flexible color customization.
 * Initializes the colors used for displaying particles.
 *
 * @param couleur Array of Vecteur3D to store colors.
 * @param Nmax Maximum number of states/colors to initialize.
 */
void init_couleurs(Vecteur3D *couleur, const int Nmax);

/**
 * TODO: Make sizes adjustable based on user-defined parameters.
 * Initializes the sizes used for displaying particles on screen.
 *
 * @param taille Array of doubles to store sizes.
 * @param Nmax Maximum number of sizes to initialize.
 */
void init_tailles(double *taille, const int Nmax);

/**
 * Displays a single particle.
 * The size is proportional to its vibration, and the rotation angles (J, MJ)
 * define its orientation along Ox and Oy axes.
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Particle to display.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_one_bod(const Body my_bod, FitParams &params);

/**
 * Displays the legend for particle states.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_state Number of distinct particle states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_legend(const vector <Body> &bod, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const double t, const int Nb_body, const int Nb_state, FitParams &params);

/**
 * Displays the color legend for particle states.
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Single particle.
 * @param Nb_state Number of states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_stat_color(Body my_bod, const int Nb_state, FitParams &params);

/**
 * Displays the particles within a defined zone.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Nb_body Total number of particles.
 * @param SIZE_real Scaling factor for size.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_zone(const vector<Body> &bod, const int Nb_body, const double SIZE_real, FitParams &params);

/**
 * Displays the entire system, including particles, states, and fields.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param fieldB Magnetic field.
 * @param fieldE Electric field.
 * @param lasers Vector of laser beams.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_of_different_type_body Number of distinct types of particles.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void Draw(const vector<Body> &bod, const vector<Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector<Laser> &lasers, const double t, const int Nb_body, const int Nb_of_different_type_body, FitParams &params);

/**************************************
Note: Templates must be implemented in the header file.
If a template is used in multiple files, it cannot be separated into .h and .cpp.
For more details: http://www.commentcamarche.net/faq/11194-les-templates-en-c
***************************************/




/************************************************************************/
/********  PARAMETERS TO DRAW MOLECULES (colors, sizes, etc.) ***********/
/************************************************************************/

/**
 * Basic colors in RGB format, normalized for OpenGL (values divided by 255).
 */
const Vecteur3D Black = Vecteur3D(0, 0, 0) / 255;
const Vecteur3D White = Vecteur3D(255, 255, 255) / 255;
const Vecteur3D Red = Vecteur3D(255, 0, 0) / 255;
const Vecteur3D Lime = Vecteur3D(0, 255, 0) / 255;
const Vecteur3D Blue = Vecteur3D(0, 0, 255) / 255;
const Vecteur3D Yellow = Vecteur3D(255, 255, 0) / 255;
const Vecteur3D Cyan_Aqua = Vecteur3D(0, 255, 255) / 255;
const Vecteur3D Magenta_Fuchsia = Vecteur3D(255, 0, 255) / 255;
const Vecteur3D Silver = Vecteur3D(192, 192, 192) / 255;
const Vecteur3D Gray = Vecteur3D(128, 128, 128) / 255;
const Vecteur3D Maroon = Vecteur3D(128, 0, 0) / 255;
const Vecteur3D Olive = Vecteur3D(128, 128, 0) / 255;
const Vecteur3D Green = Vecteur3D(0, 128, 0) / 255;
const Vecteur3D Purple = Vecteur3D(128, 0, 128) / 255;
const Vecteur3D Teal = Vecteur3D(0, 128, 128) / 255;
const Vecteur3D Navy = Vecteur3D(0, 0, 128) / 255;

/**
 * Class Param_Draw_Mol: Handles parameters for drawing molecules.
 * TODO: Implement this class if advanced molecule rendering is required.
 */
/*
class Param_Draw_Mol {
protected:
    string name;
    double mass;
    map<int, Vecteur3D> color_electronical_state;
    map<int, Vecteur3D> size_electronical_state;
    double scale_vibrational_state; // v
    double angle_rotational_state;  // J
    double angle_angular_state;     // M

public:
    Param_Draw_Mol();              // Constructor
    virtual ~Param_Draw_Mol();     // Destructor

    virtual void read(istream &flux) {
        flux >> name;
        flux >> mass;
        flux >> scale_vibrational_state;
        flux >> angle_rotational_state;
        flux >> angle_angular_state;
    }
};

istream& operator>>(istream &flux, Param_Draw_Mol &at) {
    at.read(flux);
    return flux;
}

// Default Constructor
Param_Draw_Mol::Param_Draw_Mol() {
    name = "Cs";
    mass = MCs;
};
*/

/**
 * Initializes colors for particle states.
 *
 * - Green: Ground state.
 * - Yellow: Excited state.
 * - Blue: Second excited state.
 * - Olive: Third excited state.
 * - Others: A gradient based on state index.
 */
void init_couleurs(Vecteur3D *couleur, const int Nmax) {
    couleur[fond] = Green;
    couleur[excite] = Yellow;
    couleur[excite2] = Blue;
    couleur[excite3] = Olive;
    for (int i = excite3; i < Nmax; i++) {
        couleur[i] = Vecteur3D(50 + 2 * i, 50 + 2 * i, 0) / 255;
    }
}

/**
 * Initializes sizes for particle states.
 *
 * - Size is proportional to the excitation state.
 */
void init_tailles(double *taille, const int Nmax) {
    taille[fond] = 0.01;
    taille[excite] = 0.015;
    taille[excite2] = 0.02;
    taille[excite3] = 0.03;
    for (int i = excite3; i < Nmax; i++) {
        taille[i] = 0.005 * i;
    }
}

/**
 * Displays a single particle.
 *
 * - Size is proportional to vibration.
 * - Orientation is defined by rotation angles J (around Ox) and MJ (around Oy).
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Particle to display.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_one_bod(const Body my_bod, FitParams &params) {
    int etat = my_bod.exc; // Excitation state
    if (etat < 0) etat = excite3; // Handle lost/ionized particles
    if (my_bod.get_name() == "P_bar") etat = excite3; // Special handling for "P_bar"

    int v = my_bod.v;
    int two_J = my_bod.two_J;
    int two_M = my_bod.two_M;

    const int Nb_max = 100; // Conservative limit
    if (etat >= Nb_max) etat = excite2; // Avoid array out-of-bounds

    Vecteur3D couleur_bod[Nb_max];
    double tailles_bod[Nb_max];
    init_couleurs(couleur_bod, Nb_max);
    init_tailles(tailles_bod, Nb_max);

    int two_J_max = 2 * params.LocateParam("two_J_maX")->val; // Maximum J for display
    int v_max = params.LocateParam("v_maX")->val; // Maximum vibration state

    // Rotate based on particle state
    glRotated(90. * two_J / (two_J_max + 1.1), 0., 0., 1.);
    glRotated(90. * two_M / (two_J + 1.1), 0., 1., 0.);

    double DIAM = (1.0 / 50.0);
    double length = 0.2 * (v + 1.2) / (v_max + 1.1);

    glColor3d(couleur_bod[etat].x(), couleur_bod[etat].y(), couleur_bod[etat].z());
    glutSolidSphere(tailles_bod[etat], 10, 10); // Draw sphere

    glPushMatrix();
    GLUquadricObj *quad = gluNewQuadric();
    glRotatef(90., 0.0, 1.0, 0.0); // Necessary adjustment

    gluCylinder(quad, length * DIAM, length * DIAM, 0.8 * length, 10, 10);
    glPopMatrix();
    gluDeleteQuadric(quad);

    glTranslatef(length, 0., 0.); // Position of second atom
    glColor3d(1, 1, 1);
    glutSolidSphere(tailles_bod[etat], 10, 10);
}

/**
 * Displays the legend for particle states.
 *
 * - Includes details like state population, velocity, and temperature.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_state Number of distinct particle states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_legend(const vector <Body> &bod, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const double t, const int Nb_body, const int Nb_state, FitParams &params)
{

    int Nb_max=100;
    Vecteur3D couleur_bod [Nb_max]; //  For safety because the last color is excite3 3;
    double tailles_bod [Nb_max];

    init_couleurs(couleur_bod,Nb_max);
    init_tailles(tailles_bod,Nb_max);


    int size_screen =600;     // Cela correspond à la taille dans main
    glViewport(0,0,size_screen,size_screen);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,size_screen,0,size_screen);
    glScalef(1,-1,1);
    glTranslatef(0,-size_screen,0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    Body bod_avg, bod_var;
    // vector <Stat> stat_Mol;
    // stat_molecule(bod, stat_Mol, Level,  fieldB, fieldE, Nb_body, Nb_state, params);


    // Affiche le nombre d'molecules excités dans chaque état exc
    int taille_police = 18;
    int position_y =70;
    int step_y = 1;


    position_y -= taille_police+2;



    int   i=1;
    glColor3d ( couleur_bod[i].x(),couleur_bod[i].y(),couleur_bod[i].z());
    bitmap_output(25,position_y,const_cast<char *>(" excité "),GLUT_BITMAP_TIMES_ROMAN_24); // Le const_cast<char *> est là pour éviter le warning: deprecated conversion from string constant to 'char*'
    position_y -= taille_police+step_y;
    i=0;
    glColor3d (couleur_bod[i].x(),couleur_bod[i].y(),couleur_bod[i].z());
    bitmap_output(25,position_y,const_cast<char *>(" fond "),GLUT_BITMAP_TIMES_ROMAN_24);

    Stat stat_Mol;
    stat_molecule_un_niveau(bod, stat_Mol, Level,  fieldB, fieldE, all_levels , Nb_body, params);

    int nb_fond = stat_Mol.population;

    glColor3d(1, 1, 1);    // Blanc
    position_y += taille_police+360;

    texte_output(t/nanosecond," t (ns) = ",5, position_y);
    position_y += taille_police+step_y;

    texte_output(stat_Mol.mean_v.x()/cm," vx (cm/s) = ",5, position_y);
    texte_output(stat_Mol.Temp_1D_50.x()/muK,"                                     Tx (mu K) = ",5, position_y); // Calcul la température kB T=m v^2 sur chaque axe

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_v.y()/cm," vy (cm/s) = ",5, position_y);
    texte_output(stat_Mol.Temp_1D_50.y()/muK,"                                     Ty (mu K) = ",5, position_y); // On suppose que toutes les particules ont la même masse

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_v.z()/cm," vz (cm/s) = ",5, position_y);
    texte_output(stat_Mol.Temp_1D_50.z()/muK,"                                      Tz (mu K) = ",5, position_y);

    position_y += taille_police+step_y;
    texte_output((stat_Mol.E_pot/kB)/mK/1.5/nb_fond," E_pot(mK) ",5, position_y);
    texte_output((stat_Mol.E_cin/kB)/mK/1.5/nb_fond,"                              E_cin(mK)  ",5, position_y);

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_pos.x()/MICRON," x (mu m) = ",5, position_y);
    texte_output(stat_Mol.sigma_pos.x()/MICRON,"                                         Dx (mu m) = ",5, position_y);

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_pos.y()/MICRON," y (mu m) = ",5, position_y);
    texte_output(stat_Mol.sigma_pos.y()/MICRON,"                                         Dy (mu m) = ",5, position_y);

    position_y += taille_police+step_y;
    texte_output(stat_Mol.mean_pos.z()/MICRON," z (mu m) = ",5, position_y);
    texte_output(stat_Mol.sigma_pos.z()/MICRON,"                                         Dz (mu m) = ",5, position_y);

    /******  Stats of all molecules, of N_Mol[0] and others ******/
    position_y += taille_police+step_y;
    texte_output(((stat_Mol.E_pot+stat_Mol.E_cin)/kB)/mK/3./nb_fond," E_tot(mK) ",5, position_y);
    texte_output(stat_Mol.Temp_3D_50/mK,"                               T_all (mK)= ",3, position_y);
// Attention relative à la température E_pot = 3/2 k T, E_cin =3/2 kT; E tot=3 kT

    int N_Mol_0 = params.LocateParam("N_Mol[0]")->val ;
    vector <Molecule> liste_Adresses_Mol_dans_niveau; //List of Molecule of first type
    for (int i=0; i!= N_Mol_0; i++)
        liste_Adresses_Mol_dans_niveau.push_back(bod[i]);
    stat_molecule_form_list(bod, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  Nb_body, params);
    texte_output(stat_Mol.Temp_3D_50/mK,"                                                                T[0]= ",3, position_y);

    liste_Adresses_Mol_dans_niveau.clear();
    for (int i=N_Mol_0; i!= Nb_body; i++)
        liste_Adresses_Mol_dans_niveau.push_back(bod[i]);
    stat_molecule_form_list(bod, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  Nb_body, params);
    texte_output(stat_Mol.Temp_3D_50/mK,"                                                                                    T[1]= ",3, position_y);



}


/**
 * Displays color-coded statistics for particles.
 *
 * - Assigns positions to particles based on state index.
 *
 * @tparam Body Class representing a particle.
 * @param my_bod Single particle.
 * @param Nb_state Number of states.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_stat_color(Body my_bod, const int Nb_state, FitParams &params) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW); // Select transformation stack for "point of view"
    glLoadIdentity(); // Reset transformation stack

    glTranslatef(-0.95, 0.95, 0);

    for (int i = 0; i != Nb_state; i++) {
        glPushMatrix(); // Apply transformations only to subsequent elements
        my_bod.set_pos(Vecteur3D(0, -(i + 2) / 10., 0));

        affichage_one_bod(my_bod, params);
        glPopMatrix(); // Restore previous transformations
    }
}

/**
 * Displays particles in a defined region.
 *
 * - Adjusts positions to the given scaling factor.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Nb_body Total number of particles.
 * @param SIZE_real Scaling factor for size.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void affichage_zone(const vector<Body> &bod, const int Nb_body, const double SIZE_real, FitParams &params) {
    for (int i = 0; i < Nb_body; i++) {
        Vecteur3D position(bod[i].get_pos() / SIZE_real); // Scale position

        glPushMatrix();
        glTranslatef(position.x(), position.y(), -position.z());
        affichage_one_bod(bod[i], params);
        glPopMatrix();
    }
}

/**
 * Main function to draw the entire system.
 *
 * - Includes particles, states, fields, and the coordinate system.
 *
 * @tparam Body Class representing a particle.
 * @param bod Vector of particles.
 * @param Level Vector of internal states.
 * @param fieldB Magnetic field.
 * @param fieldE Electric field.
 * @param lasers Vector of laser beams.
 * @param t Current time.
 * @param Nb_body Total number of particles.
 * @param Nb_of_different_type_body Number of distinct types of particles.
 * @param params FitParams containing relevant parameters for display.
 */
template <class Body>
void Draw(const vector<Body> &bod, const vector<Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector<Laser> &lasers, const double SIZE_real, const double t, const int Nb_body, const int Nb_of_different_type_body, FitParams &params) {
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_ACCUM_BUFFER_BIT);

    int Nb_state = Level.size();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-1.0, 1.0, -1.0, 1.0, 1, 3); // Define viewing volume

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, 0, -Z0);

    double angle = 90.0; // Default angle for gravity alignment
    double rot_axe_x = params.LocateParam("rot_axe_x")->val;
    double rot_axe_y = params.LocateParam("rot_axe_y")->val;
    double rot_axe_z = params.LocateParam("rot_axe_z")->val;
    glRotated(angle, rot_axe_x, rot_axe_y, rot_axe_z);

    affichage_zone(bod, Nb_body, SIZE_real, params); // Draw particles
    traceRepere(0.5); // Draw coordinate system
    affichage_legend(bod, Level,  fieldB, fieldE, t, Nb_body,Nb_state, params) ;     // Affiche les statistiques sur les molecules

    glutSwapBuffers(); // Avoid flickering during rendering
}

#endif

>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
