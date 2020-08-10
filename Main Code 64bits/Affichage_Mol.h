/*
  Name:
  Copyright:
  Author:
  Date: 23/10/06 10:01
  Description: Paramètre pour l'affichage des molécules avec OPEN_GL (voir aussi Affichage.h)

  A priori tout est dans le dessin d'une particule unique
affichage_one_bod

*/



#ifndef Affichage_Mol_SEEN
#define Affichage_Mol_SEEN

#include <iostream>
using namespace std;

#include "constantes_SI.h"
#include "Vecteur3D.h"
#include "Molecule.h"
#include "Field.h"
#include "Affichage.h"

/************************************************************************/
/************************** Programme *************************/
/************************************************************************/


// Choix des couleurs

//Vecteur3D(0,1,0);               // vert
//Vecteur3D(1,1,0);               // jaune
//Vecteur3D(1,0,0);               // rouge
//Vecteur3D(0,0,1);               // bleu
//Vecteur3D(1,1,1);               // blanc
//Vecteur3D(1,0,1);               // violet



void init_couleurs(Vecteur3D *couleur, const int Nmax);

// initialisation des tailles affichées sur l'écran
void init_tailles(double *taille, const int Nmax);

// dessin d'une particule unique
// Pour nous la taille est proportionnelle à la vibration
// L'angle de rotation est donné par la rotation J (selon Ox), MJ (selon 0y)
template <class Body>
void affichage_one_bod(const Body my_bod, FitParams &params);  // Affichage des points

template <class Body>
void affichage_legend(const vector <Body> &bod, const vector <Internal_state> &Level, const double t, const int Nb_body, const int Nb_state, FitParams &params);

// Affiche la légende pour l'état des particules
template <class Body>
void affichage_stat_color(Body my_bod, const int Nb_state, FitParams &params);

// Affichage des particules
template <class Body>
void affichage_zone(const vector <Body> &bod, const int Nb_body,  const double SIZE_real, FitParams &params);  // Affichage des points

// Affichage du système
template <class Body>
void Draw(const vector <Body> &bod, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const vector <Laser> &lasers, const double t, const int Nb_body, const int Nb_of_different_type_body, FitParams &params);




/**************************************
JE MET ICI L'equivalent du .cpp car "On peut pas séparer un .h et .cpp pour les template".

En effet si un symbole template est utilisé uniquement dans un .cpp (fichier source), il peut être implémenté dans ce .cpp
Sinon, il doit être implémenté dans un .hpp (header).

cf http://www.commentcamarche.net/faq/11194-les-templates-en-c
***************************************/


/************************************************************************/
/********  PARAMETER TO DRAW MOLECULES (color, size, ...) ***************/
/************************************************************************/


// List of Basic color in RGB code, but divided by 255 to fit for OPENGL
const	Vecteur3D	Black	=	Vecteur3D	(0,0,0)	/	255	;
const	Vecteur3D	White	=	Vecteur3D	(255,255,255)	/	255	;
const	Vecteur3D	Red	=	Vecteur3D	(255,0,0)	/	255	;
const	Vecteur3D	Lime	=	Vecteur3D	(0,255,0)	/	255	;
const	Vecteur3D	Blue	=	Vecteur3D	(0,0,255)	/	255	;
const	Vecteur3D	Yellow	=	Vecteur3D	(255,255,0)	/	255	;
const	Vecteur3D	Cyan_Aqua	=	Vecteur3D	(0,255,255)	/	255	;
const	Vecteur3D	Magenta_Fuchsia	=	Vecteur3D	(255,0,255)	/	255	;
const	Vecteur3D	Silver	=	Vecteur3D	(192,192,192)	/	255	;
const	Vecteur3D	Gray	=	Vecteur3D	(128,128,128)	/	255	;
const	Vecteur3D	Maroon	=	Vecteur3D	(128,0,0)	/	255	;
const	Vecteur3D	Olive	=	Vecteur3D	(128,128,0)	/	255	;
const	Vecteur3D	Green	=	Vecteur3D	(0,128,0)	/	255	;
const	Vecteur3D	Purple	=	Vecteur3D	(128,0,128)	/	255	;
const	Vecteur3D	Teal	=	Vecteur3D	(0,128,128)	/	255	;
const	Vecteur3D	Navy	=	Vecteur3D	(0,0,128)	/	255	;


/*** Class Param_Draw_Mol to decide how to draw a molecule ***/

/***********

class Param_Draw_Mol
{
protected:
    string name;
    double mass;
    map <int,Vecteur3D> color_electronical_state;
    map <int,Vecteur3D> size_electronical_state;
    double scale_vibrational_state; // v
    double angle_rotational_state;  // J
    double angle_angular_state;     // M

public:
    Param_Draw_Mol();                        //constructor
    virtual ~Param_Draw_Mol();               // Destructor

    virtual void read(istream & flux)
    {
        flux >> name;
        flux >> mass;
      //  flux >> color_electronical_state;
      //  flux >> size_electronical_state;
        flux >> scale_vibrational_state;
        flux >> angle_rotational_state;
        flux >> angle_angular_state;
    }
};

istream& operator >> (istream &flux,Param_Draw_Mol & at)
{
    at.read(flux);
    return(flux);
}

// Default Constructor
Param_Draw_Mol::Param_Draw_Mol()
{
    name = "Cs";
    mass =  MCs;
 //   color_electronical_state =  map < int,Vecteur3D > ();
//    color_electronical_state[fond] = Green; //
//    color_electronical_state = ;
//    size_electronical_state;
//    scale_vibrational_state; // v
//    angle_rotational_state;  // J
//    angle_angular_state;     // M
};

**********/

void init_couleurs(Vecteur3D *couleur, const int Nmax)
{
    couleur[fond] = Green;
    couleur[excite] = Yellow;
    couleur[excite2] = Blue;
    couleur[excite3] = Olive;
    for (int i = excite3; i < Nmax; i++)
            couleur[i] = Vecteur3D	(50+2*i,50+2*i,0)	/	255;
}

// initialisation des tailles affichées sur l'écran
void init_tailles(double *taille, const int Nmax)
{
    taille[fond] = 0.01;
    taille[excite] = 0.015;
    taille[excite2] = 0.02;
    taille[excite3] = 0.03;
    for (int i = excite3; i < Nmax; i++)
        taille[i] =  0.005*i;
}



// dessin d'une particule unique
// Pour nous la taille est proportionnelle à la vibration
// L'angle de rotation est donné par la rotation J (selon Ox), MJ (selon 0y)
template <class Body>
void affichage_one_bod(const Body my_bod, FitParams &params)  // Affichage des points
{

    int etat = my_bod.exc;
    if (etat < 0) etat = excite3;     // choose the color for the lost (or ionized) particules to plot them
    if (my_bod.get_name() == "P_bar") etat = excite3;     // Then P_bar are red


    int v = my_bod.v;
    int two_J = my_bod.two_J;
    int two_M = my_bod.two_M;

    int Nb_max = 100; // Pour etre conservatif
    if (etat >= Nb_max) etat = excite2;     // To avoid to call array[etat] outside of its broder

    Vecteur3D couleur_bod[Nb_max];
    double tailles_bod[Nb_max];
    init_couleurs(couleur_bod,Nb_max);
    init_tailles(tailles_bod,Nb_max); // On met 100 et pas Nb_max car cela buggait ?????


    int two_J_maX= params.LocateParam("N_Two_JX_out_max")->val + params.LocateParam("N_Two_JA_out_max")->val; // J_max pour l'affichage
    int v_maX=params.LocateParam("NX_out")->val; // v_max pour l'affichage


    // glRotate produces a rotation of angle	degrees	around the vector (x,y,z) for all object
    glRotated(90.*two_J/(two_J_maX+1.1),0.,0.,1.);
    glRotated(90.*two_M/(two_J+1.1),0.,1.,0.); // +1.1 pour éviter le J=0 et convertir le int en double

    double DIAM = (1.0/50.0);
    double longueur = 0.2*(v+1.2)/(v_maX+1.1); // ou  1.1*v/(v_maX+1.1) si on veut l'état v=0 petit


    // glColor3d((MJ+J)/(2.*J)+couleur_bod[etat].x(),(MJ+2.*J)/(2.*J)*couleur_bod[etat].y(),couleur_bod[etat].z());
    // Pour l'état fondamental M=-J vert sombre  M=+J jaune clair

    glColor3d(couleur_bod[etat].x(),couleur_bod[etat].y(),couleur_bod[etat].z());

    glutSolidSphere(tailles_bod[etat],10,10);

    glPushMatrix();

    GLUquadricObj *quad;
    quad = gluNewQuadric();
// Je ne comprend pas pourquoi il faut cela
    glRotatef(90.,0.0,1.0,0.0);


    gluCylinder(quad, longueur*DIAM, longueur*DIAM, 0.8*longueur, 10, 10);
    glPopMatrix();

    gluDeleteQuadric(quad);

    glTranslatef(longueur,0.,0.); //position du 2ème atome
    glColor3d (1,1,1);
    glutSolidSphere(tailles_bod[etat],10,10);
}



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
    stat_molecule_un_niveau(bod, stat_Mol, Level,  fieldB, fieldE, all, Nb_body, params);

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

// Affiche la légende pour l'état des particules
template <class Body>
void affichage_stat_color(Body my_bod, const int Nb_state, FitParams &params)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW); //la pile des transformations "point de vue" est selectionnee
    glLoadIdentity(); // initialisation de la transformation. Evite que les rotations, translations, .. s’ajoutent à la transformation courante.


    glTranslatef(-0.95,0.95,0);

    for (int i=0; i!= Nb_state; i++)
    {
        glPushMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit
        my_bod.set_pos(Vecteur3D(0,-(i+2)/10.,0));

        affichage_one_bod(my_bod, params);
        glPopMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit

    }
}



// Affichage des particules
template <class Body>
void affichage_zone(const vector <Body> &bod, const int Nb_body, const double SIZE_real, FitParams &params)  // Affichage des points
{

    for (int i = 0; i < Nb_body; i++)
    {
        Vecteur3D position ( (bod[i].get_pos())/SIZE_real); // Ramène la position

        glPushMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit
        // glTranslatef(0,0,0);
        glTranslatef(position.x(),position.y(),-position.z());

        affichage_one_bod(bod[i], params);

        glPopMatrix(); // Dépile la matrice en haut de pile et remplace la matrice courante par celle-ci.
    }
    // glFlush(); // Vide la matrice de transformation de point de vue
}

// Affichage du système
template <class Body>
void Draw(const vector <Body> &bod, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const vector <Laser> &lasers,  const double SIZE_real, const double t, const int Nb_body, const int Nb_of_different_type_body, FitParams &params)
{

    glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT|GL_ACCUM_BUFFER_BIT);

    int Nb_state = Level.size();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum (-1.0, 1.0, -1.0, 1.0, 1, 3); // glFrustum(left, right, bas, haut, proche, loin);
    // Ici la caméra est à l'origine et l'écran va de X -1,1 en Y -1,1 (Attention Z<0 mais les coeef sont >0 !)


    glMatrixMode(GL_MODELVIEW); //la pile des transformations "point de vue" est selectionnee
    glLoadIdentity(); // initialisation de la transformation. Evite que les rotations, translations, .. s’ajoutent à la transformation courante.
    //glViewport(100,200,300,300); // glViewport(GLint x,GLint y,GLSizei largeur,GLSizei hauteur) 'x' et 'y' sont les coordonnées du coin supérieur left
    glTranslatef(0,0,-Z0);
    double angle = 90; // 90 pour avoir la gravité vers le bas (Oz vers le haut), 180 (0) pour avoir Oy vers le bas (haut)

    double rot_axe_x = params.LocateParam("rot_axe_x")->val;
    double rot_axe_y = params.LocateParam("rot_axe_y")->val;
    double rot_axe_z = params.LocateParam("rot_axe_z")->val;
    glRotated(angle,rot_axe_x,rot_axe_y,rot_axe_z);  // glRotated(angle ,rot_axe_x, rot_axe_y, rot_axe_z);

    // glScalef(1.,1.,1.);            // multiplie la taille de l'image par zoom (effectue un zoom)
// On pourrait faire le zoom avec la souris et le fichier zpr.c

    affichage_zone(bod, Nb_body, SIZE_real, params);  // Affiche les particules

    traceRepere (0.5); // // Cette fonction trace un repère de taille longueur,
// centré sur le repère objet avec
// une flêche rouge le long de l'axe des x,
// une flêche verte le long de l'axe des y
// et une flêche bleue le long de l'axe des z.

    affichage_legend(bod, Level,  fieldB, fieldE, t, Nb_body,Nb_state, params) ;     // Affiche les statistiques sur les molecules

    // affichage_stat_color(bod[0],Nb_state, params);  // Affiche les statistiques sur les molecules

    glutSwapBuffers(); // A TESTER permet d'eviter un "blincking" lors de l'affichage

    // wait(1.5);   // Permet de ne pas avoir un affichage trop rapide


}





#endif




