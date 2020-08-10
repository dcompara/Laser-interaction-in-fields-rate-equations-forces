/*
  Name:
  Copyright:
  Author:
  Date: 23/10/06 10:01
  Description: Affichage de particules


Les N particules sont des Vecteur3D et sont regroupées dans une
  classe Body quelconque (pour plus de flexibilité) (bod[i] est la particule i).
      La classe Body qui doit contenir les fonctions:

         bod.get_pos  pour avoir la position
         couleur_bod[N] donne la couleur


Utilise OPEN_GL
Voici ce que j'ai compris d'OPENGL voir les autres références dans readMe_graphique.txt

glMatrixMode(GL_MODELVIEW); //la pile des transformations "point de vue" est selectionnee
glLoadIdentity(); // initialisation de la transformation. Evite que les rotations, translations, .. s’ajoutent à la transformation courante.
glRotate, glTranslatef, ... // Crée chaque fois une matrice 4*4 qui tient compte de chaque transformation et qui s'empile dans la pile de "point de vue"
fonction de dessin1 Vertex, sphère etc...
fonction de dessin2 Vertex, sphère etc...
 // Effectue les dessins, puis effectue l'ensemble des matrices puis la projection

Si on veut dépiller par exemple que dessin 2 soit sans la dernière translation. Il faut faire
glPushMatrix();
 fonction de dessin2 Vertex, sphère etc...
glPopMatrix();
fonction de dessin1 Vertex, sphère etc...

Pour la spécification d’une projection perspective :
glMatrixMode(GL_PROJECTION);
glLoadIdentity();
glFrustum(left, right, bas, haut, proche, loin);

L'unité de base est donc SIZE_real = SIZE_affichage
*/






#ifndef Affichage_SEEN
#define Affichage_SEEN



#include <iostream>


#include <time.h>                // For clock()

#include <windows.h>
#ifdef __APPLE__
#include <GL/glut.h>
#else
#include <GL/glut.h>
#endif

#include <GL/gl.h>                          // OPENGL
#include <GL/glu.h>                         // OPENGL (GLU)
#include <gl/glut.h>                       // glut bibliotheque
// #include <gl/freeglut.h>                       // glut bibliotheque
#include <GL/glext.h>

#include "Vecteur3D.h"
#include "Laser.h"                       // Classe Laser

using namespace std;

const double SIZE_SCREEN =  2.; //  L'écran est en X -1,1 en Y -1,1 et en Z -1,-3 donc SIZE_SCREEN = 2
const double Z0 = 2.; // dans OPEN_GL il faut Z >0 donc je translate de Z0;


/************************************************************************/
/************************** Programme *************************/
/************************************************************************/


// Cette fonction trace un repère de taille longueur,
// centré sur le repère objet avec
// une flêche rouge le long de l'axe des x,
// une flêche verte le long de l'axe des y
// et une flêche bleue le long de l'axe des z.

// http://www.irit.fr/~Loic.Barthe/Enseignements/TPs_OpenGL/M1_IO5/TP8/tp8.html
void traceRepere (float longueur);

// Fonction pour recadrer l'image lors de sa première apparition ou si on la bouge avec la souris
void reshape(int w,int h);

//  Permet d'afficher du texte avec glut en x,y
// essayer  ModuleFont.h
//  http://raphaello.univ-fcomte.fr/IG/Modules/Modules.htm
void bitmap_output(int x,int y,char *string,void *font);

// Procedure simplifiant l'écriture
void texte_output(const double nombre, string text_string,const int x, const int y);

// Affiche une ellipse de taille Radius3D
// couleur par défault blanc
// nb de polygônes 20 par défault pour simuler la sphère
void affichage_ellipse(const Vecteur3D Radius3D, const Vecteur3D color=Vecteur3D(1.,0.,1.), const int nb_polygone=20, const double SIZE_real = 0.005);

/************
// Pour pouvoir afficher des faisceaux lasers comme des hyperboloides de révolution
// D'après http://www-math.edu.kagoshima-u.ac.jp/~fractaro/AnschaulicheGeometrie/Hyperboloid.c
**************/

//
//int mode = 0;
//int sign = 1;
//
//GLfloat blue[4], cyan[4], orange[4], yellow[4];

double length(double v[]);
void exterior(double v1[], double v2[], double v3[]);
void cylinder(double a[], double b[], double r, int slices, GLfloat color[]);

// Affiche un hyperboloide de révolution
// D'après http://www-math.edu.kagoshima-u.ac.jp/~fractaro/AnschaulicheGeometrie/Hyperboloid.c
// hyperboloid of one sheet is expressed by an equation
// x^2/a^2 + y^2/b^2 − z^2/c^2 = 1. I.E. taille a selon x, b selon y et c selon z
// n = nb de ligne le formant
// h = ?
// r = rayon des cylindre le formant
void Hyperboloid(double a, double b, double c, int n, double h, double r);

// Affichage des lasers
void affichage_laser(const Laser my_laser, const double SIZE_real, int n, double h, double r);


#endif




