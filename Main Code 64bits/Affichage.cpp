
/*
  Name:
  Copyright:
  Author:
  Date: 23/10/06 10:01
  Description: Affichage de particules


Les N particules sont des Vecteur3D et sont regroupées dans une
  classe Body quelconque (pour plus de fléxibilité) (bod[i] est la particule i).
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

L'unité de base est donc SIZE_real = SIZE_affichage;
*/



#include "Affichage.h"


/************************************************************************/
/************************** Programme *************************/
/************************************************************************/


// Cette fonction trace un repère de taille longueur,
// centré sur le repère objet avec
// une flêche rouge le long de l'axe des x,
// une flêche verte le long de l'axe des y
// et une flêche bleue le long de l'axe des z.

// http://www.irit.fr/~Loic.Barthe/Enseignements/TPs_OpenGL/M1_IO5/TP8/tp8.html
void traceRepere (float longueur)
{

  float DIAMETRECYLINDRE = 1.0/30.0;
  float DIAMETRECONE = 1.0/15.0;

  GLUquadricObj *quad;
  quad = gluNewQuadric();


  /* Axe X */
  glColor3d(1.0,0.0,0.0);
  glPushMatrix();
  glRotatef(90.,0.0,1.0,0.0);
  gluCylinder(quad, longueur*DIAMETRECYLINDRE, longueur*DIAMETRECYLINDRE, 0.8*longueur, 10, 10);
  glTranslatef(0.0,0.0,0.8*longueur);
  gluCylinder(quad, longueur*DIAMETRECONE, 0.0, 0.2*longueur, 10, 10);
  glPopMatrix();
  /* Axe Y */
  glColor3d(0.0,1.0,0.0);
  glPushMatrix();
  glRotatef(-90.,1.0,0.0,0.0);
  gluCylinder(quad, longueur*DIAMETRECYLINDRE, longueur*DIAMETRECYLINDRE, 0.8*longueur, 10, 10);
  glTranslatef(0.0,0.0,0.8*longueur);
  gluCylinder(quad, longueur*DIAMETRECONE, 0.0, 0.2*longueur, 10, 10);
  glPopMatrix();
  /* Axe Z */
  glColor3d(0.0,0.0,1.0);
  glPushMatrix();
  glRotatef(180.,1.0,0.0,0.0);
  gluCylinder(quad, longueur*DIAMETRECYLINDRE, longueur*DIAMETRECYLINDRE, 0.8*longueur, 10, 10);
  glTranslatef(0.0,0.0,0.8*longueur);
  gluCylinder(quad, longueur*DIAMETRECONE, 0.0, 0.2*longueur, 10, 10);
  glPopMatrix();


  gluDeleteQuadric(quad);

}


// Fonction pour recadrer l'image lors de sa première apparition ou si on la bouge avec la souris
void reshape(int w,int h)
{
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,w,0,h);
  glScalef(1,-1,1);
  glTranslatef(0,-h,0);
  glMatrixMode(GL_MODELVIEW);
}


//  Permet d'afficher du texte avec GLUT en x,y
// essayer  ModuleFont.h
//  http://raphaello.univ-fcomte.fr/IG/Modules/Modules.htm
void bitmap_output(int x,int y,char *string,void *font)
{
  int len,i;
  glRasterPos2f(x,y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)
    {
      glutBitmapCharacter(font,string[i]);
    }
}

// Procedure simplifiant l'écriture
void texte_output(const double nombre, string text_string,const int x, const int y)
{
  char nb[16];
  sprintf(nb,"%d",(int) nombre);
  string snb = nb ;
  text_string +=  snb;
  char texte[256]; // Doit être assez long
  strcpy(texte, text_string.c_str());
  bitmap_output(x,y,texte,GLUT_BITMAP_HELVETICA_18);
}

// Affiche une ellipse de taille Radius3D
// couleur par défault blanc
// nb de polygônes 20 par défault pour simuler la sphère
void affichage_ellipse(const Vecteur3D Radius3D, const Vecteur3D color, const int nb_polygone, const double SIZE_real)
{
  glPushMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit
  glScalef(Radius3D.x()/SIZE_real,Radius3D.y()/SIZE_real,Radius3D.z()/SIZE_real);
  // Donne la bonne forme elliptique

  GLUquadricObj *q = gluNewQuadric();
  gluQuadricDrawStyle(q, GLU_LINE);
  glColor3d(color.x(), color.y(), color.z());    // Blanc
  double rayon=Radius3D.mag()/SIZE_real;  // Dessine le rayon
  gluSphere(q, rayon, nb_polygone, nb_polygone);     // divisée en nb_polygone polygones des 2 axes
  gluDeleteQuadric(q);
  glPopMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit
}

/************
// Affiche un hyperboloide de révolution (pas encore vraiment opérationnel)
// D'après http://www-math.edu.kagoshima-u.ac.jp/~fractaro/AnschaulicheGeometrie/Hyperboloid.c

// int mode = 0;
// int sign = 1;
// GLfloat blue[4], cyan[4], orange[4], yellow[4];
**************/



double length(double v[])
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void exterior(double v1[], double v2[], double v3[])
{
  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void cylinder(double a[], double b[], double r, int slices, GLfloat color[])
{
  double e1[3], e2[3], e3[3], va[40][3], vb[40][3], n[40][3], e4[3];
  double h, sinbeta, cosalpha, sinalpha, t, c, s;
  int m, i, m1;
  double epsilon = 0.00000001;

  for (i = 0; i < 3; i++)
    {
      e3[i] = b[i] - a[i];
    }
  h = length(e3);
  for (i = 0; i < 3; i++)
    {
      e3[i] /= h;
    }

  sinbeta = sqrt(e3[0] * e3[0] + e3[1] * e3[1]);
  if (sinbeta > epsilon)
    {
      cosalpha = e3[0] / sinbeta;
      sinalpha = e3[1] / sinbeta;
      e1[0] = -sinalpha;
      e1[1] = cosalpha;
      e1[2] = 0;
    }
  else
    {
      e1[0] = 1;
      e1[1] = e1[2] = 0;
    }
  exterior(e3, e1, e2);

  for (m = 0; m < slices; m++)
    {
      t = 2.0 * pi / slices * m;
      c = cos(t);
      s = sin(t);
      for (i = 0; i < 3; i++)
        {
          va[m][i] = a[i] + r * (c * e1[i] + s * e2[i]);
          vb[m][i] = va[m][i] + h * e3[i];
        }
    }

  for (m = 0; m < slices; m++)
    {
      t = 2.0 * pi / slices * (m + 0.5);
      c = cos(t);
      s = sin(t);
      for (i = 0; i < 3; i++)
        {
          n[m][i] = c * e1[i] + s * e2[i];
        }
    }

  for (i = 0; i < 3; i++) e4[i] = -e3[i];

  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color);

  for (m = 0; m < slices; m++)
    {
      m1 = m + 1;
      if (m1 == slices) m1 = 0;
      glBegin(GL_POLYGON);
      glNormal3dv(&n[m][0]);
      glVertex3dv(&vb[m][0]);
      glVertex3dv(&va[m][0]);
      glVertex3dv(&va[m1][0]);
      glVertex3dv(&vb[m1][0]);
      glEnd();
    }

  for (m = 0; m < slices; m++)
    {
      m1 = m + 1;
      if (m1 == slices) m1 = 0;
      glBegin(GL_POLYGON);
      glNormal3dv(&e4[0]);
      glVertex3dv(&a[0]);
      glVertex3dv(&va[m1][0]);
      glVertex3dv(&va[m][0]);
      glEnd();
    }

  for (m = 0; m < slices; m++)
    {
      m1 = m + 1;
      if (m1 == slices) m1 = 0;
      glBegin(GL_POLYGON);
      glNormal3dv(&e3[0]);
      glVertex3dv(&b[0]);
      glVertex3dv(&vb[m][0]);
      glVertex3dv(&vb[m1][0]);
      glEnd();
    }
}


// Affiche un hyperboloide de révolution
// D'après http://www-math.edu.kagoshima-u.ac.jp/~fractaro/AnschaulicheGeometrie/Hyperboloid.c
// hyperboloid of one sheet is expressed by an equation
// x^2/a^2 + y^2/b^2 − z^2/c^2 = 1. I.E. taille a selon x, b selon y et c selon z
// n = nb de ligne le formant
// h = ?
// r = rayon des cylindre le formant
void Hyperboloid(double a, double b, double c, int n, double h, double r)
{


  double phi, dphi, cosphi, sinphi, w, t, P[3], e1[3], e2[3], A1[3], B1[3], A2[3], B2[3];
  int i, m;

  dphi = 2.0 * pi / n;

  for (m = 0; m < n; m++)
    {
      phi = dphi * m;
      cosphi = cos(phi);
      sinphi = sin(phi);
      w = 1.0 / sqrt(c * c + b * b * cosphi * cosphi + a * a * sinphi * sinphi);
      t = h / (c * w);
      P[0] = a * cosphi;
      P[1] = b * sinphi;
      P[2] = 0;
      e1[0] = -a * sinphi * w;
      e1[1] = b * cosphi * w;
      e1[2] = c * w;
      e2[0] = -e1[0];
      e2[1] = -e1[1];
      e2[2] = e1[2];
      for (i = 0; i < 3; i++)
        {
          A1[i] = P[i] - t * e1[i];
          B1[i] = P[i] + t * e1[i];
          A2[i] = P[i] - t * e2[i];
          B2[i] = P[i] + t * e2[i];
        }

GLfloat blue[4], cyan[4], orange[4], yellow[4];

      cylinder(A1, B1, r, 10, yellow);
      cylinder(A2, B2, r, 10, blue);
    }
}

// Affichage des lasers
void affichage_laser(const Laser my_laser, const double SIZE_real, int /* n */, double /* h */, double /* r */)
{

//    glPushMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit
//    glLoadIdentity(); // initialisation de la transformation. Evite que les rotations, translations, .. s’ajoutent à la transformation courante.
//    // glScalef(1./SIZE_real,1./SIZE_real,1./SIZE_real); // Donne la bonne taille
//
//
//    Vecteur3D waist_position= my_laser.get_waist_pos();
//    Vecteur3D waist=my_laser.waist_size(waist_position); // wx,wy,zone de Rayleigh moyenne
//    Vecteur3D direction= my_laser.get_direction();
//
//    Vecteur3D rot=direction.cross(Vecteur3D(0.,0.,1.)); // rotation Euler  pour amener OZ sur k
//    double angle_rot= direction.angle(Vecteur3D(0.,0.,1.));
//
//    glTranslatef(waist_position.x(),waist_position.y(),waist_position.z()); // centre le laser sur son waist
//    glRotated(angle_rot ,rot.x(), rot.y(), rot.z()); // glRotated(angle ,rot_axe_x, rot_axe_y, rot_axe_z);
//
//    Hyperboloid(waist.x(), waist.y(), waist.z(), n,  h, r);
//
//    glPopMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit


  Vecteur3D waist_position= my_laser.get_waist_pos();
  Vecteur3D waist=my_laser.waist_size(waist_position); // wx,wy,zone de Rayleigh moyenne

  glPushMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit
  glLoadIdentity(); // initialisation de la transformation. Evite que les rotations, translations, .. s’ajoutent à la transformation courante.

  glScalef(10,10,10); // Donne la bonne taille
  glRotated(45,1,1,1); // glRotated(angle ,rot_axe_x, rot_axe_y, rot_axe_z);
  Hyperboloid(waist.x(), waist.y(), waist.z(), 10, 10, 0.01);
  glPopMatrix(); // Empile la matrice courante dans la pile de matrices. Permet d'effectuer la transformation uniquement sur ce qui suit
}






