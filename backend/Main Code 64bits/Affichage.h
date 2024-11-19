<<<<<<< HEAD
/*
  Name:
  Author:
  Date: 23/10/06 10:01
  Description: Displaying particles using OpenGL.

  The N particles are stored as Vecteur3D objects within a general Body class
  (to allow flexibility). For instance, bod[i] represents the i-th particle.
  The Body class must implement the following functions:

      - bod.get_pos: to retrieve the position.
      - couleur_bod[N]: provides the color.

This uses OPEN_GL.
Here is my understanding of OpenGL, see additional references in readMe_graphique.txt:

- glMatrixMode(GL_MODELVIEW): Selects the transformation stack for the "point of view".
- glLoadIdentity(): Initializes the transformation to avoid additive effects
  from previous rotations, translations, etc.
- glRotate, glTranslatef, etc.: Creates 4x4 matrices that apply transformations
  and stack them within the "point of view" stack.
- Drawing functions (Vertex, sphere, etc.): Draw shapes and apply transformations.

If you want to exclude a transformation from subsequent drawings:
- Use `glPushMatrix()` before the excluded drawing and `glPopMatrix()` after.

For perspective projection:
- glMatrixMode(GL_PROJECTION);
- glLoadIdentity();
- glFrustum(left, right, bottom, top, near, far);

The base unit is defined as SIZE_real = SIZE_display.
*/

#ifndef AFFICHAGE_SEEN
#define AFFICHAGE_SEEN

#include <iostream>
#include <ctime>                // For clock()
#include <windows.h>

#ifdef __APPLE__
#include <GL/glut.h>
#else
#include <GL/glut.h>
#endif

#include <GL/gl.h>              // OPENGL
#include <GL/glu.h>             // OPENGL (GLU)
#include <GL/glut.h>            // GLUT library
#include <GL/glext.h>

#include "Vecteur3D.h"
#include "Laser.h"              // Laser class

using namespace std;

// Constants for the screen and z-axis offset in OpenGL
const double SIZE_SCREEN = 2.0; // The screen ranges in X and Y between -1 and 1, and in Z between -1 and -3.
const double Z0 = 2.0;          // OpenGL requires Z > 0, so a Z offset is applied.

/************************************************************************/
/************************** Program Functions ***************************/
/************************************************************************/

/**
 * Draws a coordinate system centered on the object’s reference frame.
 * Axes are represented as:
 * - Red arrow along the x-axis.
 * - Green arrow along the y-axis.
 * - Blue arrow along the z-axis.
 *
 * @param length Length of each axis.
 */
void traceRepere(float length);

/**
 * Rescales the display when it first appears or when resized by the user.
 *
 * @param w New width of the window.
 * @param h New height of the window.
 */
void reshape(int w, int h);

/**
 * Displays text using GLUT at the given x, y coordinates.
 *
 * @param x X-coordinate for the text.
 * @param y Y-coordinate for the text.
 * @param string Text to display.
 * @param font Font to use for rendering.
 */
void bitmap_output(int x, int y, char *string, void *font);

/**
 * Simplifies displaying text output by combining numeric values and a string.
 *
 * @param number Number to display.
 * @param text_string Text to prepend to the number.
 * @param x X-coordinate for the text.
 * @param y Y-coordinate for the text.
 */
void texte_output(const double number, string text_string, const int x, const int y);

/**
 * Displays an ellipse with a given radius, color, and number of polygons.
 *
 * @param Radius3D Radius of the ellipse in 3D space.
 * @param color Color of the ellipse (default is white).
 * @param nb_polygone Number of polygons to approximate the ellipse (default: 20).
 * @param SIZE_real Scaling factor for the size.
 */
void affichage_ellipse(const Vecteur3D Radius3D, const Vecteur3D color = Vecteur3D(1.0, 1.0, 1.0), const int nb_polygone = 20, const double SIZE_real = 0.005);

/**
 * TODO: Improve handling of the hyperboloid geometry and its parameters.
 * Displays a hyperboloid of revolution.
 * Equation: x^2/a^2 + y^2/b^2 − z^2/c^2 = 1.
 *
 * @param a Radius along x-axis.
 * @param b Radius along y-axis.
 * @param c Radius along z-axis.
 * @param n Number of lines forming the hyperboloid.
 * @param h Height of the hyperboloid.
 * @param r Radius of the cylinders forming the hyperboloid.
 */
void Hyperboloid(double a, double b, double c, int n, double h, double r);

/**
 * Calculates the length (Euclidean norm) of a 3D vector.
 *
 * @param v 3D vector as an array.
 * @return Length of the vector.
 */
double length(double v[]);

/**
 * Computes the cross product of two 3D vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @param v3 Resulting vector (cross product of v1 and v2).
 */
void exterior(double v1[], double v2[], double v3[]);

/**
 * Draws a cylinder between two points with a given radius and color.
 *
 * @param a Starting point.
 * @param b Ending point.
 * @param r Radius of the cylinder.
 * @param slices Number of divisions for the cylinder's surface.
 * @param color Color of the cylinder.
 */
void cylinder(double a[], double b[], double r, int slices, GLfloat color[]);

/**
 * TODO: Implement better representation for laser visualization.
 * Displays laser beams.
 *
 * @param my_laser Laser object to display.
 * @param SIZE_real Scaling factor for size.
 * @param n Number of lines forming the display (optional).
 * @param h Height of the display (optional).
 * @param r Radius of the display (optional).
 */
void affichage_laser(const Laser my_laser, const double SIZE_real, int n, double h, double r);

#endif


=======
/*
  Name:
  Author:
  Date: 23/10/06 10:01
  Description: Displaying particles using OpenGL.

  The N particles are stored as Vecteur3D objects within a general Body class
  (to allow flexibility). For instance, bod[i] represents the i-th particle.
  The Body class must implement the following functions:

      - bod.get_pos: to retrieve the position.
      - couleur_bod[N]: provides the color.

This uses OPEN_GL.
Here is my understanding of OpenGL, see additional references in readMe_graphique.txt:

- glMatrixMode(GL_MODELVIEW): Selects the transformation stack for the "point of view".
- glLoadIdentity(): Initializes the transformation to avoid additive effects
  from previous rotations, translations, etc.
- glRotate, glTranslatef, etc.: Creates 4x4 matrices that apply transformations
  and stack them within the "point of view" stack.
- Drawing functions (Vertex, sphere, etc.): Draw shapes and apply transformations.

If you want to exclude a transformation from subsequent drawings:
- Use `glPushMatrix()` before the excluded drawing and `glPopMatrix()` after.

For perspective projection:
- glMatrixMode(GL_PROJECTION);
- glLoadIdentity();
- glFrustum(left, right, bottom, top, near, far);

The base unit is defined as SIZE_real = SIZE_display.
*/

#ifndef AFFICHAGE_SEEN
#define AFFICHAGE_SEEN

#include <iostream>
#include <ctime>                // For clock()
#include <windows.h>

#ifdef __APPLE__
#include <GL/glut.h>
#else
#include <GL/glut.h>
#endif

#include <GL/gl.h>              // OPENGL
#include <GL/glu.h>             // OPENGL (GLU)
#include <GL/glut.h>            // GLUT library
#include <GL/glext.h>

#include "Vecteur3D.h"
#include "Laser.h"              // Laser class

using namespace std;

// Constants for the screen and z-axis offset in OpenGL
const double SIZE_SCREEN = 2.0; // The screen ranges in X and Y between -1 and 1, and in Z between -1 and -3.
const double Z0 = 2.0;          // OpenGL requires Z > 0, so a Z offset is applied.

/************************************************************************/
/************************** Program Functions ***************************/
/************************************************************************/

/**
 * Draws a coordinate system centered on the object’s reference frame.
 * Axes are represented as:
 * - Red arrow along the x-axis.
 * - Green arrow along the y-axis.
 * - Blue arrow along the z-axis.
 *
 * @param length Length of each axis.
 */
void traceRepere(float length);

/**
 * Rescales the display when it first appears or when resized by the user.
 *
 * @param w New width of the window.
 * @param h New height of the window.
 */
void reshape(int w, int h);

/**
 * Displays text using GLUT at the given x, y coordinates.
 *
 * @param x X-coordinate for the text.
 * @param y Y-coordinate for the text.
 * @param string Text to display.
 * @param font Font to use for rendering.
 */
void bitmap_output(int x, int y, char *string, void *font);

/**
 * Simplifies displaying text output by combining numeric values and a string.
 *
 * @param number Number to display.
 * @param text_string Text to prepend to the number.
 * @param x X-coordinate for the text.
 * @param y Y-coordinate for the text.
 */
void texte_output(const double number, string text_string, const int x, const int y);

/**
 * Displays an ellipse with a given radius, color, and number of polygons.
 *
 * @param Radius3D Radius of the ellipse in 3D space.
 * @param color Color of the ellipse (default is white).
 * @param nb_polygone Number of polygons to approximate the ellipse (default: 20).
 * @param SIZE_real Scaling factor for the size.
 */
void affichage_ellipse(const Vecteur3D Radius3D, const Vecteur3D color = Vecteur3D(1.0, 1.0, 1.0), const int nb_polygone = 20, const double SIZE_real = 0.005);

/**
 * TODO: Improve handling of the hyperboloid geometry and its parameters.
 * Displays a hyperboloid of revolution.
 * Equation: x^2/a^2 + y^2/b^2 − z^2/c^2 = 1.
 *
 * @param a Radius along x-axis.
 * @param b Radius along y-axis.
 * @param c Radius along z-axis.
 * @param n Number of lines forming the hyperboloid.
 * @param h Height of the hyperboloid.
 * @param r Radius of the cylinders forming the hyperboloid.
 */
void Hyperboloid(double a, double b, double c, int n, double h, double r);

/**
 * Calculates the length (Euclidean norm) of a 3D vector.
 *
 * @param v 3D vector as an array.
 * @return Length of the vector.
 */
double length(double v[]);

/**
 * Computes the cross product of two 3D vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @param v3 Resulting vector (cross product of v1 and v2).
 */
void exterior(double v1[], double v2[], double v3[]);

/**
 * Draws a cylinder between two points with a given radius and color.
 *
 * @param a Starting point.
 * @param b Ending point.
 * @param r Radius of the cylinder.
 * @param slices Number of divisions for the cylinder's surface.
 * @param color Color of the cylinder.
 */
void cylinder(double a[], double b[], double r, int slices, GLfloat color[]);

/**
 * TODO: Implement better representation for laser visualization.
 * Displays laser beams.
 *
 * @param my_laser Laser object to display.
 * @param SIZE_real Scaling factor for size.
 * @param n Number of lines forming the display (optional).
 * @param h Height of the display (optional).
 * @param r Radius of the display (optional).
 */
void affichage_laser(const Laser my_laser, const double SIZE_real, int n, double h, double r);

#endif


>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
