/*
Name:
Copyright:
Author:
Date: 12/10/06 21:38
Description: Algorithm for the N or 1 body dynamics


see: N-body simulations of gravitational dynamics (http://arxiv.org/pdf/1105.1082.pdf)

For the N-body I use nbody_sh1WD.C:  an N-body integrator with a shared but variable time step
by Walter Dehnen

Elle est effectue les calcul sur les positions, vitesses, accélération de n particules ayant une masse
qui sont des Vecteur3D (c.f. classe Vecteur3D.h)

     bod.get_vel  pour avoir la  vitesse, il faut aussi la pos, acc, jerk
     bod.inc_vel(d_vel) et pos pour augmenter la vitesse (et la position)
     bod.set_vel(d_vel) et pos, acc, jerk pour mettre la vitesse (et la position, acc, jerk)

On utilise une fonction new_acc_pot(bod) qui calcule l'accélération via new_pot(bod) qui recalcule (et non lit ça c'est get_pot!) le potentiel
Mais si l'accélération est mise à jour ailleurs ou peut être calculé directement on utilise plutôt new_acc(bod).

Pour l'instant si on veut la force dipolaire on utilise les potentiels (Verlet_pot) et sinon l'accélération ((Verlet_acc)

Comme l'accélération dépent des forces, elle dépends des potentiels magnétique, électrique et laser qui peuvent démendre du temps

ATTENTION: Elle ne met pas à jour les potentiels une fois le mouvement fait. Il faut le faire à part.

Si la particule est trop loin (LARGE_NUMBER) on ne la fait pas évoluer
*/

#ifndef one_body_SEEN
#define one_body_SEEN

using namespace std;

#include  <iostream>
#include  <cmath>                          // to include sqrt(), etc.
#include  <cstdlib>                        // for atoi() and atof()
#include  <unistd.h>                       // for getopt()
#include  <fstream>

#include  "molecule.h"              // Pour le shift en énergie delta
#include  "Field.h"                 // Pour les champs extérieurs
#include  "laser.h"
#include "params.h"

enum N_Body_algorithmes  // Hermite = 0, // NOT USED
{
    Aucun_N_corps = -1,
    Hermite = 0, // Not used
    Verlet_acc = 1, // Calcule l'accélération directement
    Verlet_pot = 2, // Calcule l'accélération avec les potentiels
    Yoshida6_acc = 3,       // 6 th order symplectic algorithm based on Verlet_acc
    Yoshida6_pot = 4,        // 6 th order symplectic algorithm based on Verlet_pot
    Runge_Kutta_Nystrom = 5, // NOT USED
    Verlet_pot_gradient_high_order = 6, // Calcule l'accélération avec les potentiels mais à l'ordre supérieur
    Boris_Buneman = 7 // Boris_Buneman integrator
};



template <class Body> void evolve_step(const N_Body_algorithmes Algorithme_N_body,
                                       vector <Body> &bod,  const Field &fieldB, const Field &fieldE, const  vector <Laser> &my_laser,
                                       const int n, double & t,const double dt, FitParams &params);

// Juste l'accélération
template <class Body>  void  evolve_step_Aucun_N_corps(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt);


template <class Body> void evolve_step_Verlet_via_acc(vector <Body> &bod,  const Field &fieldB, const Field &fieldE, const  vector <Laser> &my_laser,
        const int n, double & t,const double dt);

// Accélération calculée avec new_acc_pot qui utilise les potentiel des champs.
// C'est plus lent car l'accélération est évaluée par pot(x+dx)-pot(x-dx)/2dx ...
template <class Body>  void  evolve_step_Verlet_via_pot(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt, FitParams &params);


// Yoshida 6 using LeapForg with acceleration
template <class Body>
void  evolve_step_Yoshida6_acc(vector <Body> &bod,  const Field &fieldB, const Field &fieldE, const  vector <Laser> &my_laser,
                               const int n, double & t,const double dt);

// Yoshida 6 using LeapForg with potentials
template <class Body>
void  evolve_step_Yoshida6_pot(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
                               const int n, double & t,const double dt, FitParams &params);


// Boris–Buneman integration scheme (notation of Journal of Computational Physics 273 (2014) 255)
template <class Body>  void  evolve_step_Boris(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt);

// Energie potentielle de la particle i
// Définie en extern (pour nous dans la classe Molecule) cf double  get_pot(const Molecule &mol, const Field &fieldB, const Field &fieldE);
template <class Body> double get_pot(vector <Body> &bod, const int i, const Field &fieldB, const Field &fieldE);


// Pour calcul a_i= Force/masse = -Grad(E_pot)/mass et met à jour l'accélération
template <class Body>  Vecteur3D new_acc_pot(vector <Body> &bod, const int i,  const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params);
// Calcul du gradient au premier ordre
template <class Body>  Vecteur3D new_acc_pot_second_order(vector <Body> &bod, const int i, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params);
// Calcul du gradient à l'ordre d'après
template <class Body>  Vecteur3D new_acc_pot_fourth_order(vector <Body> &bod, const int i, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params);

// Calculate the Coulombian acceleration (adapted from Walter Dehnen)
// Is added to accell
template <class Body>  void  add_acc_N_body_Coulomb(const vector <Body> &bod, Vecteur3D *accell, const int n);


// Calculate the Coulombian potential for the sum_j=1^n not i  qj/(r_ij)/4 pi epsilon_0
// Is used only for statistics to check the energy conservation. So it is not fully optimized
template <class Body>  double  get_coulomb_potential(const vector <Body> &bod, const Vecteur3D &ri, const int n);

// add Lorentz force for charged particles
template <class Body>  void  add_acc_Charged_particles(vector <Body> &bod, Vecteur3D *accell, const Field &fieldB, const Field &fieldE, const int n);

// kinetic energy
template <class Body> double get_kin(Body &bod);

// temps d'évolution typique des forces (ou potentiels) ~ fraction de waist/vitesse
// A améliorer avec les champs, différents lasers ...
template <class Body> double t_evol_ext(vector <Body> &bod,  const Field &fieldB, const Field &fieldE, const  vector <Laser> &my_laser, const double Temp_estime);


/**************************************
JE MET ICI L'equivalent du .cpp car "On peut pas séparer un .h et .cpp pour les template".

Si un symbole template est utilisé uniquement dans un .cpp (fichier source), il peut être implémenté dans ce .cpp
Sinon, il doit être implémenté dans un .hpp (header).

cf http://www.commentcamarche.net/faq/11194-les-templates-en-c
***************************************/



template <class Body> void evolve_step(const N_Body_algorithmes Algorithme_N_body,
                                       vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
                                       const int n, double & t,const double dt, FitParams &params)
{
    Body *old_bod = new Body[n];

    switch (Algorithme_N_body)
    {

    case Aucun_N_corps :
        evolve_step_Aucun_N_corps(bod, fieldB, fieldE, laser, n, t, dt);
        break;


    case Verlet_acc : // Accélaration calculée directement
        evolve_step_Verlet_via_acc(bod, fieldB, fieldE, laser, n, t, dt);
        break;

    case Verlet_pot : // Accélaration calculée avec gradient du potentiel
        evolve_step_Verlet_via_pot(bod, fieldB, fieldE, laser, n, t, dt, params);
        break;

    case Yoshida6_acc : // Accélaration calculée avec gradient du potentiel
        evolve_step_Yoshida6_acc(bod, fieldB, fieldE, laser, n, t, dt);
        break;

    case Yoshida6_pot : // Accélaration calculée avec gradient du potentiel
        evolve_step_Verlet_via_pot(bod, fieldB, fieldE, laser, n, t, dt, params);
        break;


    case Verlet_pot_gradient_high_order : // Accélaration calculée avec gradient du potentiel
        evolve_step_Verlet_via_pot(bod, fieldB, fieldE, laser, n, t, dt,  params);
        break;

    case Boris_Buneman : // Boris–Buneman integration scheme
        evolve_step_Boris(bod, fieldB, fieldE, laser, n, t, dt);
        break;

    default:
        /* do nothing */
        break;

    }
    t += dt;
    delete [] old_bod;
}



/*-----------------------------------------------------------------------------
 *  evolve_step  Velocity Verlet

old_acc = acc
@pos += @vel*dt + old_acc*0.5*dt*dt
new_acc = acc
@vel += (old_acc + new_acc)*0.5*dt


 *-----------------------------------------------------------------------------
*/


// Accélération calculée avec new_acc qui utilise directement les formules des gradients des champs. C'est plus rapide
// N'est pas utilisable toujours si on a pas les gradients des champs.
template <class Body>  void  evolve_step_Verlet_via_acc(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt)
{

    Vecteur3D *old_acc = new Vecteur3D[n];
    Vecteur3D *acc = new Vecteur3D[n];

// Might be improved if the acc is known before. because here we reclculate it!
    for (int i = 0; i < n ; ++i)
        old_acc[i] = new_acc(bod, i, fieldB, fieldE, laser, t); // calcul l'accélération

    add_acc_N_body_Coulomb(bod, old_acc, n); // add coulombian acceleration due to charged particles
    add_acc_Charged_particles(bod, old_acc, fieldB, fieldE, n); // add Lorentz force for charged particles

    for (int i = 0; i != n ; ++i)
    {
        if (bod[i].get_pos().mag() > LARGE_NUMBER)
            continue; // Si la particule est trop loin on la  laiss là!
        Vecteur3D d_pos = (dt*(bod[i].get_vel()) + (dt*dt*0.5)*old_acc[i]);
        bod[i].inc_pos(d_pos);   // @pos += @vel*dt + old_acc*0.5*dt*dt
    }

    for (int i = 0; i != n ; ++i)
        acc[i] = new_acc(bod, i, fieldB, fieldE, laser, t + dt);        // new_acc = acc

    add_acc_N_body_Coulomb(bod, acc, n); // add coulombian acceleration due to charged particles
    add_acc_Charged_particles(bod, acc, fieldB, fieldE, n); // add Lorentz force for charged particles

    for (int i = 0; i != n ; ++i)
    {
        Vecteur3D d_vel = (0.5*dt*(old_acc[i] + acc[i]));
        bod[i].inc_vel(d_vel);   // @vel += (old_acc + new_acc)*0.5*dt
    }

    delete [] acc;
    delete [] old_acc;

}




// Accélération calculée avec new_acc_pot qui utilise les potentieldes champs.
// C'est plus lent car l'accélération est évaluée par pot(x+dx)-pot(x-dx)/2dx ...
template <class Body>  void  evolve_step_Verlet_via_pot(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt, FitParams &params)
{

    Vecteur3D *old_acc = new Vecteur3D[n];
    Vecteur3D *acc = new Vecteur3D[n];

    for (int i = 0; i < n ; ++i)
        old_acc[i] = new_acc_pot(bod, i, fieldB, fieldE, laser, t, params); // calcul l'accélération

    add_acc_N_body_Coulomb(bod, old_acc, n); // add coulombian acceleration due to charged particles
    add_acc_Charged_particles(bod, old_acc, fieldB, fieldE, n); // add Lorentz force for charged particles

    for (int i = 0; i != n ; ++i)
    {
        if (bod[i].get_pos().mag() > LARGE_NUMBER)
            continue; // Si la particule est trop loin on la  laiss là!
        Vecteur3D d_pos = (dt*(bod[i].get_vel()) + (dt*dt*0.5)*old_acc[i]);
        bod[i].inc_pos(d_pos);   // @pos += @vel*dt + old_acc*0.5*dt*dt
    }

    for (int i = 0; i != n ; ++i)
        acc[i] = new_acc_pot(bod, i, fieldB, fieldE, laser, t + dt, params);        // new_acc = acc

    add_acc_N_body_Coulomb(bod, acc, n); // add coulombian acceleration due to charged particles
    add_acc_Charged_particles(bod, acc, fieldB, fieldE, n); // add Lorentz force for charged particles

    for (int i = 0; i != n ; ++i)
    {
        Vecteur3D d_vel = (0.5*dt*(old_acc[i] + acc[i]));
        bod[i].inc_vel(d_vel);   // @vel += (old_acc + new_acc)*0.5*dt
    }

    delete [] acc;
    delete [] old_acc;

}



// Pas d'accélération
template <class Body>  void  evolve_step_Aucun_N_corps(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt)
{
    for (int i = 0; i != n ; ++i)
    {
        if (bod[i].get_pos().mag() > LARGE_NUMBER)
            continue; // Si la particule est trop loin on la  laiss là!
        Vecteur3D d_pos = (dt*(bod[i].get_vel()));
        bod[i].inc_pos(d_pos);   // @pos += @vel*dt
    }
}



/*-----------------------------------------------------------------------------
 *  evolve_step Yoshida6

d = [0.784513610477560e0, 0.235573213359357e0, -1.17767998417887e0,
1.31518632068391e0]
for i in 0..2 do leapfrog(dt*d[i]) end
leapfrog(dt*d[3])
for i in 0..2 do leapfrog(dt*d[2-i]) end

 *-----------------------------------------------------------------------------
*/


// Yoshida 6 using LeapForg with acceleration
template <class Body>
void  evolve_step_Yoshida6_acc(vector <Body> &bod,  const Field &fieldB, const Field &fieldE, const  vector <Laser> &my_laser,
                               const int n, double & t,const double dt)
{
    double  d[4] =
    {
        0.784513610477560e0, 0.235573213359357e0, -1.17767998417887e0, 1.31518632068391e0
    };
    for (int i = 0; i <= 3 ; ++i)                           // for i in 0..2 do leapfrog(dt*d[i])
        evolve_step_Verlet_via_acc(bod, fieldB, fieldE, my_laser, n, t, dt*d[i]);
    for (int i = 0; i <= 2 ; ++i)
        evolve_step_Verlet_via_acc(bod, fieldB, fieldE, my_laser, n, t, dt*d[2-i]); // for i in 0..2 do leapfrog(dt*d[2-i])
}


// Yoshida 6 using LeapForg with potentials
template <class Body>
void  evolve_step_Yoshida6_pot(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
                               const int n, double & t,const double dt, FitParams &params)
{
    double  d[4] =
    {
        0.784513610477560e0, 0.235573213359357e0, -1.17767998417887e0, 1.31518632068391e0
    };
    for (int i = 0; i <= 3 ; ++i)                           // for i in 0..2 do leapfrog(dt*d[i])
        evolve_step_Verlet_via_pot(bod, fieldB, fieldE, laser, n, t, dt*d[i], params);
    for (int i = 0; i <= 2 ; ++i)
        evolve_step_Verlet_via_pot(bod, fieldB, fieldE, laser, n, t, dt*d[2-i], params); // for i in 0..2 do leapfrog(dt*d[2-i])

}





/*-----------------------------------------------------------------------------

SIMULATIONS FOR ION TRAPS – METHODS AND NUMERICAL IMPLEMENTATION

This is the book which summarise several of those ones


Velocity Verlet Method is energy conserving, but not applicable for this
problem, since it does not include a magnetic field.
Sometimes the Gear method is sufficient for this case.
However for uniform magnetic field Eq (57) of the book gives
v_new  = v + dt q/m ( E + (v+v_new)^B/2)  to be compared with the (non magnetic one Eq (53) v_new = v + dt q/m E
r_new = r + v_new dt
Sometimes the notations are with the times steps
v =v_(n-1/2) and v_new = v_(n+1/2) where r=r_n and r_new = r_(n+1)



But, An energy conserving integrator for
magnetic fields has been invented by J. P. Boris in 1970.
has a weak point however;
it requires fine resolution of the Larmor angular frequency
Omega = q B/m typically dt Omega <~ 0.1

Theta[dt] = (qB/m)dt is the rotation angle

Boris–Buneman integration scheme
(cf http://arxiv.org/pdf/1211.3664.pdf or http://www.particleincell.com/blog/2011/vxb-rotation/ or
 Journal of ComputationalPhysics 273 (2014) 255 with space charged)
we note t= Theta[dt/2] = (qB/m)dt/2
or
Theta[dt] = 2 Tan^{-1} (dt/2 qB/m)

Drift: x' - x  = v dt/2

Kick: v- = v + (qE_tot/m) (dt/2) where E_tot = E_ext + E_int (self field)
v' = v- + v-^t
v+ = v- + v' ^(2t/(1+t^2)) is the rotation (axis B) of v- with the angle Theta[dt]. Recall that Sin[theta]=2t/(1+t^2) and 1-Cos[theta] = 2t^2/(1+t^2) where t = tan [theta/2]
v_final = v+ + (qE_tot/m) (dt/2)

Drift:  x_final - x'  = v_final dt/2



Spreiter and Walter: Taylor expansion algorithm works in the opposite limit (dt Omega >>  1).
it is the magnetic field generalization of the velocity-Verlet algorithm.
It is  unstable when dt Omega < 1 (cf Journal of Computational Physics 228.7 (2009)).
 In addition it does not exactly conserve energy, which could be a problem if used in codes where long-term
particle tracking is necessary.

vB(r,v,dt) = v + sin[Theta] (n_B ^ v) + (1-Cos[theta])(n_B ^(n_B ^ v)) is the rotation of v around B (cf Eq (40) of PHYSICAL REVIEW E 77, 066401 (2008)
where I recall that Theta[dt] = (-qB/m)dt is the rotation angle

n_B is the unit vector along B and ^is the cross product
The algorithm is (2a cf PHYSICAL REVIEW E 77, 066401 (2008) but only for magnetic field)
v1 = vB(r,v,dt/2)
r1 = r + v1 dt
v2 = vB(r1,v1,dt/2)

 so it is (I think):
v1 = vB(r,v,dt)
r1 = r + v1 dt
to be compare with Boris (in pure magnetic field)
r' =  r + v dt
v_final = vB(r',v,dt)


We have the generalization is Eq (83) (85) of the book where the space charge is also taken into account
It is also well describe in  "Cooling of highly charged ion Penning Trap (HITRAP)" cf Eq (3.60) Thesis of Giancarlo Maero


Finally  a Cyclotronic integrator has been developped
(cf Journal of Computational Physics 228.7 (2009))
but Cyclotronic integrator does not ensure energy conservation
in non-uniform magnetic field conditions, while the Boris scheme does


 *-----------------------------------------------------------------------------*/



// Boris–Buneman integration scheme (notation of Journal of Computational Physics 273 (2014) 255)
template <class Body>  void  evolve_step_Boris(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt)
{

    Vecteur3D *acc = new Vecteur3D[n];

//  Drift: x' - x  = v dt/2
// In fact I group with the other Drift to make a single one at dt (cf book)
    for (int i = 0; i != n ; ++i)
    {
        if (bod[i].get_pos().mag() > LARGE_NUMBER)
            continue; // Si la particule est trop loin on la  laiss là!
        Vecteur3D d_pos = bod[i].get_vel()*dt; // dt step (because of the merging with the last Drift)
        bod[i].inc_pos(d_pos);   // @pos += @vel*dt + old_acc*0.5*dt*dt
    }

// calcul Zeeman and Stark acceleration
    for (int i = 0; i < n ; ++i)
        acc[i] = new_acc(bod, i, fieldB, fieldE, laser, t);
    add_acc_N_body_Coulomb(bod, acc, n); // add coulombian acceleration due to charged particles


// Kick:
    for (int i = 0; i != n ; ++i)
    {

        double q = bod[i].get_charge();
        double m = bod[i].get_mass();
        Vecteur3D r =  bod[i].get_pos();
        if (r.mag() > LARGE_NUMBER)
            continue;

        Vecteur3D v = bod[i].get_vel();
        Vecteur3D E = fieldE.get_Field(r) ; // E is the local field due to external field (the other fields of space charge are in acc[i])
        Vecteur3D B = fieldB.get_Field(r);
        double B_norm = B.mag();
        Vecteur3D Phi = (B/B_norm)*tan(B_norm*(q/m)*(dt/2.)); // It is t. The simple formula is   Vecteur3D Phi = B*(q/m)*(dt/2.).

        v = v + ((q*E/m) + acc[i])*(dt/2);  // v = v + (qE_tot/m) (dt/2)


        Vecteur3D w = v + v.cross(Phi); // w = v + v * Phi where Phi = dt/2  q B/m
        Vecteur3D s =  2.*Phi/(1.+Phi.mag2());  // s = (2 Phi/(1+Phi^2))

        v = v + w.cross(s) + ((q*E/m) + acc[i])*(dt/2.); // v = v + w * s  + (qE/m) (dt/2)

        bod[i].set_vel(v);
    }

// Drift:  x_final - x'  = v_final dt/2
//    for (int i = 0; i != n ; ++i)
//    {
//        if (bod[i].get_pos().mag() > LARGE_NUMBER) continue; // Si la particule est trop loin on la  laiss là!
//        Vecteur3D d_pos = bod[i].get_vel()*dt/2.;
//        bod[i].inc_pos(d_pos);
//    }

    delete [] acc;
}


// Spreiter and Walter Taylor = magnetized generalization of the velocity-Verlet algorithm.
// I use notation of algorithm 2a of PHYSICAL REVIEW E 77, 066401 (2008). But because it was for B sole. I make the use of the Lorentz transform to have
// B replace by (B - v * E /v^2) so the force q v*B becomes the usual  if v orthogonal to E
template <class Body>  void  evolve_step_Spreiter_Walter(vector <Body> &bod, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        const int n, double & t,const double dt)
{

    Vecteur3D *acc = new Vecteur3D[n];


    for (int i = 0; i != n ; ++i)
    {
        double q = bod[i].get_charge();
        double m = bod[i].get_mass();
        Vecteur3D r =  bod[i].get_pos();
        if (r.mag() > LARGE_NUMBER)
            continue;

        Vecteur3D v = bod[i].get_vel();
        Vecteur3D E = fieldE.get_Field(r) ; // E is the local field due to external field (the other fields of space charge are in acc[i])
        Vecteur3D B = fieldB.get_Field(r);
        double Bnorm = B.mag();
        Vecteur3D n_B;  // unit vector pointing in the direction of the magnetic field
        if (Bnorm !=0)
            n_B = B/Bnorm; // on the other case B is zero so the effect will be nul


        double theta = (-q*B/m)*dt;    // Eq (36) of Phys Rev E


        Vecteur3D v1 = v + sin(theta)*n_B.cross(v) + (1.-cos(theta))*n_B.cross(n_B.cross(v)); // Eq (40)

    }

    delete [] acc;
}




// Pour calcul a_i= Force/masse = -Grad(E_pot)/mass;
// Calcul et met a jour l'accélération en utilisant le potentiel
template <class Body>  Vecteur3D new_acc_pot(vector <Body> &bod, const int i, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params)
{
    if (params.LocateParam("Choix_algorithme_N_corps")->val == Verlet_pot)
        return new_acc_pot_second_order(bod, i, fieldB, fieldE, laser, t,  params);

    if (params.LocateParam("Choix_algorithme_N_corps")->val == Verlet_pot_gradient_high_order)
        return new_acc_pot_fourth_order(bod, i, fieldB, fieldE, laser, t,  params);

    return Vecteur3D(0.,0.,0.);
}


// Grad calculé par (E(r+eps)-E(r-eps))/(2*eps) pour chaque axe
template <class Body>  Vecteur3D new_acc_pot_second_order(vector <Body> &bod, const int i, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params)
{
    Body &my_bod = bod[i]; // Body my_bod = bod[i]; is wrong. We need to take a real copy of bod[i] in order to modify its position

    Vecteur3D ai(0.,0.,0.);       // Accélération de la particule i (pour masse = 1, voir à la fin)

    double epsilon = params.LocateParam("choix_epsilon")->val; // paramètre petit pour le calcul du gradient. Erreur en epsilon^2

    Vecteur3D epsx(epsilon,0.,0.);

    my_bod.inc_pos(epsx); // position + epsx
    double Eplus = new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(-2.*epsx); // position - epsx
    double Eminus= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(epsx);  // remet la particule à sa place
    double aix=-((Eplus-Eminus)/(2.*epsilon));
    ai += Vecteur3D(aix,0.,0.); // Force selon x


    Vecteur3D epsy(0.,epsilon,0.);
    my_bod.inc_pos(epsy); // position + epsy
    Eplus = new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(-2.*epsy); // position - epsy
    Eminus= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(epsy);  // remet la particule à sa place
    ai += Vecteur3D(0.,-((Eplus-Eminus)/(2.*epsilon)),0.); // Force selon y


    Vecteur3D epsz(0.,0.,epsilon);
    my_bod.inc_pos(epsz); // position + epsz
    Eplus = new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(-2.*epsz); // position - epsz
    Eminus= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(epsz);  // remet la particule à sa place
    ai += Vecteur3D(0.,0.,-((Eplus-Eminus)/(2.*epsilon))); // Force selon z

    Vecteur3D acc;
    acc = gravity  + ai/my_bod.get_mass();

    my_bod.set_acc(acc);
    return acc; // ACCELERATION
}

// Calcul du gradient par
// (1/12 E(r-2 eps)-2/3 E(r-eps) + 2/3 E(r+eps ) - 1/12 E(r-2 eps))/(eps) pour chaque axe
// cf http://en.wikipedia.org/wiki/Finite_difference_coefficients
// http://en.wikipedia.org/wiki/Five-point_stencil
template <class Body>  Vecteur3D new_acc_pot_fourth_order(vector <Body> &bod, const int i, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, FitParams &params)
{
    Body &my_bod = bod[i];

    Vecteur3D ai(0.,0.,0.);       // Accélération de la particule i (pour masse = 1, voir à la fin)

    double epsilon = params.LocateParam("choix_epsilon")->val; // paramètre petit pour le calcul du gradient. Erreur en epsilon^2

    double E2,E1,Em1,Em2;

    Vecteur3D eps(epsilon,0.,0.);

    my_bod.inc_pos(-2.*eps); // position - 2 eps
    Em2 = new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(eps); // position - eps
    Em1= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(2.*eps);  //  position + eps
    E1= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(eps);  //  position + 2*eps
    E2= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(-2.*eps);  //  position
    double ac=((1./12.)*Em2 - (2./3.)*Em1 + (2./3.)*E1 - (1./12.)*E2)/(epsilon);
    ai += Vecteur3D(-ac,0.,0.); // Force selon x (- gradient)


    eps = Vecteur3D(0.,epsilon,0.);

    my_bod.inc_pos(-2.*eps); // position - 2 eps
    Em2 = new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(eps); // position - eps
    Em1= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(2.*eps);  //  position + eps
    E1= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(eps);  //  position + 2*eps
    E2= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(-2.*eps);  //  position
    ac=((1./12.)*Em2 - (2./3.)*Em1 + (2./3.)*E1 - (1./12.)*E2)/(epsilon);
    ai += Vecteur3D(0.,-ac,0.);  // Force selon y


    eps = Vecteur3D(0.,0.,epsilon);

    my_bod.inc_pos(-2.*eps); // position - 2 eps
    Em2 = new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(eps); // position - eps
    Em1= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(2.*eps);  //  position + eps
    E1= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(eps);  //  position + 2*eps
    E2= new_pot(bod, i, fieldB, fieldE, laser, t,  params);
    my_bod.inc_pos(-2.*eps);  //  position
    ac=((1./12.)*Em2 - (2./3.)*Em1 + (2./3.)*E1 - (1./12.)*E2)/(epsilon);
    ai += Vecteur3D(0.,0.,-ac);  // Force selon y



    Vecteur3D acc;
    acc = gravity  + ai/my_bod.get_mass();

    my_bod.set_acc(acc);
    return acc; // ACCELERATION
}


/*-----------------------------------------------------------------------------
 *  get_acc_jerk_coll  --  calculates accelerations and jerks, and as side
 *                         effect also calculates the time scale coll_time for
 *                         significant changes in local configurations to occur.
 *
 *  force on i due to j/m_i  =  a_ji = q_i q_j r_ji/(4 pi epsilon_0 ||r_ji||^3) where r_ji = r_i - r_j (vector from j to i)
 *
 * Modifications by Walter (but adapted to Coulombian case by Daniel Comparat) to improve performance:
 *
 *  - using of vector
 *  - (re-) using of register variables.
 *  - avoiding divisions
 *  - computing forces rather than accelerations:
 *    the loop to compute accelerations looks like this (with obvious notation):
 *
 *    for(i=0; i!=n; ++i)
 *        a[i] = 0.                                 // reset accelerations to 0
 *    for(i=0; i!=n; ++i) {
 *        qi = q[i];                                // get charge and pos of ith
 *        ri = r[i];                                // body into register
 *        ai = 0.;                                  // to hold force due to j>i
 *        for(j=i+1; j!=n; ++j) {                   // loop pairs j>i
 *           rji  =  ri- r[j] ;                      // distance vector
 *           rji *= qi*q[j] / |rji|^3               // mutual force
 *           ai  += rji;                            // add: force_i due to j>i
 *           a[j]-= rji;                            // add: force_j due to i<j
 *        }
 *        a[i] = (a[i] + ai)/m[i];                    // acceleration
 *    }
 *
 *-----------------------------------------------------------------------------
 */


// Calculate the Coulombian acceleration (adapted from Walter Dehnen)
// Is added to accell
template <class Body>  void  add_acc_N_body_Coulomb(const vector <Body> &bod, Vecteur3D *accell, const int n)
{

    Vecteur3D *acc = new Vecteur3D[n];

    for (int i = 0; i < n ; i++)
        acc[i] = Vecteur3D(0.,0.,0.);

    for (int i = 0; i != n ; ++i)
    {
        const double qi = bod[i].get_charge();                    // qi
        if (qi == 0.)
            continue;

        const Vecteur3D ri = bod[i].get_pos();                  // ri
        register Vecteur3D ai = Vecteur3D(0.,0.,0.);        // ai due to j>i
        for (int j = i+1; j != n ; ++j)
        {
            // pre-compute some auxiliary quantites
            Vecteur3D rji= -bod[j].get_pos();
            rji += ri; // rji = distance vector
            register double
            r2    = rji.mag2(),                   // rji^2
            pr2   = 1./r2;                       // 1 / rji^2

            // add the {j (i)} contribution to the {i (j)} values of
            // force
            pr2 *= bod[j].get_charge()*qi*sqrt(pr2);  // mi*mj / |rji|^3  (for Newton, here qi qj/ 4pi epsilon_0 r^3)
            rji *= pr2;                          // mutual force = mi*mj*da
            ai   += rji;               // force     at i due to j>i
            acc[j] -= rji;               // force     at j due to i<j
        }
        // add forces due to j>i and divide by mass to get accelerations
        register double tmp = kCoulomb/bod[i].get_mass(); // kCoulomb=1/(4 pi epsilon_0)
        acc[i]   += ai;
        acc[i]   *= tmp;
    }

    for (int i = 0; i < n ; i++)
        accell[i] += acc[i];

    delete [] acc;

}


// Calculate the Coulombian field for the sum_j=1^n not i
template <class Body>  Vecteur3D  get_coulomb_field(const vector <Body> &bod, const int i, const Vecteur3D &ri)
{
    Vecteur3D Ei=Vecteur3D(0.,0.,0.);

    for (int j = 0; j != bod.size() ; ++j)
    {
        // pre-compute some auxiliary quantites
        Vecteur3D rji= -bod[j].get_pos(); // -rj
        rji += ri; // rji =  ri-rj = distance vector
        double r2    = rji.mag2();                   // rji^2
        if (r2 < VERY_SMALL_NUMBER)
            continue; // this mean that rji=0 and so i = j

            double pr2   = 1./r2;                       // 1 / rji^2

        pr2 *= bod[j].get_charge()*sqrt(pr2);  //  qj/  r^3
        rji *= pr2;
        Ei += rji;          //  qj vec(rji)/rji^3
    }

    return kCoulomb*Ei; // kCoulomb=1/(4 pi epsilon_0)

}




// Calculate the Coulombian potential for the sum_j=1^n not i qj/(r_ij)/4 pi epsilon_0
// Is used only for statistics to check the energy conservation. So it is not fully optimized
template <class Body>  double  get_coulomb_potential(const vector <Body> &bod, const Vecteur3D &ri, const int n)
{

    double pot_i = 0.;

    for (int j = 0; j != n ; ++j)
    {
        double rji2 = (ri - bod[j].get_pos()).mag2(); // r_ij^2
        if (rji2 < VERY_SMALL_NUMBER)
            continue; // this mean that rji=0 and so i = j
        pot_i += kCoulomb*bod[j].get_charge()/sqrt(rji2); // kCoulomb=1/(4 pi epsilon_0)
    }
    return pot_i;

}


// add Lorentz force for charged particles
template <class Body>  void  add_acc_Charged_particles(vector <Body> &bod, Vecteur3D *accell, const Field &fieldB, const Field &fieldE, const int n)
{
    for (int i = 0; i < n ; i++)
    {
        double q = bod[i].get_charge();
        double m = bod[i].get_mass();
        Vecteur3D r =  bod[i].get_pos();
        Vecteur3D v = bod[i].get_vel();
        Vecteur3D E = fieldE.get_Field(r);
        Vecteur3D B = fieldB.get_Field(r);
        Vecteur3D crossedvalue = v.cross(B);
        Vecteur3D force = q*(E + crossedvalue)/m;
        accell[i] += q*(E + v.cross(B))/m; // Lorentz force
    }

}



// kinetic energy
template <class Body> double get_kin(Body &bod)
{
    double  ekin = 0.;
    const double m = bod.get_mass();                    // mi

    const Vecteur3D v = (bod.get_vel());                    // vi
    ekin += m * v*v;

    ekin *= 0.5;

    return ekin;
}

// temps d'évolution typique des forces (ou potentiels) ~ 0.1 waist/vitesse. Le 0.1 est empirique et dépend des paramètres comme la largeur de la transition
//  A améliorer avec les champs, différents lasers (pas le seul premier) et Temp_estime = Temp ?...
template <class Body> double t_evol_ext(vector <Body> &bod,  const Field &fieldB, const Field &fieldE, const  vector <Laser> &my_laser, const double Temp_estime)
{
    double t_evol;
    double waist =  my_laser[0].get_waist().mag(); // Waist (Attention juste du premier laser, en supposant qu'il représente bien les autres waits)
    double vitesse_therm = sqrt(kB*Temp_estime/bod[0].get_mass()); // vitesse thermique en m/s 1/2 m vx^2 = 1/2 kB T
    t_evol = min(0.1*waist/vitesse_therm,0.001); // au max 1 ms
    return t_evol;
}



#endif


