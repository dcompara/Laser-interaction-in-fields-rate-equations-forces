#include "scaling_velocities.h"
#include "Atome.h"
#include "Internal_state.h"
#include "Vecteur3D.h"
#include "Laser.h"                       // Classe Laser
#include "Molecule.h"
#include "one_body.h"
#include "Stat_Molecule.h"

// fonction qui applique la constante lambda(t) aux vitesses(t).
//INPUTS: -coupling_efficiency is in percentage.
//OUTPUT: -file_scal: the aimed temperature and the rescaled one are outputs.
void scaling_velocities( const double t, vector <Molecule> &Mol, const int nb_mol,const vector <Internal_state> &Level, const Field &fieldB,const Field &fieldE, FitParams &params,double coupling_efficiency, ofstream & file_scal)
{


    Stat stat_Mol; //for having the T temperature of the cloud (Ti= stat.Temp_1D_100)
    vector <Molecule> liste_Adresses_Mol; //List of Molecule of first type
    for (int i=0; i!= nb_mol; i++)
    {
        liste_Adresses_Mol.push_back(Mol[i]);
    }

    stat_molecule_form_list(Mol, liste_Adresses_Mol, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);

//Parameters of the scaling:
    //temperature that we want to reach :
    double Temp_ini_x  = params.LocateParam(string("Temp_ini_x[0]"))->val;
    double Temp_ini_y  = params.LocateParam(string("Temp_ini_y[0]"))->val;
    double Temp_ini_z  = params.LocateParam(string("Temp_ini_z[0]"))->val;
    Vecteur3D Temp_0(Temp_ini_x,Temp_ini_y,Temp_ini_z);  // (Tx0,Ty0,Tz0)


    Vecteur3D Temp_instant; //temperature instantanée du nuage, à temps t

    Temp_instant= stat_Mol.Temp_1D_50;

// coefficient LAMBDA for rescaling
    Vecteur3D ratio_T;
    Vecteur3D factor1;
    double coeff=coupling_efficiency;

    ratio_T = Temp_0/Temp_instant;
    ratio_T= ratio_T-1.;
    factor1= coeff*ratio_T;
    factor1= factor1+1.;
    factor1= abso(factor1);
    //BERENDSEN constant:
    Vecteur3D Lambda=racine(factor1);


    vector <Molecule> liste_Adresses_Mol_dans_niveau; //List of Molecule of first type
    for (int i=0; i!= nb_mol; i++)
    {
        Vecteur3D tot;
        tot= Hadamard(Mol[i].get_vel(),Lambda);
        Mol[i].set_vel(tot); //on multiplie les vitesses selon x,y et z par Lambda_x,y,z(t)
        liste_Adresses_Mol_dans_niveau.push_back(Mol[i]);
    }

    stat_molecule_form_list(Mol, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);

    Temp_instant = stat_Mol.Temp_1D_50; // on recalcule T

    file_scal << t << " " << Temp_0<< " " << Temp_instant << " " ;
    file_scal << endl;

    liste_Adresses_Mol.clear(); // erase the vector:


}
//DEUX FONCTIONS POUR PERMETTRE DE RESCALER A UN TEMPS t+dt:
// v'(t+dt)=v(t+dt)*Lambda(t)
//Fonction qui calcule lambda à t (il faut appliquer le facteur de rescalng Lambda aux vitesses à t+dt_scal)
Vecteur3D calcul_lambda( const double t, vector <Molecule> &Mol, const int nb_mol,const vector <Internal_state> &Level, const Field &fieldB,const Field &fieldE, FitParams &params,double coupling_efficiency, ofstream & file_scal)
{

    Stat stat_Mol; //for having the T temperature of the cloud (Ti= stat.Temp_1D_100)
    vector <Molecule> liste_Adresses_Mol; //List of Molecule of first type
    for (int i=0; i!= nb_mol; i++)
    {
        liste_Adresses_Mol.push_back(Mol[i]);

    }
    stat_molecule_form_list(Mol, liste_Adresses_Mol, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);

//Parameters of the scaling:
    //temperature that we want to reach :
    double Temp_ini_x  = params.LocateParam(string("Temp_ini_x[0]"))->val;
    double Temp_ini_y  = params.LocateParam(string("Temp_ini_y[0]"))->val;
    double Temp_ini_z  = params.LocateParam(string("Temp_ini_z[0]"))->val;
    Vecteur3D Temp_0(Temp_ini_x,Temp_ini_y,Temp_ini_z);  // (Tx0,Ty0,Tz0)


    Vecteur3D Temp_instant; //temperature instantanée du nuage
    Temp_instant= stat_Mol.Temp_1D_50;
cout << " Temp_1D_50 avant " << Temp_instant << endl;
    // coefficient LAMBDA for rescaling
    Vecteur3D ratio_T;
    Vecteur3D factor1;
    double coeff=coupling_efficiency;

    ratio_T = Temp_0/Temp_instant;
    ratio_T= ratio_T-1.;
    factor1= coeff*ratio_T;
    factor1= factor1+1.;
    factor1= abso(factor1);
    //BERENDSEN constant:
    Vecteur3D Lambda=racine(factor1);


    return Lambda;
}
//fonction qui applique Le coefficient Lambda aux vitesses, à un temps t. (il faut l'appliquer à t+dt_scal).
void rescaling_velocities_after_dt ( const double t, vector <Molecule> &Mol, const int nb_mol,const vector <Internal_state> &Level, const Field &fieldB,const Field &fieldE, FitParams &params,double coupling_efficiency, ofstream & file_scal, Vecteur3D Lambda)
{
    Stat stat_Mol; //for having the T temperature of the cloud (Ti= stat.Temp_1D_100)
 vector <Molecule> liste_Adresses_Mol_dans_niveau; //List of Molecule of first type
    for (int i=0; i!= nb_mol; i++)
    {
        Vecteur3D tot;
        tot= Hadamard(Mol[i].get_vel(),Lambda);
        Mol[i].set_vel(tot);
        liste_Adresses_Mol_dans_niveau.push_back(Mol[i]);
    }

    stat_molecule_form_list(Mol, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);

Vecteur3D    Temp_instant = stat_Mol.Temp_1D_50; // on recalcule T
    cout << " Temp_1D_50 apres " << Temp_instant << endl;
    file_scal << t << " " << Temp_instant << " " ;
    file_scal << endl;

}
