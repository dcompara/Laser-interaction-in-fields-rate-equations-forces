/*
  Name:
  Copyright:
  Author:
  Date: 23/10/06 10:01
  Description: Paramètres (température, positions, dispersion, moyenne des molécules) (utilisé dans Affichage_Mol.h)


*/

#include "Stat_Molecule.h"
#include <math.h>   // Pour log


// sort le tableau nexc du nombre de molécule par niveau nexc[n° Level]
void stat_molecule_nb_par_etat(const vector <Molecule> &Mol, Internal_state *Levels, int *nexc, const int Nb_Mol, const int Nb_state, FitParams &params)
{
    for (int i=0; i!= Nb_state; i++)  // Nombre de particules excitées dans chaque état exc
        nexc[i]=0;

    for (int i=0; i!= Nb_state; i++)
    {
        for (int j=0; j!= Nb_Mol; j++)
        {
            if (Mol[j].is_equal(Levels[i])) nexc[i]++;  // Calcul le nombre de particules excitées dans chaque état j
        }
    }

    return;
}

// Paramètres statistiques d'une assemblée de Molecule;
// Le dernier niveau stat_Mol[nb_Levels] donne la statistique globale (tout niveaux confondus)
void stat_molecule(const vector <Molecule> &Mol,  vector <Stat> &stat_Mol, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const int Nb_Mol, const int Nb_state, FitParams &params)
{
    stat_Mol.clear();
    Stat stat_Mol_current;
    for (int i=0; i!= Nb_state; i++)
    {
        stat_molecule_un_niveau(Mol,stat_Mol_current, Level, fieldB, fieldE, i, Nb_Mol, params); // Statistique pour le niveau i
        stat_Mol.push_back(stat_Mol_current);
    }
    stat_molecule_un_niveau(Mol, stat_Mol_current, Level, fieldB, fieldE, all , Nb_Mol, params); // Statistiques pour toutes les molécules
    stat_Mol.push_back(stat_Mol_current);
}

// Paramètres statistiques d'une assemblée de Molecule; Retourne le nb de molécules n dans l'état étudié
int stat_molecule_un_niveau(const vector <Molecule> &Mol, Stat &stat, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const int num_niveau_etudie, const int Nb_Mol, FitParams &params)
{

    int num_manifold_not_studied = (int) params.LocateParam("num_manifold_not_studied")->val;


    /***************** Liste des molécules de l'état choisi leurs positions et vitesse  *************************/

    vector <Molecule> liste_Adresses_Mol_dans_niveau; //Liste des molécules dans ce niveau

    for (int i=0; i!= Nb_Mol; i++)
        if ( (num_niveau_etudie == all || Mol[i].is_equal(Level[num_niveau_etudie])) &&  Mol[i].exc != num_manifold_not_studied) // NE PREND QUE LES MOL DANS l'état spécifique
            liste_Adresses_Mol_dans_niveau.push_back(Mol[i]);


    int n = liste_Adresses_Mol_dans_niveau.size();
    stat_molecule_form_list(Mol, liste_Adresses_Mol_dans_niveau, stat, Level, fieldB, fieldE,  Nb_Mol, params);

// erase the vector:
    liste_Adresses_Mol_dans_niveau.clear();

    return n;
}


// Paramètres statistiques d'une assemblée de Molecule; Retourne le nb de molécules n dans l'état étudié
int stat_molecule_form_list(const vector <Molecule> &Mol,  vector <Molecule> &liste_Adresses_Mol_dans_niveau, Stat &stat, const vector <Internal_state> &Level,  const Field &fieldB, const Field &fieldE, const int Nb_Mol, FitParams &params)
{
    int n = liste_Adresses_Mol_dans_niveau.size();
    stat.population = n;
    if (n == 0) return n; // Ne calcul rien si aucune molécules dans l'état spécifique

    double **pos = new double*[3]; // Tableau des positions
    double **vel = new double*[3]; // Tableau des vitesses
    double **E_cinetique_Temp = new double*[3]; // Tableau des energies cinétiques (sans le mouvement du centre de masse)
    for (int j=0; j<3; j++)   // 3 axes: x, y ou z
    {
        pos[j] = new double[n];
        vel[j] = new double[n];
        E_cinetique_Temp[j] = new double[n];
    }

    stat.E_cin = 0.;   // Energie totale
    stat.E_pot = 0.;

    int nb =0;
    double mass = 0.; // Mass on the particles. But BE CAREFUL if several mass and global motion, the temperature is not well defined!)

    for ( vector <Molecule>::iterator i=liste_Adresses_Mol_dans_niveau.begin(); i != liste_Adresses_Mol_dans_niveau.end(); ++i) // boucle sur les molécules de l'état considéré
    {
        Molecule Mol_i = *i;
        mass = Mol_i.get_mass();

        for (int j=0; j<=2; j++)  // 3 axes: x, y ou z
        {
            pos[j][nb] =  Mol_i.get_pos().get(j);
            vel[j][nb] =  Mol_i.get_vel().get(j);
            E_cinetique_Temp[j][nb] = 0.5* mass * vel[j][nb]* vel[j][nb] ; // 1/2 m v_j^2
        }

        stat.E_pot += get_pot(Mol_i, fieldB, fieldE);
        stat.E_pot += 0.5*Mol_i.get_charge()*get_coulomb_potential(Mol, Mol_i.get_pos(), Nb_Mol); // add colombian interaction will be 1/2 sum_i qi sum_j=1^n not i qj/(r_ij)/4 pi epsilon_0
        // HERE we have to sum over all level because the coulomb field exist for all states

        stat.E_cin += get_kin(Mol_i);
        nb++;
    }



    /*****************  Statistiques sur leurs positions et vitesse  *************************/


    for (int j=0; j<=2; j++)  // 3 axes: x, y ou z
    {
        stat.mean_pos.set(j,gsl_stats_mean (pos[j], 1 ,  n));
        stat.mean_v.set(j,gsl_stats_mean (vel[j], 1 ,  n));
        stat.sigma_pos.set(j,sqrt(gsl_stats_variance_m (pos[j], 1 ,  n, stat.mean_pos.get(j)))); // Sigma = sqrt(variance)
        stat.sigma_v.set(j,sqrt(gsl_stats_variance_m (vel[j], 1 ,  n, stat.mean_v.get(j))));
    }


    // We remove the center of mass velocity to calculate 1/2 kB T is average of 1/2 m (v-<v>)^2
    nb =0;
    for (vector <Molecule>::iterator i=liste_Adresses_Mol_dans_niveau.begin(); i != liste_Adresses_Mol_dans_niveau.end(); ++i) // boucle sur les molécules de l'état considéré
    {
        for (int j=0; j<=2; j++)  // 3 axes: x, y ou z
        {
            double mean_v = stat.mean_v.get(j);
            E_cinetique_Temp[j][nb] +=  0.5 * mass * ( -2.*mean_v*vel[j][nb] + mean_v*mean_v); // 1/2 m (v-<v>)^2 = 1/2 m (v^2 - 2 v <v> + <v>^2)
        }
        nb++;
    }


//    for (int j=0; j<=2; j++)  // 3 axes: x, y ou z
//    {
//        stat.Temp_cin.set(j,2.*gsl_stats_mean (E_cinetique_Temp[j], 1 ,  n)/k_Boltzmann); //  1/2 kB T = 1/2 m <v - <v> >^2
//        qsort(E_cinetique_Temp[j], n, sizeof(double), &compare_doubles); // Ordonne E_cinetique[j]
//        double median_j = gsl_stats_median_from_sorted_data (E_cinetique_Temp[j], 1, n);
//        stat.Temp_1D_50.set(j,median_j /(0.227468*k_Boltzmann)); // 1D Boltzmann   proba density distribution propto 1/sqrt(pi kB T E) e^{-  E/ KB T} ;  E =1/2 m v_i^2 --> median = 0.227468 k T
//    }


    for (int j=0; j<=2; j++)  // 3 axes: x, y ou z
    {
        stat.Temp_cin.set(j,2.*gsl_stats_mean (E_cinetique_Temp[j], 1 ,  n)/k_Boltzmann); //  1/2 kB T = 1/2 m <v - <v> >^2
        qsort(vel[j], n, sizeof(double), &compare_doubles); // Ordonne vel[j]
        double v25 = gsl_stats_quantile_from_sorted_data (vel[j], 1, n, 0.25); // v25 = velocity of the molecule that have the 25th velocity
        double v75 = gsl_stats_quantile_from_sorted_data (vel[j], 1, n, 0.75); // v75 = velocity of the molecule that have the 25th velocity
        double dv = (v75-v25)/(2.*0.67449); // 0.67449 is the quantil at 75% of a Normal (Gaussian distribution) e(-x^2/2). THis  is =\Phi^{-1}(0.75)   where  \Phi(x)= P[X\leq x], lorsque X suit la loi normale centrée réduite \scriptstyle \mathcal N(0,1).
        stat.Temp_1D_50.set(j,dv*dv*mass/k_Boltzmann);
    }


    stat.Temp_3D_50 = stat.Temp_1D_50.mag()/sqrt(3.); // T_3D^2 = (Tx^2+Ty^2+Tz^2)/3




    /*****************  Delete table  *************************/

    for (int j=0; j<3; j++)
    {
        delete [] pos[j];
        delete [] vel[j];
        delete [] E_cinetique_Temp[j];
    }
    delete [] pos;
    delete [] vel;
    delete [] E_cinetique_Temp;

    return n;
}


// Statistique sur les trajectoires
// distance Min, moyenne (somme) et max
void stat_traj_molecule(const vector <Molecule> &Mol, double & dist_min, double & dist_sum, double & dist_max, int & nb_call, const int Nb_Mol)
{
    for (int i=0; i!= Nb_Mol; i++)
    {
        double dist_mol = (Mol[i].get_pos()).mag(); // Distance de la molécule au centre (0,0,0)

        if ( dist_mol < dist_min)  // position
            dist_min = dist_mol;

        if ( dist_mol > dist_max)  // position
            dist_max = dist_mol;

        dist_sum +=  dist_mol;
    }

    nb_call++;

    cout << "min (mm) " << dist_min/mm << "avg " << dist_sum/nb_call/Nb_Mol/mm << " max " << dist_max/mm << endl;

    return;

}
