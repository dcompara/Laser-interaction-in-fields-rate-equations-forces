/*
 Fichier de sortie des données moléculaires: transition, probablilitées ...
*/

#include "sortie_donnees.h"


//Sauve l'état du générateur de nombre aléatoire pour ne jamais reprendre le même
void save_random_generator(gsl_rng * r, const char *nom_file)
{
    FILE *fp;
    fp=fopen(nom_file, "wb");    // Fichier contenant les données du générateur de nombre aléatoire
    int ran_gen_result = gsl_rng_fwrite (fp, r); // This function writes the random number state of the random number generator r
    if (ran_gen_result == GSL_EFAILED)
        cout << " Problem writing to the file: Data/random_gen.txt" << endl;
    fclose(fp);
}



void Sortie_donnee_example(ofstream & file_out,  vector <Molecule> &Mol,  vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params,  DataCards &data, const int number_photons)
{
    //set_pot_all_mol(Mol, fieldB, fieldE, laser, t, nb_mol, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour uen sortie
// ATTENTIION THIS DOES NOT WORK FOR THE POTENTIALS


// SOrtie des paramètres scannés
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                file_out << p.name ;
                file_out << p.val_t0 << " " ;
            }
        }
    }

    file_out<< setprecision(8);

    const int num_niveau_etudie = (int) params.LocateParam("num_niveau_etudie")->val; // numéro du niveau étudié pour faire des stats. -1 pour toutes les molécules

    Stat stat_Mol;
    stat_molecule_un_niveau(Mol, stat_Mol, Level,  fieldB, fieldE, num_niveau_etudie, Mol.size(), params);


    for (int i = 0; i < nb_mol; i++)
    {
        set_pot_mol(Mol, i, fieldB, fieldE, laser, t, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour une sortie
        file_out  << t << " ";
        file_out << stat_Mol.Temp_3D_50 << " ";
        file_out  << stat_Mol.population << " ";
        file_out  << Mol[i].get_pos() << " ";
        file_out  << Mol[i].get_vel() << " ";
        file_out  << Mol[i].deg_number << " ";
        cout << endl;
    }
    return;
}



void Sortie_donnee(ofstream & file_out,  vector <Molecule> &Mol,  vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params,  DataCards &data, const int number_photons)
{
    //set_pot_all_mol(Mol, fieldB, fieldE, laser, t, nb_mol, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour uen sortie
// ATTENTIION THIS DOES NOT WORK FOR THE POTENTIALS


// SOrtie des paramètres scannés
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
//                 file_out << p.name ;
//                 file_out << p.val_t0 << " " ;
            }
        }
    }

    file_out<< setprecision(8);

    const int i = (int) params.LocateParam("num_niveau_etudie")->val; // numéro du niveau étudié pour faire des stats. -1 pour toutes les molécules

    Stat stat_Mol;
    stat_molecule_un_niveau(Mol, stat_Mol, Level,  fieldB, fieldE, i, Mol.size(), params);


    for (int i = 0; i < nb_mol; i++)
    {
        set_pot_mol(Mol, i, fieldB, fieldE, laser, t, params); //Met à jour de tous les potentiels (gravité, PAS dipolaire, magnétique, electrique, ...) avec la nouvelle position pour une sortie
    }

//    double Temp_ini_z  = params.LocateParam("Temp_ini_z[0]")->val;
//    const double vitesse_therm_z = sqrt(kB*Temp_ini_z/Mol[0].get_mass()); // vitesse thermique en m/s 1/2 m vx^2 = 1/2 kB T



//    file_out  << t << " ";
//    file_out << number_photons << " ";
//    file_out << stat_Mol.Temp_3D_50 << " ";
//    file_out  << stat_Mol.population << " ";
//    file_out  << stat_Mol.Temp_1D_50.z() << " ";


//  file_out  << (stat_Mol.E_pot/kB)/mK/1.5/nb_mol << "  ";
    //  file_out  << (stat_Mol.E_cin/kB)/mK/1.5/nb_mol << "  ";
    //  file_out << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";

//    cout  << (stat_Mol.E_pot/kB)/mK/1.5/nb_mol << "  ";
//    cout  << (stat_Mol.E_cin/kB)/mK/1.5/nb_mol << "  ";
//    cout << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";
//
//    cout << endl;

    for (int i = 0; i < nb_mol; i++)
    {
//        double z_init = 3.*(i-nb_mol/2.)*(params.LocateParam("size_z[0]")->val)/nb_mol;
//        double vz_init = 3.*(i-nb_mol/2.)*vitesse_therm_z/nb_mol;
//        file_out  << z_init << " ";
//        file_out  << Mol[i].get_pos().z() << " ";
//        file_out  << vz_init << " ";
//        file_out  << Mol[i].get_vel().z() << " ";
//        file_out << t << " ";
//        file_out  << Mol[i].get_pos() << " ";
//        file_out  << Mol[i].get_vel() << " ";
//        file_out  << Mol[i].deg_number << " ";
//    file_out << endl;


        /*****   CALCUL of parameters for the dipoles (in Debye)or diagonalization *******/


        MatrixXcd d[3];
        SelfAdjointEigenSolver<MatrixXcd> es; // eigenstates and eigenvalues
        Diagonalization(Level, Mol[i], fieldB, fieldE, params, es, d);

        Vecteur3D v;
        v = Mol[i].get_vel();
        Vecteur3D Bfield,Efield;
        Bfield= fieldB.get_Field(Mol[i].get_pos());
        Efield= fieldE.get_Field(Mol[i].get_pos());
        double B = Bfield.mag();
        double E = Efield.mag();
        double v_perp= (v.cross(Bfield)).mag()/B;

        for (int j=0; j< (int) Level.size(); j++) //  we scan over the levels to calculate the parameter
        {
            double tripletness = 0.; //This is the parameter we want to calculate (here the triplet character)

            for (int j0=0; j0< (int)  Level.size(); j0++) //  | j> =  sum_|j>_O   0_<j | j>  |j>_0  so we scan over |j>_0 hat is the order in the Level file
                // 0_<j | j>  is given by   es.eigenvectors()(j0,j) . This is the coefficient of the |j> level (ordered in Energy) on the |j>_0 Level (the order in the Level file). We round it to 100%
            {
                // file_out << B << " " << v_perp << " " << i << " " << j  << "  " << j0 << " " << Level[j0].two_M << " " << abs(round(100.*es.eigenvectors()(j0,j)))/100 << endl;
                if (Level[j0].v == 2) // If the state is triplet (2S+1=3 so S=1 coded in v) we look on the decomposition, |0_<i | i>|^2 , and sum them
                {
                    tripletness += abs( pow((es.eigenvectors()(j0,j)),2) ); // sum_|triple, j>_O   |0_<j | j>|^2.
                }
            }
            // file_out << endl;
            // PARAMETER THAT GIVE THE TRIPLETNESS OF THE STATE //

// Level[j].write_Level_B(file_out);




             file_out << B << " " << E << " " << v_perp << " " << j << " " << Level[j].Energy_cm << " " << tripletness << endl;
        }
    }



    /******  The new dipole are given by d[polar] = evec^dag.d0[polar].evec that is
    d[n_polar+1] =  (es.eigenvectors().adjoint())*d0[n_polar+1]*(es.eigenvectors());
    With d0[q+1]_ij = 0_<i | d^(q) | j>_0
    d[q+1]_ij = <i | d^(q) | j> = sum_|j>_O    <i| d^(q) |j>_0   0_<j | j>
    ********/




    /**  Stat for specific states v  **/

//    vector <Molecule> liste_Adresses_Mol_dans_niveau; //List of Molecule of first type
//    int nb_Mol_in_this_state = 0;
//    for (int i=0; i!= nb_mol; i++)
//        if (Mol[i].v == 0) // Molecules not in v=1,2,3 of X state
//        {
//            liste_Adresses_Mol_dans_niveau.push_back(Mol[i]);
//            nb_Mol_in_this_state++;
//        }
//    stat_molecule_form_list(Mol, liste_Adresses_Mol_dans_niveau, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);
//    file_out << nb_Mol_in_this_state << " ";
//    file_out  << stat_Mol.Temp_1D_50.z() << " ";
//    file_out  << stat_Mol.population << " ";
//    file_out << stat_Mol.E_pot/kB/mK/1.5/nb_mol << "  ";
//    file_out << stat_Mol.E_cin/kB/mK/1.5/nb_mol << "  ";
//    file_out << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";
//    liste_Adresses_Mol_dans_niveau.clear(); // erase the vector:
//
//    cout << " t " << t << " photons = " << (double) number_photons ;
//    cout << " Epot " << stat_Mol.E_pot/kB/mK/1.5/nb_mol << "  ";
//    cout << " Ekin " <<  stat_Mol.E_cin/kB/mK/1.5/nb_mol << "  ";
//    cout << " E " << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";
//    cout << endl;

    /**  Stat for specific states of "best" molecules in the sens of position  **/

    /****

    vector <Molecule> liste_Adresses_Mol; //List of Molecule of this type
    //  double size_limite = params.LocateParam("size_x[0]")->val;
    double size_limite = 0.02;
    int nb_Mol_in_this_state = 0;
    double niveau_moyen = 0.;
    for (int i=0; i!= nb_mol; i++)
        if (Mol[i].get_pos().mag() < size_limite) // Molecules  within initial size (in x)
        {
            liste_Adresses_Mol.push_back(Mol[i]);
            nb_Mol_in_this_state++;
            niveau_moyen += Mol[i].exc;
        }
    niveau_moyen = niveau_moyen/nb_Mol_in_this_state;
    stat_molecule_form_list(Mol, liste_Adresses_Mol, stat_Mol, Level, fieldB, fieldE,  nb_mol, params);
    file_out  << t << " ";
    file_out << nb_Mol_in_this_state << " ";
    file_out << niveau_moyen << " ";
    // file_out  << stat_Mol.sigma_pos.mag() << " ";
    // file_out << (stat_Mol.E_pot+stat_Mol.E_cin)/kB/mK/1.5/nb_mol << "  ";
    file_out << stat_Mol.E_cin/kB/mK/1.5/nb_mol << "  ";
    // cout << " N " <<  nb_Mol_in_this_state <<  " T " << stat_Mol.E_cin/kB/mK/1.5/nb_mol << endl;
    liste_Adresses_Mol.clear(); // erase the vector:

    ***/


    // file_out << endl;
    return;
}


// Toutes à la suites en temps
void Sortie_donnee_pop_vJ(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, const int N_two_JX, FitParams &params)
{
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                // file_out << p.name ;
                file_out <<  "  " << p.val_t0 << " " ;
            }
        }
    }


    int **pX=new int*[NX];
    for (int vX=0; vX<NX; vX++)
    {
        pX[vX]=new int[N_two_JX];

    }

    for (int vX = 0; vX < NX; vX++)
        for (int two_JX = 0; two_JX < N_two_JX; two_JX++)
            pX[vX][two_JX] = 0;

    for (int i = 0; i < nb_Mol; i++) // Calcul des populations dans vX,JX
        if (Mol[i].exc == 0)
            pX[Mol[i].v][Mol[i].two_J]++;



    // cout << "    time t =  " << t << endl << endl ;
    // file_out << "    time t =  " << t << endl << endl ;

    file_out << t  << " ";

    for (int vX = 0; vX < NX; vX++)
        for (int two_JX = 0; two_JX < N_two_JX; two_JX++)
        {
            file_out << pX[vX][two_JX] << " ";
        }


    file_out <<  endl ;
    for (int vX = 0; vX < NX; vX++)
    {
        delete[] pX[vX];
    }
    delete[] pX;


    return;

}



// Sortie des populations dans l'état vX à la suite les unes des autres en temps
void Sortie_donnee_pop_v(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, FitParams &params, int number_photons)
{
    // SOrtie des paramètres scannés
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                // file_out << p.name ;
                file_out <<  "  " << p.val_t0 << " " ;
            }
        }
    }

    int *pX=new int[NX];

    for (int vX = 0; vX < NX; vX++)
        pX[vX] = 0;

    for (int i = 0; i < nb_Mol; i++) // Calcul des populations dans vX
        if (Mol[i].exc == 0)
            pX[Mol[i].v]++;

    cout << "    time t =  " << t << endl << endl ;
    file_out << t   << " ";
    file_out << (double) number_photons/nb_Mol << " ";


    for (int vX = 0; vX < NX; vX++)
    {
        //  cout << "  pop[vX="<< vX << "] = " << pX[vX] << endl;
        file_out <<  pX[vX] << " ";
    }

    file_out <<  endl ;

    delete[] pX;

    return;
}




void Sortie_rate_example(ofstream & file_rate, const  vector <double> &rate,  vector <Internal_state> &Level, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int N_Mol, const double t,  FitParams &params)
{
    file_rate<< setprecision(12);
    file_rate << " time t = " << t << endl;

    int nb_rate = rate.size();
    for (int i = 0; i < nb_rate; i++)
    {
        file_rate  << " " << i;
        file_rate  << " " << rate[i];
        int n_mol= reaction_list[i].n_mol;
        int n_laser = reaction_list[i].n_laser;

        file_rate << " " << n_mol ;
        file_rate << " " << n_laser;
        file_rate <<  " " << (reaction_list[i].final_internal_state).two_M ;
        file_rate <<  " " <<  Mol[reaction_list[i].n_mol].two_M << endl;
        file_rate << endl ;
    }
}



void Sortie_rate(ofstream & file_rate, const  vector <double> &rate,  vector <Internal_state> &Level, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int N_Mol, const double t,  FitParams &params)
{
    file_rate<< setprecision(12);
//   file_rate << " time t = " << t << endl;

    const int nb_levels=Level.size(); // Level.size()

 //         1S00     3S1-1    3S11    3S10
 // LEvel   3           4       5       6
    double rate_level_i_vers_level_3[nb_levels] =  {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double rate_level_i_vers_level_4[nb_levels] =  {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double rate_level_i_vers_level_5[nb_levels] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
     double rate_level_i_vers_level_6[nb_levels] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double rate_level_i_total[nb_levels] =  {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double Ecm_i[nb_levels] =  {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    int nb_rate = rate.size();
    double B, E, current_rate ;
    Vecteur3D r,v,k_laser,Bfield,Efield;
    Internal_state Internal_state_in,Internal_state_out;

    for (int i = 0; i < nb_rate; i++)
    {
        current_rate = rate[i];
        int n_mol= reaction_list[i].n_mol;
        int n_laser = reaction_list[i].n_laser;
        //reaction is {n_mol; n_laser; quant_axis; pol_vector; k_eff_laser; final_internal_state;}.
        // For stimlated emission it is { n_mol, n_las, k, k, k, Internal_state_out}
        // For spontaneous emission { num_mol, -1, axe_quant, Vecteur3D(d[0](i,j),d[1](i,j),d[2](i,j)), Vecteur(k,0,0) , Level[j]}

        // file_rate << " " << n_mol ;
        // file_rate << " " << n_laser;

        r = Mol[n_mol].get_pos();
        v= Mol[n_mol].get_vel();
        Bfield= fieldB.get_Field(r);
        B = Bfield.mag();
        Efield= fieldE.get_Field(r);
        B = Bfield.mag();
        E= Efield.mag();
        k_laser = reaction_list[i].k_eff_laser;

        Internal_state_in = Mol[n_mol] ; //  état interne de la molecule
        Internal_state_out = reaction_list[i].final_internal_state ; //  état interne de la molecule après la réaction

//        file_rate  << B;
//
//        file_rate  << " " << i;
//        file_rate  << " " << rate[i];
//
//
//        file_rate <<  "  " << Internal_state_in.deg_number;
//        file_rate << "  " << Internal_state_in.Energy_cm ;
//
//        file_rate <<  "  " << Internal_state_out.deg_number;
//        file_rate << "  " << Internal_state_out.Energy_cm ;


        int i_in = Internal_state_in.deg_number;
        int i_out = Internal_state_out.deg_number;
        Ecm_i[i_in] =  Internal_state_in.Energy_cm;

        if (i_out == 3)   rate_level_i_vers_level_3[i_in] +=  current_rate;
        if (i_out == 4)   rate_level_i_vers_level_4[i_in] +=  current_rate;
        if (i_out == 5)   rate_level_i_vers_level_5[i_in] +=  current_rate;
        if (i_out == 6)   rate_level_i_vers_level_6[i_in] +=  current_rate;
        rate_level_i_total[i_in] += current_rate;


//        file_rate <<  " " << (reaction_list[i].final_internal_state).two_M ;
//        file_rate <<  " " <<  Mol[reaction_list[i].n_mol].two_M << endl;

 //     file_rate << endl ;
    }


    for (int i = 0; i < (int) Level.size() ; i++)
    {

        file_rate <<  " " << B ;
        file_rate <<  " " << E ;
        file_rate <<  " " << v.mag() ;
        file_rate <<  " " << i ;
        file_rate <<  " " << Ecm_i[i] ;
        file_rate <<  " " << Level[i].Energy_cm  ;
         file_rate <<  " " << rate_level_i_vers_level_3[i] ;
        file_rate <<  " " << rate_level_i_vers_level_4[i] ;
        file_rate <<  " " << rate_level_i_vers_level_5[i] ;
         file_rate <<  " " << rate_level_i_vers_level_6[i] ;
        file_rate <<  " " << rate_level_i_total[i] ;
        file_rate << endl ;
    }

}

// Sortie de la variation temporelle de l'intensité  laser
void Sortie_laser_intensity(ofstream & file_out, const vector <Laser> &laser, FitParams &params, int num_laser)
{
    Laser my_laser = laser[num_laser];
    file_out<< setprecision(8);

    for (double t_ps = -100e3; t_ps < 100e3; t_ps++)
    {
        double I_shape = my_laser.intensity_t_nanosecond(t_ps/1000.); // Intensité laser façonnée
        file_out << t_ps/1000. << " " << I_shape << endl;
        // cout << t_ns << " " << I_shape << endl;
    }

    return;
}




void Sortie_laser_spectrum(ofstream & file_out, const vector <Laser> &laser, FitParams &params, int num_laser)
{
    Laser my_laser = laser[num_laser];
    double Energy_transition_laser_cm = cm/my_laser.get_lambda(); // Energie de la transition laser en cm^-1

    file_out<< setprecision(8);
    double intensity0=intensity_Convolution_linewidth(1., 0., 0., my_laser.get_Gamma_Laser(), my_laser.get_type_laser(),num_laser, 0., 0., params); // intensity at resonance

    for (double E_cm = 41100; E_cm < 41200; E_cm = E_cm + 0.01)
    {
        double I_shape = my_laser.transmission_spectrum(E_cm); // Intensité laser façonnée
        double delta =(E_cm - cm/my_laser.get_lambda())*Conv_Ecm_delta  ;// detuning de la transition en s^-1
        double I_laser =  intensity_Convolution_linewidth(1., delta, 0., my_laser.get_Gamma_Laser(), my_laser.get_type_laser(),num_laser,Energy_transition_laser_cm, E_cm, params)/intensity0; //
        file_out << E_cm << " " << I_shape << " " << I_laser << endl;
    }

    return;
}



