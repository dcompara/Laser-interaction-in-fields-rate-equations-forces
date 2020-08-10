/*
 Fichier de sortie des donn�es mol�culaires: transition, probablilit�es ...
*/

#include "sortie_donnees.h"


//Sauve l'�tat du g�n�rateur de nombre al�atoire pour ne jamais reprendre le m�me
void save_random_generator(gsl_rng * r, const char *nom_file)
{
    FILE *fp;
    fp=fopen(nom_file, "wb");    // Fichier contenant les donn�es du g�n�rateur de nombre al�atoire
    int ran_gen_result = gsl_rng_fwrite (fp, r); // This function writes the random number state of the random number generator r
    if (ran_gen_result == GSL_EFAILED)
        cout << " Problem writing to the file: Data/random_gen.txt" << endl;
    fclose(fp);
}


/* Cr�ation du fichier de sortie des donn�es totales  */
void sortie_fichier(ofstream &  /* file2 */, Molecule /* Mol[] */)
{

    /* A mettre si on veut une num�rotation automatique (ainsi que la partie finale)


    ifstream file_in("Data/file_num.txt"); // Ouvre le fichier contenant le futur num�ro de fichier
    int num;
    file_in >> num  ;
    file_in.close();

    char snum[256];
    sprintf(snum,"%d",num);
    string ssnum = snum ;
    string nom_fichier = "Data/param_Ryd_" + ssnum + ".txt" ;

    ofstream file_out(nom_fichier.c_str()); // cr�er le fichier de donn� de sortie param_Ryd_N� .txt contenant i, exc[i],  x[i], y[i], z[i], pot[i], pot_ji (plus proche voisin), theta_ji, d_ji

    */


    /* A mettre si on veut une num�rotation automatique (ainsi que la partie comment�e initiale )
        file_out.close();
        ofstream file_num_out("Data/file_num.txt"); // cr�er le fichier contenant le futur num�ro de fichier
        file_num_out << ++num;
        file_num_out.close();
    */
}



void Sortie_donnee_etat_int_simple(ofstream & file_out, const vector <Molecule> &Mol, const vector <Laser> &my_laser, const double t, FitParams &params)
{
    file_out<< setprecision(8);
    for(int i=0; i<(int) Mol.size(); i++)
    {
        file_out << "   " << t << " " << i << " " << Mol[i].exc  << "   " << Mol[i].two_J/2.  << "   " << Mol[i].two_N/2. << "   " << Mol[i].two_M/2. << "   " << Mol[i].Energy0_cm << endl ;
    }
    // cout << endl << "    time t =  " << t << endl  ;
    file_out << endl;

    return;
}





void Sortie_test_debug(ofstream & file_out,  vector <Molecule> &Mol,  const vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params,  DataCards &data, const int number_photons)
{

//    int number_matrix_lines = 2480;
//    fieldB.Calculate_Derivative_Matrix(number_matrix_lines, "Data/Na/MagneticField2D_derivative.dat");

    double step_r = 0.001;
    double step_z = 0.001;
    double nb_steps = 100;

    // for (double  x= -nb_steps*step_r; x < nb_steps*step_r; x+=step_r)
    //for (double  y= -nb_steps*step_r; y < nb_steps*step_r; y+=step_r)
    for (double  x= 0.003; x < 0.005; x+=step_r)

        for (double  z= -nb_steps*step_z; z < nb_steps*step_z; z+=step_z)
        {
            double y=0.;
            file_out << sqrt(x*x+y*y) << " " <<  x << " " << y << " " << z << " "  << fieldB.get_Field(Vecteur3D(x,y,z)) << " " << fieldB.get_grad_field_F2(Vecteur3D(x,y,z)) << endl;
        }
    file_out << endl;

    return;
}





void Sortie_donnee(ofstream & file_out,  vector <Molecule> &Mol,  vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int nb_mol,  FitParams &params,  DataCards &data, const int number_photons)
{
    //set_pot_all_mol(Mol, fieldB, fieldE, laser, t, nb_mol, params); //Met � jour de tous les potentiels (gravit�, PAS dipolaire, magn�tique, electrique, ...) avec la nouvelle position pour uen sortie
// ATTENTIION THIS DOES NOT WORK FOR THE POTENTIALS


// SOrtie des param�tres scann�s
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les param�tres
        {
            Param &p = **i;
            if (p.is_scanned == true)
            {
                // file_out << p.name ;
                // file_out << p.val_t0 << " " ;
            }
        }
    }

    file_out<< setprecision(8);

    const int i = (int) params.LocateParam("num_niveau_etudie")->val; // num�ro du niveau �tudi� pour faire des stats. -1 pour toutes les mol�cules

    Stat stat_Mol;
    stat_molecule_un_niveau(Mol, stat_Mol, Level,  fieldB, fieldE, i, Mol.size(), params);


    for (int i = 0; i < nb_mol; i++)
    {
        set_pot_mol(Mol, i, fieldB, fieldE, laser, t, params); //Met � jour de tous les potentiels (gravit�, PAS dipolaire, magn�tique, electrique, ...) avec la nouvelle position pour une sortie
    }

//    double Temp_ini_z  = params.LocateParam("Temp_ini_z[0]")->val;
//    const double vitesse_therm_z = sqrt(kB*Temp_ini_z/Mol[0].get_mass()); // vitesse thermique en m/s 1/2 m vx^2 = 1/2 kB T



//    file_out  << t << " ";
//    file_out << number_photons << " ";
//    file_out << stat_Mol.Temp_3D_50 << " ";
//    file_out  << stat_Mol.population << " ";
//    file_out  << stat_Mol.Temp_1D_50.z() << " ";


//  Attention relative � la temp�rature E_pot = 3/2 k T, E_cin =3/2 kT; E tot=3 kT

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


        MatrixXd d[3] ;
        SelfAdjointEigenSolver<MatrixXd> es; // eigenstates and eigenvalues
        Vecteur3D v;
        v = Mol[i].get_vel();

        Vecteur3D Bfield;
        Bfield= fieldB.get_Field(Mol[i].get_pos());
        double B = Bfield.mag();
        double v_perp= (v.cross(Bfield)).mag()/B;
        Diagonalization_Energy_dipole(Level, B,  v_perp, es, d);
        for (int j=0; j< Level.size(); j++) //  we scan over the levels to calculate the parameter
        {
            double tripletness = 0.; //This is the parameter we want to calculate (here the triplet character)

            for (int j0=0; j0< Level.size(); j0++) //  | j> =  sum_|j>_O   0_<j | j>  |j>_0  so we scan over |j>_0 hat is the order in the Level file
                // 0_<j | j>  is given by   es.eigenvectors()(j0,j) . This is the coefficient of the |j> level (ordered in Energy) on the |j>_0 Level (the order in the Level file). We round it to 100%
            {
                // file_out << B << " " << v_perp << " " << i << " " << j  << "  " << j0 << " " << Level[j0].two_M << " " << abs(round(100.*es.eigenvectors()(j0,j)))/100 << endl;
                if (Level[j0].two_Omega == 2) // If the state is triplet (2S+1=3 so S=1) we look on the decomposition, |0_<i | i>|^2 , and sum them
                {
                    tripletness +=  pow((es.eigenvectors()(j0,j)),2); // sum_|triple, j>_O   |0_<j | j>|^2
                }
            }
            // file_out << endl;
            // PARAMETER THAT GIVE THE TRIPLETNESS OF THE STATE //
            file_out << B << " " << v_perp << " " << j << " " << Level[j].Energy_cm << " " << tripletness << endl;
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


void Sortie_donnee_electrons(ofstream & file_out,  vector <Molecule> &Mol,  const vector <Internal_state> &Level, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t, const int number_mol,  FitParams &params,  DataCards &data, const int number_photons)
{
    Vecteur3D r,v;
    r = Mol[number_mol].get_pos();
    v = Mol[number_mol].get_vel();
    Vecteur3D E = fieldE.get_Field(r);
    double V = fieldE.get_electric_potential(r);
    double Vcoulomb = get_coulomb_potential(Mol, r, Mol.size());

    file_out<< setprecision(8);
    file_out  << t << " ";
    file_out << r << "  ";
    file_out << E.mag() << "  ";
    file_out <<  get_coulomb_field(Mol, number_mol, r).mag() << " ";
    file_out << V << "  ";
    file_out << Vcoulomb << "  ";

    file_out << endl;
    return;
}


// Toutes � la suites en temps
void Sortie_donnee_pop_vJ(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, const int N_two_JX, FitParams &params)
{
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les param�tres
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



// Sortie des populations dans l'�tat vX � la suite les unes des autres en temps
void Sortie_donnee_pop_v(ofstream & file_out, const vector <Molecule> &Mol, const int nb_Mol, const double t, const int NX, FitParams &params, int number_photons)
{
    // SOrtie des param�tres scann�s
    if( ((int) params.LocateParam("is_param_scanned_out")->val) == ((int) true) )
    {
        for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les param�tres
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


void Sortie_rate(ofstream & file_rate, const  vector <double> &rate,  vector <Internal_state> &Level, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int N_Mol, const double t,  FitParams &params)
{
    file_rate<< setprecision(12);
//   file_rate << " time t = " << t << endl;

    int nb_levels=Level.size(); // Level.size()

    double rate_level_i_vers_level_2[nb_levels] =  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double rate_level_i_vers_level_3[nb_levels] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double rate_level_i_total[nb_levels] =  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double Ecm_i[nb_levels] =  {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    int nb_rate = rate.size();
    double B, current_rate ;
    Vecteur3D r,v,k_laser,Bfield;
    Internal_state Internal_state_in,Internal_state_out;

    for (int i = 0; i < nb_rate; i++)
    {
        // file_rate  << " " << i;
        // file_rate  << " " << rate[i];
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
        k_laser = reaction_list[i].k_eff_laser;
//        file_rate  << " " << r ;
//        file_rate  << " " << v.mag(); ;
//        file_rate  << " " << B ;
//        file_rate  << " " << k_laser.mag()/(2*pi*100.) ;  // k = 2*pi*100.*Energy_transition_cm;
//        double E = fieldE.get_Field(r).mag();
        Internal_state_in = Mol[n_mol] ; //  �tat interne de la molecule
//        double Energy_in = Internal_state_in.Energy_cm;
//        double Energy_in = Internal_state_in.Energy0_cm + (Internal_state_in.Energy_Shift_B_cm(B) + Internal_state_in.Energy_Shift_E_cm(E));
//        file_rate  << " " <<  Internal_state_in.deg_number ;

        Internal_state_out = reaction_list[i].final_internal_state ; //  �tat interne de la molecule apr�s la r�action
//        double Energy_out = Internal_state_out.Energy_cm;
//
//
//  double Energy_out = Internal_state_out.Energy0_cm + (Internal_state_out.Energy_Shift_B_cm(B) + Internal_state_out.Energy_Shift_E_cm(E));
//        double Energy_transition_laser_cm = cm/laser[n_laser].get_lambda(); // Energie de la transition laser en cm^-1


//        file_rate <<  " " << Bfield.z();
//        file_rate <<  " " << Energy_out- Energy_in - Energy_transition_laser_cm;
//        file_rate << "  " << Energy_in ;

//        file_rate <<  "  " << Internal_state_in.deg_number;
//        file_rate <<  "  " << Internal_state_out.deg_number;
//        file_rate << "  " << Internal_state_in.Energy_cm ;


        int i_in = Internal_state_in.deg_number;
        int i_out = Internal_state_out.deg_number;
        Ecm_i[i_in] =  Internal_state_in.Energy_cm;

        if (i_out == 2)   rate_level_i_vers_level_2[i_in] +=  current_rate;
        if (i_out == 3)   rate_level_i_vers_level_3[i_in] +=  current_rate;
        rate_level_i_total[i_in] += current_rate;


        /*****   CALCUL of parameters for the dipoles (in Debye)or diagonalization,
         because the Internal States are not the correct one we need to diagonalized in order to find
         the proper one


         for all 21 levels              j        Energy_j         Energy_J-Energy_out     TRIPLET_CHARACTER

            *******/

        /****
                if ( (params.LocateParam("is_Levels_Lines_Diagonalized")->val) )
                {
                    int level_in = Internal_state_in.deg_number;
                    int level_out = Internal_state_out.deg_number;

                    file_rate  << " " << level_in ; // is the number for the level of the molecule
                    file_rate  << " " << level_out ; // is the number for the level of the molecule

                    file_rate  << endl;
                    MatrixXd d[3] ;
                    SelfAdjointEigenSolver<MatrixXd> es; // eigenstates and eigenvalues

                    double v_perp= (v.cross(Bfield)).mag()/B;
                    Diagonalization_Energy_dipole(Level, B, v_perp, es, d);


        //            for (int i=0; i<Level.size(); i++)
        //                file_rate << " i  "  << i  << " E_i " <<  round(100.*(es.eigenvalues()(i)))/100. << endl;
        //
        //
        //
        //            for (int i=0; i<Level.size(); i++)
        //                for (int j=0; j<Level.size(); j++)
        //                {
        //                    file_rate << " i  "  << i  << " j " << j  << " evec_{ij} " <<  round(100.*(es.eigenvectors())(i,j))/100. << endl;
        //                }



                    for (int j=0; j< Level.size(); j++) //  we scan over the levels to calculate the parameter
                    {
                        double param = 0.; //This is the parameter we want to calculate (here the triplet character)

                        for (int j0=0; j0< Level.size(); j0++) //  | j> =  sum_|j>_O   0_<j | j>  |j>_0  so we scan over |j>_0 hat is the order in the Level file
                            // 0_<j | j>  is given by   es.eigenvectors()(j0,j) . This is the coefficient of the |j> level (ordered in Energy) on the |j>_0 Level (the on in the Level file)
                        {
                            if (Level[j0].two_Omega == 2) // If the state is triplet (2S+1=3 so S=1) we look on the decomposition, |0_<i | i>|^2 , and sum them
                            {
                                param +=  pow((es.eigenvectors()(j0,j)),2); // sum_|triple, j>_O   |0_<j | j>|^2
                            }
                        }
                        //                       file_rate << " " << j << " " << es.eigenvalues()(j)<< " " << es.eigenvalues()(j) - es.eigenvalues()(level_out) << " " << abs(round(10.*param))/10. << endl;
                        //          file_rate << " " << es.eigenvalues()(j) << " " << abs(round(100.*param))/100. ;
                    }
                }
        ****/

//        file_rate <<  " " << (reaction_list[i].final_internal_state).two_M ;
//        file_rate <<  " " <<  Mol[reaction_list[i].n_mol].two_M << endl;

//        file_rate << endl ;
    }





    MatrixXd d[3] ;
    SelfAdjointEigenSolver<MatrixXd> es; // eigenstates and eigenvalues
    Diagonalization_Energy_dipole(Level, B,  v.mag(), es, d);
    for (int i = 0; i <Level.size() ; i++)
    {

        file_rate <<  " " << B ;
        file_rate <<  " " << v.mag() ;
        file_rate <<  " " << i ;
        file_rate <<  " " << Ecm_i[i] ;
        file_rate <<  " " << Level[i].Energy_cm  ;
        file_rate <<  " " << rate_level_i_vers_level_2[i] ;
        file_rate <<  " " << rate_level_i_vers_level_3[i] ;
        file_rate <<  " " << rate_level_i_total[i] ;
        file_rate << endl ;
    }

}



void Sortie_laser_spectrum(ofstream & file_out, const vector <Laser> &laser, FitParams &params, int num_laser)
{
    Laser my_laser = laser[num_laser];
    double Energy_transition_laser_cm = cm/my_laser.get_lambda(); // Energie de la transition laser en cm^-1

    file_out<< setprecision(8);
    double intensity0=intensity_Convolution_linewidth(1., 0., 0., my_laser.get_Gamma_Laser(), my_laser.get_type_laser(),num_laser, 0., 0., params);

    for (int E_cm = 0; E_cm < 20000; E_cm++)
    {
        double I_shape = my_laser.transmission_spectrum(E_cm); // Intensit� laser fa�onn�e
        double delta =(E_cm - 0.01/my_laser.get_lambda())*Conv_Ecm_delta  ;// detuning de la transition en s^-1
        double I_laser =  intensity_Convolution_linewidth(1., delta, 0., my_laser.get_Gamma_Laser(), my_laser.get_type_laser(),num_laser,Energy_transition_laser_cm, E_cm, params)/intensity0; //
        file_out << E_cm << " " << I_shape << " " << I_laser << endl;
    }

    return;
}

// Debug. Gives state, potential, ...
void Sortie_debug(ofstream & file_rate, const  vector <double> &rate, const vector <type_codage_react> &reaction_list, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const int n_reac, const double t,  FitParams &params)
{
    // t	rate	nlas	z	B	delta	E_fin	E_in	2J_fin	2J_in	2Mfin	2M_in

    file_rate<< setprecision(12);
    //  file_rate << " time t = " ;
    file_rate << t << " ";

    int i=  n_reac;

//   file_rate  <<  " rate[" << i << "] = ";
    file_rate << rate[i] << " ";
    int n_mol= reaction_list[i].n_mol;
    int n_laser = reaction_list[i].n_laser;
//   file_rate << " nMol = " ;
    //  file_rate<< n_mol << " ";
//   file_rate << " nlas = ";
    file_rate << n_laser << " ";

    Vecteur3D r;
    r = Mol[n_mol].get_pos();
    //  file_rate << " pos = " ;
    file_rate << r.z() << " ";

    double B = fieldB.get_Field(r).mag();
    double E = fieldE.get_Field(r).mag();
    Internal_state Internal_state_in = Mol[n_mol] ; //  �tat interne de la molecule
    double Energy_in = Internal_state_in.Energy0_cm + (Internal_state_in.Energy_Shift_B_cm(B) + Internal_state_in.Energy_Shift_E_cm(E));
    Internal_state Internal_state_out = reaction_list[i].final_internal_state ; //  �tat interne de la molecule apr�s la r�action
    double Energy_out = Internal_state_out.Energy0_cm + (Internal_state_out.Energy_Shift_B_cm(B) + Internal_state_out.Energy_Shift_E_cm(E));
    double Energy_transition_laser_cm = cm/laser[n_laser].get_lambda(); // Energie de la transition laser en cm^-1


//    file_rate << " B " ;
    file_rate << B << " ";
//   file_rate << " detuning_cm " ;
    file_rate << abs(Energy_out- Energy_in) - Energy_transition_laser_cm << " ";
//   file_rate << " Efin0 " ;
//   file_rate << Internal_state_out.Energy0_cm<< " ";
    //  file_rate << "E_in0 ";
    //  file_rate << Internal_state_in.Energy0_cm << " ";
    //  file_rate << " Efin " ;
    file_rate << Energy_out << " ";
    //  file_rate  << "E_in ";
    file_rate << Energy_in << " " ;
    //  file_rate << " 2J_fin ";
    file_rate << (reaction_list[i].final_internal_state).two_J << " ";
    //  file_rate << " 2J_in ";
    file_rate << Mol[reaction_list[i].n_mol].two_J << " ";
    file_rate << Mol[reaction_list[i].n_mol].two_M << " " ;

    file_rate << endl ;

}


// collisional cross section for charge exchange Ps Pbar.

double Cross_section_Ps(double v,  const int n)
{
    double s1= 1.32e-16;
    double s2= 1.12e-15;
    double v_electron = ALPHA * C/(2. * n);
    double kv = v/v_electron;
    double sigma = n*n*n*n*(s1/(kv * kv) + s2)/(1.+pow(kv/1.8,20)); //  Cf PRA 94, 022714 (2016) fitted

    double B = 4.5*tesla;  // MAGNETIC FIELD (PUT HERE BY HAND !!!!!)

    double E_max_field_ionization = (-QE/(4.*pi*EPSILON0 *(2.*A0)*(2.*A0)))/(16.*n*n*n*n) ; // Maximum field for field iniastion. 1/16n^4 can be replaced by 1/9 n^4  in pure electric field

    if ( v*B < E_max_field_ionization)
    {
        return sigma ;
    }
    else
    {
        return 0. ; // because the particles is field ionized
    }
}



