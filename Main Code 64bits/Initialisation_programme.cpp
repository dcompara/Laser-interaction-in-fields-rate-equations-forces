
#include "Initialisation_programme.h"

/************************************************************************/
/******************* INITIALISATION DU PROGRAMME ************************/
/************************************************************************/


// Initialisation des niveaux, des transitions, des forces de raies et des durées de vie
// retourne le nombre de niveaux
int   initialisation_trans_mol(const char *nom_file_Levels, const char *nom_file_Lines, vector <Internal_state> &Level, FitParams &params)
{
    int nb_Levels = Level_List(nom_file_Levels, Level, params);
    cout << "Nb de niveaux = " << nb_Levels << endl << endl;

    cout << " Only Zeeman effect no Stark effect " << endl;
    int nb_raies =  Line_List(nom_file_Lines, Level, params);
    cout << "Nb de raies  (en absorption) = " << nb_raies  << endl << endl;

    initialisation_Gamma(Level, params); // Initialisation des durées de vies

    return nb_Levels;
}

//Initialise l'état du générateur de nombre aléatoire
void set_random_generator(gsl_rng  * r, int nombre_seed, const char *nom_file)
{
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;  // "Mersenne Twister" generator by default. Plus rapide que RANLUX. cf. chapitre 17.12 de gsl_ref.
    r = gsl_rng_alloc (T);


//    if (nombre_seed>=0)
//        gsl_rng_set(r,nombre_seed);  //  Pour débugage on initialise le générateur toujours de la même façon
//    else
//    {
//        FILE *fp;
//        fp=fopen("const char *nom_file", "rb");    // Fichier contenant les données du générateur de nombre aléatoire
//        int ran_gen_result = gsl_rng_fread (fp, r); // This function writes the random number state of the random number generator r
//        if (ran_gen_result == GSL_EFAILED)
//            cout << " Problem writing to the file: Data/random_gen.txt" << endl;
//        fclose(fp);
//    }
//    return;
}




// Initialise la position et la vitesse (gaussienne) des molecules
void Init_Molecule(const gsl_rng * r, vector <Molecule> &Mol, const Field fieldB, const Field fieldE, const int Nb_of_type_of_different_Mol, FitParams &params, DataCards data)
{
// Initialize 1st type of Mol from  Mol[0] to Mol[Nom_Mol[0]-1]
// then second type from Mol[Nom_Mol[0]] to Mol[Nom_Mol[0]+Nom_Mol[1]-1] ...
    Mol.clear();
    int nb_of_total_molecules_up_to_now = 0;
    Molecule Mol_i; // Current molecule
    for (int i=0; i < Nb_of_type_of_different_Mol ; i++)
    {
        std::ostringstream oss;
        oss << "[" << i << "]";
        std::string num = oss.str(); // num = "[i]"

        int Nb_Mol = (int)  params.LocateParam(string("N_Mol")+ num)->val;
        string Nom_Mol = data.SParam(string("Nom_Mol")+ num); // nom de la molécule

        int proc[3];   // choix de l'initialisation si potentiel linéaire ou quadratique selon x (0), y(1) et z (2)
        proc[0] = (int) params.LocateParam(string("Procedure_init_x")+ num)->val;
        proc[1] = (int) params.LocateParam(string("Procedure_init_y")+ num)->val;
        proc[2] = (int) params.LocateParam(string("Procedure_init_z")+ num)->val;

        Vecteur3D sigma_pos;
        sigma_pos.set(0,params.LocateParam(string("size_x")+ num)->val);
        sigma_pos.set(1,params.LocateParam(string("size_y")+ num)->val);
        sigma_pos.set(2,params.LocateParam(string("size_z")+ num)->val);

        Vecteur3D const_pos_initial;
        const_pos_initial.set(0,params.LocateParam(string("offset_x")+ num)->val);
        const_pos_initial.set(1,params.LocateParam(string("offset_y")+ num)->val);
        const_pos_initial.set(2,params.LocateParam(string("offset_z")+ num)->val);

        Vecteur3D v0;  // Initial velocity added to the random one
        v0.set(0,params.LocateParam(string("v0_x")+ num)->val);
        v0.set(1,params.LocateParam(string("v0_y")+ num)->val);
        v0.set(2,params.LocateParam(string("v0_z")+ num)->val);

        for (int i=nb_of_total_molecules_up_to_now; i<nb_of_total_molecules_up_to_now+Nb_Mol; i++)
        {

            double masse_mol=MCs; // Default
            double charge = 0.;
// List of molecules (a map could be better!)
            if (Nom_Mol == "H")
            {
                Mol_i.set_name("H");
                masse_mol = MH;
            }
            if (Nom_Mol == "Cs")
            {
                Mol_i.set_name("Cs");
                masse_mol = MCs;
            }
            if (Nom_Mol == "Rb87")
            {
                Mol_i.set_name("Rb87");
                masse_mol = MRb87;
            }
            if (Nom_Mol == "Rb85")
            {
                Mol_i.set_name("Rb85");
                masse_mol = MRb85;
            }
            if (Nom_Mol == "Na")
            {
                Mol_i.set_name("Na");
                masse_mol = MNa;
            }
            if (Nom_Mol == "Cs2")
            {
                Mol_i.set_name("Cs2");
                masse_mol = MCs2;
            }
            if (Nom_Mol == "CO")
            {
                Mol_i.set_name("CO");
                masse_mol = MCO;
            }
            if (Nom_Mol == "NH")
            {
                Mol_i.set_name("NH");
                masse_mol = MNH;
            }
            if (Nom_Mol == "BaF")
            {
                Mol_i.set_name("BaF");
                masse_mol = MBaF;
            }
            if (Nom_Mol == "Li6Cs")
            {
                Mol_i.set_name("Li6Cs");
                masse_mol = MLi6Cs;
            }
            if (Nom_Mol == "Li7Cs")
            {
                Mol_i.set_name("Li7Cs");
                masse_mol = MLi7Cs;
            }
            if (Nom_Mol == "Rb87Cs")
            {
                Mol_i.set_name("Rb87Cs");
                masse_mol = MRb87Cs;
            }
            if (Nom_Mol == "Rb85Cs")
            {
                Mol_i.set_name("Rb85Cs");
                masse_mol = MRb85Cs;
            }
            if (Nom_Mol == "Ps")
            {
                Mol_i.set_name("Ps");
                masse_mol = MPs;
            }
            if (Nom_Mol == "C2minus")
            {
                Mol_i.set_name("C2minus");
                masse_mol = MC2moins;
                charge = QE;
            }
            if (Nom_Mol == "Ps_minus")
            {
                Mol_i.set_name("Ps_minus");
                masse_mol = 3.*ME;
                charge = QE;
            }
            if (Nom_Mol == "P_bar")
            {
                Mol_i.set_name("P_bar");
                masse_mol = Mproton;
                charge = QE;
            }

            Mol_i.set_mass(masse_mol); // Initialise la masse
            Mol_i.set_charge(charge);

            double Temp_ini_x  = params.LocateParam(string("Temp_ini_x")+ num)->val;
            double Temp_ini_y  = params.LocateParam(string("Temp_ini_y")+ num)->val;
            double Temp_ini_z  = params.LocateParam(string("Temp_ini_z")+ num)->val;


            Vecteur3D pos,velocity;

            const double vitesse_therm_x = sqrt(kB*Temp_ini_x/Mol_i.get_mass()); // vitesse thermique en m/s 1/2 m vx^2 = 1/2 kB T
            const double vitesse_therm_y = sqrt(kB*Temp_ini_y/Mol_i.get_mass()); // vitesse thermique en m/s 1/2 m vx^2 = 1/2 kB T
            const double vitesse_therm_z = sqrt(kB*Temp_ini_z/Mol_i.get_mass()); // vitesse thermique en m/s 1/2 m vx^2 = 1/2 kB T


            Vecteur3D Temp_ini(Temp_ini_x,Temp_ini_y,Temp_ini_z);
            Vecteur3D sigmav(vitesse_therm_x,vitesse_therm_y,vitesse_therm_z);
            velocity.gaussian_initialisation (r,sigmav);

            for (int j=0; j<=2; j++)  // 3 axes: x, y ou z
            {
// Exemple Taille de l'échantillion dans le potentiel U(x,y,z) = U(0) + Zeeman_x |x| + Zeeman_y |y| + Zeeman_z z^2
// La distribution est exp(−U(x,y,z)/(k T)) est exp(−|x/sigma_x |)  exp(−|y/sigma_y |)  exp(−z^2/(2 sigma_z^2) )
// const Vecteur3D sigmaSAMPLE(kB*Temp_ini_x/Zeeman_x, kB*Temp_ini_y/Zeeman_y, sqrt(kB*Temp_ini_z/(2.*Zeeman_z)));

                double delta_prime; // dérivée du shift permettant de calculer le moment magnétique effectif
                double sigma_r; // Taille thermique de l'échantillon

                // We addapte the position for different model (some of them also implies velocity modification)
                switch( proc[j] )
                {

                case -1: // ordonne les positions des molécules au départ (selon un axe mettre les autres axes aléatoires)
                {
                    pos.set(j, 3.*(i-Nb_Mol/2.)*sigma_pos.get(j) /Nb_Mol );
                    velocity.set(j, 3.*(i-Nb_Mol/2.)*(sigmav.get(j)) /Nb_Mol );
                    break;
                }

                case 0: // Taille fixe (aléatorie gaussien)
                {
                    pos.set(j, gsl_ran_gaussian (r,sigma_pos.get(j)));
                    break;
                }

                case 1: // Potentiel linéaire --> Laplace. On rappelle que  gsl_ran_laplace (const gsl_rng * r, double a) returns a random variate:  p(x) dx = {1 \over 2 a}  \exp(-|x/a|) dx
                {
                    delta_prime = fabs(delta_field_shift_B_E(Mol, i, fieldB.get_F1().get(j),0.));
                    sigma_r = kB*Temp_ini.get(j)/(HBAR * delta_prime);
                    pos.set(j, gsl_ran_laplace (r,sigma_r));
                    break;
                }

                case 2: // Potentiel quadratique --> gaussien
                {
                    delta_prime = fabs(delta_field_shift_B_E(Mol, i, fieldB.get_F2().get(j),0.)); // dérivée du shift permettant de calculer le moment magnétique effectif
                    sigma_r = sqrt(kB*Temp_ini.get(j)/(2.*HBAR * delta_prime));
                    pos.set(j,gsl_ran_gaussian (r,sigma_r));
                    break;
                }

                case 3: // potential (voltage) elec. quadratique e (electric field linear) but for CHARGED particles (for instance in a Paul trap)--> Gaussien
                {
                    double m_omega2 = - Mol_i.get_charge()* fieldE.get_F1().get(j); //  m omega_j^2  in the potential q V(r_j)=1/2 m omega^2 r_j^2 so  m omega_j^2 = - q F1.j()     (E = - Grad pot)
                    sigma_r = sqrt(kB*Temp_ini.get(j)/m_omega2); // because 1/2 m omega^2 sigma_r^2 = 1/2 k_B T --> sigma_r^2=kB T/(m omega^2))
                    // sqrt(3.*kB) is able to give the proper E_pot
                    pos.set(j,gsl_ran_gaussian (r,sigma_r));
                    break;
                }

                case 4: // random gaussian in r and pure gaussian (non random) in v
                {
                    pos.set(j,gsl_ran_gaussian (r,sigma_pos.get(j)));
                    double vj =  gsl_cdf_gaussian_Pinv((1+i)/(Nb_Mol+1.),sigmav.get(j) ); // x^2 = m v^2/(kB T) ; x = v/sigma_v so gaussian distribution is 1/(sqrt(2 pi)) *  exp (-x^2/2)
                    velocity.set(j,  vj); // We cut in N+1 interval with equiproba to have a molecule.
                    // So, we put a molecule exactly at x_(i+1) (i=0,...,N-1) where int_(x_(i)^x_(i+1))  e^(-x^2/2) dx =1/(N+1)
                    break;
                }

                case 5: // 5 effusive beam. Meaning as in case 0 but we keep only the positive velocities
                {
                    double vj;
                    do
                    {
                        vj =  gsl_ran_gaussian (r,sigmav.get(j));
                    }
                    while(vj<0);
                    velocity.set(j,  vj);
                    pos.set(j, gsl_ran_gaussian (r,sigma_pos.get(j)));
                    break;
                }

                }
            }

            velocity = velocity + v0;
            pos = pos + const_pos_initial;  // Initial velocity and position added to the random one


            Mol_i.set_pos(pos);      // Initialise la position (gaussienne) des molecules
            Mol_i.set_vel(velocity);   // Initialise la vitesse (gaussienne) des molecules
            Mol.push_back(Mol_i); // Mol[i]=Mol_i;
        }
        nb_of_total_molecules_up_to_now += Nb_Mol;
    }

    return;
}


// Initialise les lasers
void Init_Laser(vector <Laser> &laser, const int Nb_lasers, FitParams &params, const char *nom_file_spectrum, const char *nom_file_intensity)
{
    laser.clear();
    for (int i=0; i < Nb_lasers ; i++)
    {
        std::ostringstream oss;
        oss << "[" << i << "]";
        std::string num = oss.str(); // num = "[i]"

        Vecteur3D waist_position;
        waist_position.set(0,params.LocateParam(string("waist_pos_x") + num)->val);
        waist_position.set(1,params.LocateParam(string("waist_pos_y") + num)->val);
        waist_position.set(2,params.LocateParam(string("waist_pos_z") + num)->val);

        Vecteur3D Direction_laser;
        Direction_laser.set(0,params.LocateParam(string("direction_x") + num)->val);
        Direction_laser.set(1,params.LocateParam(string("direction_y") + num)->val);
        Direction_laser.set(2,params.LocateParam(string("direction_z") + num)->val);

        Vecteur3D waist_laser;
        waist_laser.set(0,params.LocateParam(string("waist") + num)->val);
        waist_laser.set(1,params.LocateParam(string("waist") + num)->val);

        double Gamma_L_MHz  = params.LocateParam(string("Gamma_L_MHz") + num)->val * params.LocateParam("scale_Gamma")->val;
        double Power = params.LocateParam(string("Power") + num)->val * params.LocateParam("scale_Power")->val; // Parametre multiplicatif de la puissance des lasers
        double Energie_cm = params.LocateParam(string("Energie_cm") + num)->val + params.LocateParam("Offset_Detuning_cm")->val; // Parametre multiplicatif de la puissance des lasers

        Vecteur3D polarisation_laser; // Polarisation -1,0,+1: sigma-,pi,sigma+
        polarisation_laser.set(0,params.LocateParam(string("Pol_circulaire_right_sm") + num)->val);
        polarisation_laser.set(2,params.LocateParam(string("Pol_circulaire_left_sp") + num)->val);

        double polar_angle_degree = params.LocateParam(string("polar_angle_degree") + num)->val;

        int type_laser = (int) params.LocateParam(string("type_laser") + num)->val;

        int coherent_avec_laser_num = (int) params.LocateParam(string("coherent_avec_laser_num") + num)->val;

        Laser current_Laser;
        current_Laser.set_waist_pos(waist_position);
        current_Laser.set_direction(Direction_laser);
        current_Laser.set_waist(waist_laser);
        current_Laser.set_energy_cm(Energie_cm);
        current_Laser.set_Gamma_Laser_MHz(Gamma_L_MHz);
        current_Laser.set_Power(Power);
        current_Laser.set_polarisation(polarisation_laser);
        current_Laser.set_polar_angle(polar_angle_degree);
        current_Laser.set_type_laser(type_laser);
        current_Laser.set_coherent_avec_laser_num(coherent_avec_laser_num);
        string nom_file_spectrum_short = (string(nom_file_spectrum)).substr(0,string(nom_file_spectrum).length() -4); // enlève le .dat
        string nom_file_spectrum_long = nom_file_spectrum_short + num + ".dat";
        const char *nom_file_spectrum = (nom_file_spectrum_long).c_str(); // Nom_file[numero_laser].dat
        current_Laser.read_Spectrum(nom_file_spectrum); // Lit le fichier du spectre laser (rien si le ficheir n'existe pas)
        string nom_file_intensity_short = (string(nom_file_intensity)).substr(0,string(nom_file_intensity).length() -4); // enlève le .dat
        string nom_file_intensity_long = nom_file_intensity_short + num + ".dat";
        const char *nom_file_intensity = (nom_file_intensity_long).c_str(); // Nom_file[numero_laser].dat
        current_Laser.read_Intensity(nom_file_intensity); // Lit le fichier du spectre laser (rien si le ficheir n'existe pas)
        laser.push_back(current_Laser); // laser[i] = current_Laser;
    }

    return;
}


// Modification temporelle des paramètres
void Modif_Param(const double t,  FitParams &params)
{
    for (ParamIterator i=params.begin(); i != params.end(); ++i) // boucle sur les paramètres
    {
        Param &p = **i;
        // cout << " t " << t << "  " << p << endl << endl;
        if( p.is_time_dependent == false ) // Pas de boucle pour ce paramètre
            continue;
        else
            p.val = p.val_t0*exp(-t/p.tau); // On modifie le paramètre
    }
    return;
}



/************************************************************************/
/******************* INITIALISATION DES CHAMPS ************************/
/************************************************************************/

// Initialise les champs,  zero par défaut.
void Init_Field(Field &fieldB, Field &fieldE, FitParams &params, const char *nom_Magn_file, const char *nom_Elec_file)
{

    /******** Magnetic Field *************/

    double B_x,B_y,B_z;
    B_x  = params.LocateOrAddNonConstParam("B_x")->val;
    B_y  = params.LocateOrAddNonConstParam("B_y")->val;
    B_z  = params.LocateOrAddNonConstParam("B_z")->val;
    Vecteur3D B_Field =  Vecteur3D(B_x, B_y, B_z);

    double grad_B_x,grad_B_y,grad_B_z;
    grad_B_x  = params.LocateOrAddNonConstParam("grad_B_x")->val;
    grad_B_y  = params.LocateOrAddNonConstParam("grad_B_y")->val;
    grad_B_z  = params.LocateOrAddNonConstParam("grad_B_z")->val;
    Vecteur3D Grad_B =  Vecteur3D(grad_B_x, grad_B_y, grad_B_z);

    double grad_grad_B_x, grad_grad_B_y, grad_grad_B_z;
    grad_grad_B_x  = params.LocateOrAddNonConstParam("grad_grad_B_x")->val;
    grad_grad_B_y  = params.LocateOrAddNonConstParam("grad_grad_B_y")->val;
    grad_grad_B_z  = params.LocateOrAddNonConstParam("grad_grad_B_z")->val;
    Vecteur3D grad_grad_B =  Vecteur3D(grad_grad_B_x, grad_grad_B_y, grad_grad_B_z);

    int n_value_B = (int) params.LocateOrAddNonConstParam("n_value_B")->val;
    double Bn_x, Bn_y, Bn_z,r_cut;
    Bn_x  = params.LocateOrAddNonConstParam("Bn_x")->val;
    Bn_y  = params.LocateOrAddNonConstParam("Bn_y")->val;
    Bn_z  = params.LocateOrAddNonConstParam("Bn_z")->val;
    r_cut = params.LocateOrAddNonConstParam("r_cut")->val;
    Vecteur3D Bn =  Vecteur3D(Bn_x, Bn_y, Bn_z);

    fieldB.set_F0(B_Field);
    fieldB.set_F1(Grad_B);
    fieldB.set_F2(grad_grad_B);
    fieldB.set_n(n_value_B);
    fieldB.set_Fn(Bn);

//    int type_field = params.LocateOrAddNonConstParam("type_field")->val;

    fieldB.type_field = Magnetic_Field;
    fieldE.type_field = Electric_Field;

    int type_field_read_B = params.LocateOrAddNonConstParam("type_field_read_B")->val;

    int Nb_bobines = params.LocateOrAddNonConstParam("Nb_bobines")->val;
    double gap_bobines = params.LocateOrAddNonConstParam("gap_bobines")->val;
    double courant_bobines = params.LocateOrAddNonConstParam("courant_bobines")->val;
    double rayon_bobines = params.LocateOrAddNonConstParam("rayon_bobines")->val;
    int is_Helmholtz = params.LocateOrAddNonConstParam("is_Helmholtz")->val;

    int type_of_field_for_internal_state_shift  = params.LocateOrAddNonConstParam("type_of_field_for_internal_state_shift")->val;

    fieldB.Nb_bobines = Nb_bobines;
    fieldB.gap_bobines = gap_bobines;
    fieldB.courant_bobines = courant_bobines;
    fieldB.rayon_bobines = rayon_bobines;
    fieldB.is_Helmholtz = is_Helmholtz;
    fieldB.type_field_read = type_field_read_B;


    if (fieldB.type_field_read >= 2) // We read the files
    {
        if(nom_Magn_file !=NULL) // At the first initialization we give the name of the file and we have to read the file and initialize the INterpolator
        {
            if (fieldB.Read_Matrix(nom_Magn_file)==0)
            {
                cerr << "Empty Magnetic field file " << nom_Magn_file << endl;
                exit(1);
            }
            cout << "We read Magnetic field file " ;
            cout << "Nb lignes " << fieldB.Read_Matrix(nom_Magn_file) << endl; // If no file this does not do anything and so do not modify the matrix_field
            fieldB.Init_Interpolator(); // Give the size, offset and spacing of the grid used
        }
    }




    /******** Electric Field *************/

    double E_x,E_y,E_z;
    E_x  = params.LocateOrAddNonConstParam("E_x")->val;
    E_y  = params.LocateOrAddNonConstParam("E_y")->val;
    E_z  = params.LocateOrAddNonConstParam("E_z")->val;
    Vecteur3D E_Field =  Vecteur3D(E_x, E_y, E_z);

    double grad_E_x,grad_E_y,grad_E_z;
    grad_E_x  = params.LocateOrAddNonConstParam("grad_E_x")->val;
    grad_E_y  = params.LocateOrAddNonConstParam("grad_E_y")->val;
    grad_E_z  = params.LocateOrAddNonConstParam("grad_E_z")->val;
    Vecteur3D Grad_E =  Vecteur3D(grad_E_x, grad_E_y, grad_E_z);

    double grad_grad_E_x, grad_grad_E_y, grad_grad_E_z;
    grad_grad_E_x  = params.LocateOrAddNonConstParam("grad_grad_E_x")->val;
    grad_grad_E_y  = params.LocateOrAddNonConstParam("grad_grad_E_y")->val;
    grad_grad_E_z  = params.LocateOrAddNonConstParam("grad_grad_E_z")->val;
    Vecteur3D grad_grad_E =  Vecteur3D(grad_grad_E_x, grad_grad_E_y, grad_grad_E_z);

    fieldE.set_F0(E_Field);
    fieldE.set_F1(Grad_E);
    fieldE.set_F2(grad_grad_E);

    int n_value_E = (int) params.LocateOrAddNonConstParam("n_value_E")->val;
    double En_x, En_y, En_z;
    En_x  = params.LocateOrAddNonConstParam("En_x")->val;
    En_y  = params.LocateOrAddNonConstParam("En_y")->val;
    En_z  = params.LocateOrAddNonConstParam("En_z")->val;
    Vecteur3D En =  Vecteur3D(En_x, En_y, En_z);

    fieldE.set_n(n_value_E);
    fieldE.set_Fn(En);

    int type_field_read_E = params.LocateOrAddNonConstParam("type_field_read_E")->val;
    fieldE.type_field_read = type_field_read_E;


    if (fieldE.type_field_read >= 2) // We read the files
    {
        if(nom_Elec_file == NULL) // At the first initialization we give the name of the file and we have to read the file and initialize the INterpolator
        {
            cout << " because no Elec file name !" << endl;
        }

        if(nom_Elec_file != NULL) // At the first initialization we give the name of the file and we have to read the file and initialize the INterpolator
        {
            if (fieldE.Read_Matrix(nom_Elec_file)==0)
            {
                cerr << "Empty Electric field file " << nom_Elec_file << endl;
                exit(1);
            }
            fieldE.Read_Matrix(nom_Elec_file);
            fieldE.Init_Interpolator(); // Give the size, offset and spacing of the grid used
        }
    }

}

