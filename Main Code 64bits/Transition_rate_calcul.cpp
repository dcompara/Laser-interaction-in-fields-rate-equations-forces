/*
  Name: Algorithme de calcul des tauxf
  Author: Daniel Comparat
  Date: 16/12/08


REMARQUES:
On pourrait traiter la (pré ou photo)-dissociation comme un parametre phénoménologue ainsi que les transition a deux photons aussi
la desexcitation dans le continuuum peut être prise en compte avec Sum_FC != 1 la partie non 1 est le continuum


ATTENTION
Unités des fichier (cm^-1)
*/



#include "Transition_rate_calcul.h"
#include "diagonalization.h"

#include <sstream> //


/************************************************************************/
/************************** Interaction laser - 2 niveaux *************/
/************************************************************************/



// Taux d'emission spontanée d'une transition (sans prendre en compte l'ordre des états:
// Voir aussi  Internal_state::Einstein_desex_rate() qui lui est zéro (VERY_SMALL_NUMBER en fait) si l'ordre est initial en bas!
double Gamma_spon(const double dip_Debye, const double energie_cm)
{
    return dip_Debye*dip_Debye*Spol_Debye_A_s*energie_cm*energie_cm*energie_cm;
    //  Spol_Debye_A_s = (8e6*pi*pi*C*C*C*Debye*Debye)/(3.*EPSILON0*C*C*C*HBAR); // 3.13618932*10^-7 = Conversion HonlLondon (en Debye) en A Einstein (s^-1) si energie en cm^-1
}


// Forme  gaussienne Gamma = FWHM
inline double Gauss(const double I_tot, const double delta, const double Gamma)
{
    double sigma = Gamma/(2.* sqrt(2.*log(2))); // Gamma = FWHM  = 2 \sqrt{2 \ln(2)} \sigma  = 2,3548 \sigma
    return (I_tot/(sigma * sqrt(2.*pi))) * exp(-delta*delta/(2.*sigma*sigma));
}

// Forme laser plate
// retourne l'intensité (spectrale) au décalage delta souhaité
// Gamma = FWHM
inline double intensity_flat(const double I_tot, const double delta, const double Gamma)
{
    if (abs(delta) < Gamma/2.)
        return I_tot/Gamma; // Gamma = FWHM
    else
        return 0.;
}


// Forme laser lorentzienne Gamma = FWHM
// I(omega) = I*Gamma_L/(2.*pi)/(Gamma_L^2/4 + delta^2)
inline double Lorentz(const double I_tot, const double delta, const double Gamma)
{
    return I_tot*(Gamma/(2.*pi))/((Gamma/2.)*(Gamma/2.) + delta*delta);
}



/**

There is several possible approximation for the Voigt profile (Lorentzian covlution with Gaussian) that we could use when dealing with gaussian laser shape.
CF RAPID COMPUTATION OF THE VOIGT PROFILE S. R. DRAYNI   or J.L. Kielkopf, J. Opt. Soc. Am. 63 , 987 (1973) or "Empirical fits to the Voigt line width: A brief review"
"Simple empirical analytical approximation to the Voigt profile" or "the Voigt Profile as a Sum of a Gaussian and a Lorentzian Functions"
  or http://arxiv.org/pdf/0805.2274.pdf or "Determination of the Gaussian and Lorentzian content of experimental line shapes"
We shall use the simple Pseudo-Voigt-Profil : Voigt ~ η·Lorentz + (1-η)·Gauss
We use the result by "Extended pseudo-Voigt function for approximating the Voigt profile" cf http://www.crl.nitech.ac.jp/~ida/research/reprints/expv.pdf

eta = 1.36603 (G_L/G) - 0.47719 (G_L/G)^2 + 0.11116(G_L/G)^3
where G = ([G_G^5 + 2.69269 G_G^4 G_L + 2.42843 G_G^3 G_L^2 + 4.47163 G_G^2 G_L^3 + 0.07842 G_G G_L^4 + G_L^5]^1/5
for G_G and G_L beeing the FWHM  of the functions: G_G = 2 sqrt (ln 2) sigma_G  and G_L = 2 gamma_L

The worst error is 1% and occurs at the center and for G_G = G_L
**/


// Approximation of the Voigt profile [Lorentz*Gaussian](omega) where Lorentz and Gaussian FWHM values of G_L and G_G and centered on 0.
// From "Extended pseudo-Voigt function for approximating the Voigt profile"
double Pseudo_Voigt(double delta, const double G_L, const double G_G)
{
    double Gfive = pow(G_G,5) + 2.69269*pow(G_G,4)*G_L + 2.42843*pow(G_G,3)*pow(G_L,2) + 4.47163*pow(G_G,2)*pow(G_L,3) + 0.07842*G_G*pow(G_L,4) + pow(G_L,5);
    double G = pow(Gfive,1./5.);
    double eta = 1.36603*(G_L/G) - 0.47719*(G_L/G)*(G_L/G) + 0.11116*(G_L/G)*(G_L/G)*(G_L/G);
    double fG = Gauss(1., delta, G);
    double fL = Lorentz(1., delta, G);
    double Voigt = (1-eta)*fG + eta*fL;
    return Voigt;
}

// comb spectum = gaussian laser but with individual Lorentzian comb lines
// the comb line are positioned at nu_offset + n*nu_repetition
// The intensity is taken at the detuning of the nearest (lorentzian) comb line
// The linewidth is due to the spontaneous emission of the 2 levels system drived by the comb laser
double comb_shape(double I,double delta, const double G_L, const double G_G, const int n_las, const double Energy_transition_cm, FitParams &params)
{
    double IGaussien=I*Pseudo_Voigt(delta, G_L, G_G); // INtensity of the broadband gaussian spectral shape

    std::ostringstream oss;
    oss << "[" << n_las << "]";
    std::string num = oss.str(); // num = "[i]"

// Should be put on laser Class!
// The position of the n^th comb line is nu_offset_MHz + n nu_repetition_MHz and its linewidth is nu_individual_comb_line_MHz + (natural linewidth of the transition)
    double nu_offset_MHz   = params.LocateParam(string("nu_offset_MHz") + num)->val;
    double nu_repetition_MHz = params.LocateParam(string("nu_repetition_MHz") + num)->val;
    double nu_individual_comb_line_MHz = params.LocateParam(string("nu_individual_comb_line_MHz") + num)->val;
    double nu_transition_MHz = 100* C *Energy_transition_cm/MHz;

    int n_comb = double ((nu_transition_MHz - nu_offset_MHz)/nu_repetition_MHz + 0.5) ; // number of the closest comb line. Simple assumption to calcultae the rate. We do not sum over all lines byut just take the closest one!
    double delta_comb =2*pi*MHz*(nu_transition_MHz - (n_comb*nu_repetition_MHz + nu_offset_MHz)); // detuning from the comb line
    double Gamma_line = 2*pi*MHz*nu_individual_comb_line_MHz;// linewidth of individual line
    double I_transition = Lorentz(Gamma_line*pi/2.,  delta_comb, Gamma_line); // Lorentzian normalized to 1 at the center

    return IGaussien*I_transition;
}



// intensity (irradiance) of a pseudo_BBR (pseudo_Black Body Radiation) that is a BBR spectrum but with a gaussian beam
//I(w)=rho(w)*c, et I(w)=Itot*cste*w^3/(exp(hb w/kT) -1) avec Itot= (sigmaT^4)/4
inline double pseudo_BBR_intensity(const double I, const double Energy_transition_cm, const double T)
{
    double integrated_spectrum= pow((kB* pi *T)/(4.*100. * hPlanck *C),4) /15.;//1/4 of (kB^4 Pi^4 T^4)/(15*100^4 h^4 c^4)
    double convolutionFactor = 2* pi * 100. * C; //L(omega) convolué à I(omega) équation B7 Molecular Cooling via Sisyphus
    return (I/integrated_spectrum/convolutionFactor)*Energy_transition_cm*Energy_transition_cm*Energy_transition_cm/(exp(hPlanck*100*C*Energy_transition_cm/(k_Boltzmann*T))-1.); // the spectral desnity is proportional to omega^3/(exp^(hbar omega/kB T) -1) and integrate x^3/(exp-x -1) = pi^4/15
}




// Calcul de l'intensité locale [Lorentz*I](omega) with a natural Cauchy-Lorents linewitdh with FWHM Gamma_atomic = Gamma_spon_fond+Gamma_spon_exc.
// delta = omega_Laser - omega_atomic
double intensity_Convolution_linewidth(const double I, const double delta, const double Gamma_atomic, const double Gamma_Laser, const int laser_type, const int n_las,  const double Energy_transition_cm, const double Energy_transition_laser_cm, FitParams &params)
{
    if (laser_type == lorentzien || laser_type == field_ionization) // for the field ionization it is here just to select roughly the region where the power is not too big
        return Lorentz(I, delta, Gamma_Laser + Gamma_atomic); // Convolution of 2 Lorentzian
    if (laser_type == gaussien || laser_type == pseudo_collimated_top_hat  )
        return I*Pseudo_Voigt(delta, Gamma_atomic, Gamma_Laser); // Voigt = Lorentz convoluted by Gaussian
    if (laser_type == comb)
        return comb_shape(I,delta, Gamma_atomic, Gamma_Laser, n_las, Energy_transition_cm, params); // Voigt = Lorentz convoluted by Gaussian
    if (laser_type == pseudo_BBR)
        return pseudo_BBR_intensity(I, Energy_transition_cm, Energy_transition_laser_cm); // BBR spectra. Temperature is in Energy_transition_laser_cm in Liste_Param
    return 0.;
}



/*************
The rates are explained in the appendix of the article
"Molecular cooling via Sisyphus processes."
*********/

// Excitation rate
double rate_two_level(const double I_loc, const double dip_Debye)
{
    double dipole = dip_Debye*Debye;
    double rate_0 = (dipole/HBAR)*(dipole/HBAR)*I_loc*pi/(epsilon0*C);
    return rate_0;
}

// Rate for (photo-)ionization
double rate_ionization(const double I, const double dip_Debye)
{
    // rate I sigma/h nu = 8 pi^2 alpha <r>^2  I/E_hartree = 8 pi^2 alpha (<d>/e)^2  I/(E_hartree)
    // In principle the cross section is energy dependent but here it is chosen
    // TODO (daniel#4#): treat energy dependent cross section. May be by few energy levels in the continuum and interpolating
    double Gamma = 8.*pi*pi*ALPHA*(dip_Debye*Debye/QE)*(dip_Debye*Debye/QE)*I/hartreeJ; // QE^2=e^2
    return Gamma;
}

// Field ionization rate
double rate_field_ionization(const Laser& my_laser, const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE)
{
    Molecule my_mol = Mol[n_mol];
    int n = Mol.size();

    Vecteur3D r,E,E_coulomb;
    r = my_mol.get_pos();
    E = fieldE.get_Field(r);
    // double Vcoulomb = get_coulomb_potential(Mol, r, n);

    E_coulomb = get_coulomb_field(Mol, n_mol, r); // SI on veux ajouter le champ coulombien
    // get_coulomb_field(const vector <Body> &bod, const int i, const Vecteur3D &ri)
    E = E + E_coulomb;

    double Gamma;
    Gamma = (4900-E.mag())*1e+7;
    return max(Gamma,0.);
}


//If we use BBR, the ionization cross section is given by another formula
// According to (23) Bezuglov 2009 NJP 11 2009 013052
double rate_ionization_BBR(const double I, const Laser& my_laser, const Molecule &my_mol,const double Energy_transition_cm)
{
    double Temp_BBR=cm/my_laser.get_lambda(); // Température codée dans ListParam par energie_cm
    double N_BBR=1/(1.-exp(-hPlanck*100*C*Energy_transition_cm/(k_Boltzmann*Temp_BBR)));
    double n = double (my_mol.exc)-2.47;
    double l = double ((my_mol.deg_number)-(my_mol.two_J))/1000.; //# is 1000l+2j
    //  double rate_ioni_BBR=((k_Boltzmann*Temp_BBR/hbar_Planck)/(pi2*C*C*C))*(2.80/pow(n,7./3.)+(2.09*l*l)/pow(n,11./3.))*log(N_BBR)*conv_ua_smoins1;
    //  return rate_ioni_BBR*(4.*I/(sigmaSB*pow(Temp_BBR,4)));
    //return rate_ioni_BBR;
    return 0.;
}

double rate_ionization_BBR_bis(const double I, const Laser& my_laser, const Molecule &my_mol)
{
    double Temp_BBR_2=cm/my_laser.get_lambda(); // Température codée dans ListParam par energie_cm
    double n_bis = double (my_mol.exc)-2.47;
    double C_log = 1/(1-exp(-157890./(Temp_BBR_2*n_bis*n_bis)));
    double l_bis = double ((my_mol.deg_number)-(my_mol.two_J))/1000.; //# is 1000l+2j
    double C_1 = 14423./pow(n_bis,7./3.);
    double C_2 = (10770.*l_bis*l_bis)/pow(n_bis,11./3.);
    double rate_ioni_BBR_bis = Temp_BBR_2*(C_1 + C_2)*log(C_log);
    return rate_ioni_BBR_bis*(I/(sigmaSB*pow(Temp_BBR_2,4)));
    //return rate_ioni_BBR_bis;
}
// Calculation of the excitation rate for a laser (or a class of coherent laser) for a given transition (bound or ionization)
// Add this rate only is is_rate_recorded is true
double rate_excitation(vector <type_codage_react> &reaction_list, vector <double> &rate,const Laser& my_laser, const vector <Molecule> &Mol, const int n_mol, const int n_las, const Vecteur3D k, const double Itot_loc,
                       const double Itot, const double dipole_debye, const double delta,  const double Energy_transition_cm, const Field &fieldB, const Field &fieldE, const int is_bound_transition,  const Internal_state Internal_state_out, const bool is_rate_recorded) // calculate rate and reaction_list for this transition
{
    Molecule my_mol = Mol[n_mol];
    double rate_exc = 0.; // rate of excitation

    if ((bool) is_bound_transition) // transition bound-bound (0=false)
        rate_exc = rate_two_level(Itot_loc, dipole_debye);
    else // bound_free transition (potential photo_ionization)
        if (my_laser.get_type_laser() == pseudo_BBR) //If we study BBR
            rate_exc= rate_ionization_BBR(Itot, my_laser, my_mol,Energy_transition_cm);
        else if (delta>0 ) // Test if the laser transition is above the ionization continuum. O.K. also for photodetachment
        {
            if (my_laser.get_type_laser() == field_ionization)
                rate_exc = rate_field_ionization(my_laser, Mol, n_mol, fieldB, fieldE);
            else
                rate_exc = rate_ionization(Itot, dipole_debye); // photoionization
        }


    if (rate_exc < SMALL_NUMBER_RATE)
        return 0; // It is not useful to calculate too small rates. So if the rate is too small we do not calculate it

    if (is_rate_recorded == true)
    {
        type_codage_react reaction = { n_mol, n_las, k, {k(0),k(1),k(2)}, k, Internal_state_out};   // {n_mol; n_laser; quant_axis; pol_vector; k_eff_laser; final_internal_state;}. IN this case the quant_axis; pol_vector are not important (we put k but any other values will do )
        reaction_list.push_back( reaction ); // transition de la molécule n_mol avec l'impulstion k_tot vers l'état out
        rate.push_back(rate_exc);
    }
    return rate_exc;
}



/************************************************************************/
/************************** Fonctions transition Molecule ***************/
/************************************************************************/


// Calcul et ajoute de tous les taux d'emission spontanée de la molécule.
int rates_molecule_spon(vector <Internal_state> &Level, vector <type_codage_react> &reaction_list, vector <double> &rate,
                        const Molecule &my_mol, const Field &fieldB, const Field &fieldE, const int num_mol,
                        FitParams &params,  MatrixXcd &H,  MatrixXcd &E0_cm, MatrixXcd &Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[])
{
    Vecteur3D axe_quant,r,v,k_spon;
    r = my_mol.get_pos();
    v = my_mol.get_vel();
    Internal_state Internal_state_in = my_mol ; //  état interne de la molecule

    axe_quant = fieldB.get_Field(r)+ fieldE.get_Field(r); // A priori les deux champs doivent avoir le même axe on choisi. Cela permetra d'orienter l'émission spontanée
// TODO (Daniel#2#): CAREFUL THE QUANTIZATION AXIS IS given by the sum of the fields!!!! Should not be like this should be the dominant field!!! Check elsewhere also

    Vecteur3D Bfield;
    Bfield= fieldB.get_Field(r);
    double B = Bfield.mag();
    double E = fieldE.get_Field(r).mag();
    double Energy_transition_cm =0.;

    double Gamma=0.; // Spontaneous emission rate
    double k=0.; // Spontaneous emission lenght wave vector


    /**********  WE DIAGONALIZED THE ENERGY LEVELS so the transition dipole moment are not constant and need to be calculated
    BUt we do not know the real k vector of the sponteneous emission so we do not correct for the recoil energy level **********/

    if ( (params.LocateParam("is_Levels_Lines_Diagonalized")->val) )
    {
        Diagonalization(Level, my_mol, fieldB, fieldE, params, H, E0_cm, Zeeman_cm_B, d0, d);
        // diagonalized the Hamiltionian for B field and v velocity and give the eigenvectors and eigenvalues and  update all Level[n].Energy_cm

        int i = my_mol.deg_number; // The molecules is in the Level number n_level_in.// so Level[ # = deg_number] shall be the Level itself// So in the Level file the deg_number is the Level number (START FROM 0)
        Internal_state_in = Level[i]; // Internal_state_in = my_mol ; //  état interne de la molecule

        for (int j=0; j<i; j++) // We scan only for levels that are below level i in energy. So we do a transition i-->j
        {
            Gamma =0.;
            Energy_transition_cm = Level[i].Energy_cm - Level[j].Energy_cm;
            k = 2*pi*100.*Energy_transition_cm;
            for(int n_polar = -1; n_polar <= 1; n_polar++) // Gamma_ij propto |<i | d | j>|^2  = Sum_q |<i | d_q | j>|^2
            {
                complex<double> dipole_ij = d[n_polar+1](i,j);  //   dipole_ij  =   <i | d_q | j>
                // <i|d_q|j> is coded here as d[q+1][i][j] , i = line, j = column ; that is for a i<-->j transition (with E_i> E_j)
                double dip_square= norm(dipole_ij); // norm(z)=|z|^2
                Gamma +=  dip_square * Spol_Debye_A_s * Energy_transition_cm*Energy_transition_cm*Energy_transition_cm;
            }

            if (Gamma < SMALL_NUMBER_RATE)
                continue; // Pour une transition trop faible inutile de la calculer

            k_spon =  Vecteur3D(k,0.,0.) ; // Will be used to give the proper photon recoil
            rate.push_back(Gamma); // rate[nb_rate]=Gamma
            type_codage_react reaction = { num_mol, -1, axe_quant, {d[0](i,j),d[1](i,j),d[2](i,j)}, k_spon, Level[j]}; // {n_mol; n_laser; quant_axis; pol_vector; k_eff_laser; final_internal_state;}.
            reaction_list.push_back (reaction);
            // numéro laser -1=spon_emission pour émission spontanée.
            // for the spontaneous emission we give  the polarization vector to determine the emitted photon (not the  impulsion photon as for stimulated or absorption)
            // So we give the polarization vector (dipole transition for polar -1, 0, +1 ) in the local quantization axis frame axis (the one where OZ = axe_quant)

        }
    }

    /********** THE ENERY LEVELS ARE not diagonalized and so the transition dipole moment are constant **********/

    else
    {
        for( int i = 0; i < (int) my_mol.liste_raies.size(); i++ ) // Boucle sur les transitions accessibles
        {
            Internal_state Internal_state_out = *(my_mol.liste_raies[i].first) ; //  état interne de la molecule après la réaction
            double Energy_in = Internal_state_in.Energy0_cm + (Internal_state_in.Energy_Shift_B_cm(B) + Internal_state_in.Energy_Shift_E_cm(E));
            double Energy_out = Internal_state_out.Energy0_cm + (Internal_state_out.Energy_Shift_B_cm(B) + Internal_state_out.Energy_Shift_E_cm(E));
            Energy_transition_cm = Energy_in - Energy_out; // TODO (dc#4#): On ne prend pas le déplacement lumineux. ...Il faudrait revoir sérieusement cet affaire de Energie mise à jour

            if (Energy_transition_cm < 0.)
                continue; // Etat en dessous de l'autre donc pas de photon spontané

            double dipole_debye = my_mol.liste_raies[i].second; // dipole de la transition en debye
            Gamma= Gamma_spon(dipole_debye, Energy_transition_cm);

            if (Gamma < SMALL_NUMBER_RATE)
                continue; // Pour une transition trop faible inutile de la calculer

            k = 2*pi*100.*Energy_transition_cm;
            k_spon =  Vecteur3D(k,0.,0.) ; // Will be used to give the proper photon recoil

            rate.push_back(Gamma); // rate[nb_rate]=Gamma
            type_codage_react reaction = { num_mol, -1, axe_quant, {axe_quant(0),axe_quant(1),axe_quant(2)}, k_spon, Internal_state_out};  // {n_mol; n_laser; quant_axis; pol_vector; k_eff_laser; final_internal_state;}.
            // IN this case the pol_vector, is  not important (the photon will be calculated by the Delta_M given by the Internal_states) so we put axe_quant but any other values will do )
            reaction_list.push_back (reaction); // numéro laser -1=spon_emission pour émission spontanée. Le vecteur impulsion sera calculé plus tard
        }
    }

    return rate.size();
}




// Calcul of rate (if no interference between laser) between level in and out for a given laser and for a given molecule. This add the light shift effect to the dipolar potential (delta_pot_dipolaire)
// update the detuning (delta), the polarization vector also updated), Gamma_spot_total, sqrt of the local intensity of the laser
int rates_single_molecule_laser_level(const int n_las, double dipole, double &delta, double &Gamma_spon_tot, double sqrt_intensity_loc[],
                                      vector <Internal_state> &Level, Internal_state &Internal_state_in, Internal_state &Internal_state_out, vector <type_codage_react> &reaction_list,
                                      vector <double> &rate, const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const Laser &my_laser, const double t, double &delta_pot_dipolaire,
                                      FitParams &params,  MatrixXcd d[], bool is_rate_recorded, int is_bound_transition, const int n_level_in, const int n_level_out,
                                      const double Gamma_in)
{
    Molecule my_mol= Mol[n_mol];

    if ( (params.LocateParam("is_Levels_Lines_Diagonalized")->val) )
    {
        Internal_state_out = Level[n_level_out];
        is_bound_transition = Internal_state_out.bound_level; // To know if the transition is towards continuum (or annihilation) for instance
    }


    Vecteur3D k,v,r;
    k = my_laser.wave_vector();
    r = my_mol.get_pos();
    v = my_mol.get_vel();

    double Gamma_Las = my_laser.get_Gamma_Laser();  // Largeur spectrale FWHM  (en s^-1)

    // we treat absorption and emission  together.
    double Energy_transition_cm = abs(Internal_state_out.Energy_cm - Internal_state_in.Energy_cm)  + k*v/Conv_Ecm_delta; // Il faut mettre le detuning  car c'est la fréquence vue par les particules et donc par le spectre du laser.
    double Energy_transition_laser_cm = cm/my_laser.get_lambda(); // Energie de la transition laser en cm^-1
    delta =(Energy_transition_laser_cm  - Energy_transition_cm)*Conv_Ecm_delta  ;// detuning de la transition en s^-1. delta = omega_L - k.v - omega_transition
    double I_laser_tot = my_laser.intensity_lab_axis(r) * my_laser.intensity_t_nanosecond(t/ns) *  my_laser.transmission_spectrum(Energy_transition_cm); // Intensité laser à la position de la molécule au temps t

    Vecteur3D axe_quant;
    complex<double> dipole_vector[3]= {0.,0.,0.};
    axe_quant = fieldB.get_Field(r)+ fieldE.get_Field(r); // A priori les deux champs doivent avoir le même axe on choisi. Cela permetra d'orienter l'émission spontanée

    double rate_exc=0.;
    double dipole_debye_trans_including_laser =0.;

//  dipole_q_ij <i|d_q|j> is coded here as d[q+1][i][j] (complex), i = line, j = column ; that is for a i<-->j transition (WITH E_i> E_j)
    if ( (params.LocateParam("is_Levels_Lines_Diagonalized")->val) )
    {
        for(int q = -1; q <= 1; q++) // scan over the polarization d^(q) sur -1, 0, +1  ;d[0]=dsigma-, d[1] =dpi , d[2]=dsigma+ for polarization in absorption ; d0[q+1]_ij = 0_<i | d^(q) | j>_0
        {
            // We choose the proper dipole because d[q+1][i][j] (WITH E_i> E_j to get m_i = m_j+q
            if (Internal_state_in.Energy_cm > Internal_state_out.Energy_cm) // EMISSION  up=i=n_level_in   low=j=n_level_out
                dipole_vector[q+1] =  d[q+1](n_level_in,n_level_out);
            else   // ABSORPTION up=i=n_level_out   low=j=n_level_in
                dipole_vector[q+1] =  d[q+1](n_level_out,n_level_in);
        }
        dipole =  sqrt(norm(dipole_vector[0]) + norm(dipole_vector[2]) + norm(dipole_vector[2]) );  // dipole of the transition.  norm(z)=|z|^2

        Gamma_spon_tot =   Gamma_in +  Gamma_Level_from_diagonalized_dipole(Level, d, n_level_out);
    }
    else // No diagonalization but fixed transition strength and a give Delta M transition
    {
        int two_M_up,two_M_low; // Projection du moment angulaire
        if (Internal_state_in.Energy_cm > Internal_state_out.Energy_cm) // EMISSION up=i=n_level_in   low=j=n_level_out
        {
            two_M_up = Internal_state_in.two_M;
            two_M_low = Internal_state_out.two_M;
        }
        else   // ABSORPTION
        {
            two_M_low = Internal_state_in.two_M;
            two_M_up = Internal_state_out.two_M;
        }
        dipole_vector[(two_M_up-two_M_low)/2+1]=dipole; // dipole_vector[M_up-M_low+1] est la composante du champ ayant la bonne polarisation pour la transition que ce soit en absorption ou émission
        Gamma_spon_tot = Internal_state_in.one_over_lifetime + Internal_state_out.one_over_lifetime; // Spontaneous decay on the transition Sum_k Gamma_k
    }


    //  the rate is proportional to  dipole*dipole*I noted dipole*dipole*Itot_loc in rate_excitation function
    // but in reality comes from |<i| d.E|j>|^2 =  |sum_q  E^(q) <i| d_(q) |j>|^2 = E^2  |sum_q  eps^(q) <i| d_(q) |j> |^2.
    // So destructive interference can appears and so the effectif_dipole_local has to be included in the dipole


    dipole_debye_trans_including_laser = effectif_dipole_local(dipole_vector, axe_quant, my_laser);  // (asbolute value of) effectif dipole d.e_laser = sum_p d_p epsilon^p
    double eps_pol = dipole_debye_trans_including_laser/dipole; // Part of the laser amplitude that excited the transition


    sqrt_intensity_loc[n_las] =  eps_pol*sqrt(intensity_Convolution_linewidth(I_laser_tot, delta, Gamma_spon_tot, Gamma_Las, my_laser.get_type_laser(),n_las,
                                 Energy_transition_cm, Energy_transition_laser_cm, params)); // Proportional to the Rabi Frequency -d.E (we include the laser polarisation here and not in the dipole)


    if (dipole_debye_trans_including_laser < SMALL_DIPOLE_DEBYE )  // Pour une transition trop faible inutile de calculer le taux!
        return rate.size();

    if (my_laser.get_coherent_avec_laser_num() == -1)  // Laser seul, pas d'interférence avec les autres, on le calcul maintenant pour ne pas le recalculer ensuite
    {
        double Itot_loc =  sqrt_intensity_loc[n_las]* sqrt_intensity_loc[n_las];
        //    if ( (Itot_loc < SMALL_NUMBER)  && !((bool) is_bound_transition) )
        // TEST FOR SPEED IF NEEDED
        //        return rate.size();; // Pour une transition trop faible inutile de calculer le taux! I removed this because for continuum the detuning is big and so the "intensity" is small but should be calculated
        rate_exc = rate_excitation(reaction_list, rate,my_laser, Mol, n_mol, n_las, k, Itot_loc, I_laser_tot, dipole, delta,Energy_transition_cm, fieldB, fieldE, is_bound_transition, Internal_state_out, is_rate_recorded); // calculate rate and reaction_list for this transition. If we do not wat to calculate the rate we just use reaction_list[0] which is *reaction_list
        delta_pot_dipolaire += sgn(Internal_state_out.Energy_cm - Internal_state_in.Energy_cm) * rate_exc * delta/Gamma_spon_tot;
    }
    return rate.size();
}




// Forced_rate for any given purpose (repumoing, Sisyphus, ...)
int forced_rates_single_molecule_laser_level(const int n_las, double dipole, double &delta, double &Gamma_spon_tot, double sqrt_intensity_loc[],
        vector <Internal_state> &Level, Internal_state &Internal_state_in, Internal_state &Internal_state_out, vector <type_codage_react> &reaction_list,
        vector <double> &rate, const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const Laser &my_laser, const double t, double &delta_pot_dipolaire,
        FitParams &params,  MatrixXcd d[], bool is_rate_recorded, int is_bound_transition, const int n_level_in, const int n_level_out,
        const double Gamma_in)
{
    Molecule my_mol= Mol[n_mol];


    Vecteur3D k,v,r;
    k = my_laser.wave_vector();
    r = my_mol.get_pos();
    v = my_mol.get_vel();


    double Energy =  Internal_state_in.Energy_cm;

    double rate_exc=0.;

    Vecteur3D B;

    B = fieldB.get_Field(r); // A priori les deux champs doivent avoir le même axe on choisi. Cela permetra d'orienter l'émission spontanée

    // Transfer MF = +1 --> -1 si B>1T
    if(my_mol.two_M==1 && B.mag()>1 )
        rate_exc = 1e6;

    // Transfer MF = -1 --> +1 si B<0.1T
    if(my_mol.two_M==-1 && B.mag()<0.2 )
       rate_exc = 1e6;

    type_codage_react reaction = { n_mol, n_las, k, {k(0),k(1),k(2)}, k, Internal_state_out};   // {n_mol; n_laser; quant_axis; pol_vector; k_eff_laser; final_internal_state;}. IN this case the quant_axis; pol_vector are not important (we put k but any other values will do )
    reaction_list.push_back( reaction ); // transition de la molécule n_mol avec l'impulstion k_tot vers l'état out
    rate.push_back(rate_exc);



    return rate.size();
}





// Calcul de tous les taux de la molécule  pour tous les lasers
// On les ajoutes aux taux des autres
//  return nb_rate

int rates_molecule(vector <Internal_state> &Level, vector <type_codage_react> &reaction_list,  vector <double> &rate, const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
                   const double t, double &delta_pot_dipolaire, FitParams &params,
                   MatrixXcd &H, MatrixXcd &E0_cm, MatrixXcd &Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[], bool is_rate_recorded)
{
    Molecule my_mol = Mol[n_mol];
    delta_pot_dipolaire = 0.;
    double Gamma_spon_tot=0.;
    double dipole_debye =0.;
    double rate_exc=0;
    double sqrt_intensity_loc[100], delta[100]; // intensities ~sqrt(I), détuning
    int is_bound_transition = 1 ; // Trick to treat the ionization or annihilation (in this case it is 0).
    Internal_state Internal_state_in,Internal_state_out; // State for transitions

    Vecteur3D r,v,Bfield;
    r = my_mol.get_pos();
    v = my_mol.get_vel();
    Bfield= fieldB.get_Field(r);
    double B = Bfield.mag();
    double E = fieldE.get_Field(r).mag();

    const int Nb_laser = laser.size();

    if (Nb_laser>=100)
        cerr << "Nb de laser trop grand: modifier la taille du tableau et recompiler " << endl;


    /**********  WHEN WE DIAGONALIZED THE ENERGY LEVELS,the transition dipole moment are not constant and need to be calculated (we do not treat interference between lasers) **********/

    if ( (params.LocateParam("is_Levels_Lines_Diagonalized")->val) )
    {
        Diagonalization(Level, my_mol, fieldB, fieldE, params, H, E0_cm, Zeeman_cm_B, d0, d); // calcul of the new dipole transition for the new levels.

        int n_level_in = my_mol.deg_number; //The molecules is in the Level number n_level_in.// so Level[ # = deg_number] shall be the Level itself// So in the Level file the deg_number is the Level number (START FROM 0)
        Internal_state_in = Level[n_level_in]; // Internal_state_in = my_mol ; //  état interne de la molecule
        double Gamma_in = Gamma_Level_from_diagonalized_dipole(Level, d, n_level_in); // Initial spontaneous decay rate of the state

        for (int n_las = 0 ; n_las < Nb_laser; n_las++) // On scan sur les lasers pour calculer l'intensité sur la transtion // If we switch we add Nb_laser to the number
        {
            Laser my_laser = laser[n_las];
            Vecteur3D k;
            k = my_laser.wave_vector();
            double m=my_mol.get_mass();


            /*** We assume that there is no level crossing with the initial level so that absorption and emission are well define by the energy ordering . But the upper or lower level ordering may change due to hbar k modification  ****/

            /***** Stimulated emission: v--> v-HBAR*k/m *****/


            if (Internal_state_in.bound_level == 1) // If the level is not a bound level it can not cahneg (continnumm, dead or else)
            {
                my_mol.set_vel(v-HBAR*k/m); // new volocity for the diagonalization
                Diagonalization(Level, my_mol, fieldB, fieldE, params, H, E0_cm, Zeeman_cm_B, d0, d);// calcul of the new dipole transition for the new levels (after the  emission of photons). The energy levels order may have changed
                my_mol.set_vel(v); // put back as normal value
            }
            for( int n_level_out = 0; n_level_out < n_level_in; n_level_out++ )  // Emission so we scan only on level below
            {
                rates_single_molecule_laser_level( n_las, dipole_debye, delta[n_las],  Gamma_spon_tot, sqrt_intensity_loc, Level, Internal_state_in, Internal_state_out,
                                                   reaction_list, rate, Mol, n_mol, fieldB, fieldE, my_laser, t, delta_pot_dipolaire,
                                                   params, d, is_rate_recorded, is_bound_transition, n_level_in,  n_level_out, Gamma_in);
            }

            /***** Absorption: v--> v+HBAR*k/m *****/


            if (Internal_state_in.bound_level == 1) // If the level is not a bound level it can not cahneg (continnumm, dead or else)
            {
                my_mol.set_vel(v+HBAR*k/m); // new volocity for the diagonalization
                Diagonalization(Level, my_mol, fieldB, fieldE, params, H, E0_cm, Zeeman_cm_B, d0, d);// calcul of the new dipole transition for the new levels (after the  emission of photons). The energy levels order may have changed
                my_mol.set_vel(v); // put back as normal value
            }
            for( int n_level_out = n_level_in+1; n_level_out < (int)  Level.size(); n_level_out++ ) // Absorption  so we scan only on level above
            {
                is_bound_transition = Level[n_level_out].bound_level;
                rates_single_molecule_laser_level(n_las, dipole_debye, delta[n_las], Gamma_spon_tot, sqrt_intensity_loc, Level, Internal_state_in, Internal_state_out,
                                                  reaction_list, rate, Mol, n_mol, fieldB, fieldE, my_laser, t, delta_pot_dipolaire, params, d,
                                                  is_rate_recorded, is_bound_transition, n_level_in,  n_level_out, Gamma_in);
            }
        }

// TODO (Daniel#1#): Do interference with diagonalization

        /** Calcul of rate due to interference between lasers **/

        /*   for( int n_level_out = 0; n_level_out < (int)  Level.size(); n_level_out++ )
          {

              if (Internal_state_in.bound_level == 1) // If the level is not a bound level it can not cahneg (continnumm, dead or else)
              {
                  rates_single_molecule_lasers_interference_level(Level[n_level_out], reaction_list, rate, Mol, n_mol, fieldB, fieldE, laser, sqrt_intensity_loc,
                          Internal_state_in, Internal_state_out, delta, dipole_debye, delta_pot_dipolaire, Gamma_spon_tot, t, is_bound_transition, is_rate_recorded);
              }
          }
        */

    }

    /********** THE ENERGY LEVELS ARE not diagonalized and so the transition dipole moment are constant **********/

    else
    {
        Internal_state_in = my_mol ; //  état interne de la molecule
        Internal_state_in.Energy_cm  = Internal_state_in.Energy0_cm + (Internal_state_in.Energy_Shift_B_cm(B) + Internal_state_in.Energy_Shift_E_cm(E));

// If looking with debugger here this create Segmentation fault. this is due (but not understood !) by the copy my_mol=Mol[n_mol] and look at my_mol.liste_raies.size()
        for( int n_trans = 0; n_trans < (int) my_mol.liste_raies.size(); n_trans++ ) // Boucle sur les transitions accessibles
        {
            dipole_debye = my_mol.liste_raies[n_trans].second ; // dipole de la transition en debye
            Internal_state_out = *(my_mol.liste_raies[n_trans].first) ; //  état interne de la molecule après la réaction
            is_bound_transition = Internal_state_out.bound_level;

            if (dipole_debye < SMALL_DIPOLE_DEBYE)
                continue; // Pour une transition trop faible inutile de calculer le taux! (les "dipole" sont des "petite section pour la photoinizaiton donc on fait exception)

            Internal_state_out.Energy_cm = Internal_state_out.Energy0_cm + (Internal_state_out.Energy_Shift_B_cm(B) + Internal_state_out.Energy_Shift_E_cm(E));
// The dipolar potential is not included in the shift for the transition. ...
// TODO (Daniel#4#): This avoids accumulation (not so clear what dipolar detuning need to be included, if many laser and wavelength), but for some case it may be good to have it (to be checked)
//

            for (int n_las = 0 ; n_las < Nb_laser; n_las++) // On scan sur les lasers pour calculer l'intensité sur la transtion // If we switch we add Nb_laser to the number
            {
                Laser my_laser = laser[n_las];

                //  Calcul les taux dans le cas où il n'y a pas d'interférence lasers
                if (my_laser.get_coherent_avec_laser_num() == -1)  // Laser seul, pas d'interférence avec les autres
                {
                    rates_single_molecule_laser_level( n_las, dipole_debye, delta[n_las], Gamma_spon_tot, sqrt_intensity_loc, Level,  Internal_state_in, Internal_state_out, reaction_list,
                                                       rate, Mol, n_mol, fieldB, fieldE, laser[n_las], t, delta_pot_dipolaire, params, d, is_rate_recorded, is_bound_transition);
// Calcul of rate (if no interferece between lasers)  between level in and out for a given laser and for a given molecule. This add the light shift effect to the dipolar potential (delta_pot_dipolaire)
// update the detuning (delta), the polarization vector also updated), Gamma_spot_total, sqrt of the local intensity of the laser

                    if ( (params.LocateParam("is_forced_rates")->val) ) // Forced rates if exists
                        forced_rates_single_molecule_laser_level( n_las, dipole_debye, delta[n_las], Gamma_spon_tot, sqrt_intensity_loc, Level,  Internal_state_in, Internal_state_out, reaction_list,
                                rate, Mol, n_mol, fieldB, fieldE, laser[n_las], t, delta_pot_dipolaire, params, d, is_rate_recorded, is_bound_transition);

                }
                else
                    rates_single_molecule_laser_level( n_las, dipole_debye, delta[n_las], Gamma_spon_tot, sqrt_intensity_loc, Level,  Internal_state_in, Internal_state_out, reaction_list,
                                                       rate, Mol, n_mol, fieldB, fieldE, laser[n_las], t, delta_pot_dipolaire, params, d, false, is_bound_transition);
// Mainly needed to calculate sqrt of the local intensity of the laser for the interference case below
            }

            /** Calcul of rate due to interference between lasers **/
            rates_single_molecule_lasers_interference_level(Level, reaction_list, rate, Mol, n_mol, fieldB, fieldE, laser, sqrt_intensity_loc,
                    Internal_state_in, Internal_state_out, delta, dipole_debye, delta_pot_dipolaire, Gamma_spon_tot, t, is_bound_transition, is_rate_recorded);
        }
    }
    return rate.size();
}



// Calcul of rate (useful only if interference between lasers) and for a given molecule on  given transition
// update the dipolar shift
void rates_single_molecule_lasers_interference_level(vector <Internal_state> &Level, vector <type_codage_react> &reaction_list,  vector <double> &rate,
        const vector <Molecule> &Mol, const int n_mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser,
        double sqrt_intensity_loc[], Internal_state &Internal_state_in, Internal_state &Internal_state_out,  double delta[], double dipole_debye,
        double &delta_pot_dipolaire, double Gamma_spon_tot, double t, bool is_bound_transition,  bool is_rate_recorded)
{
    const int Nb_laser = laser.size();
    Molecule my_mol = Mol[n_mol];
    Vecteur3D r,v;
    r = my_mol.get_pos();
    v = my_mol.get_vel();

    /** Il y a des interférences: on va calculer l'interférence ainsi que l'absorption
    Itot_loc c'est proportionel à |Ω|^2= Σ j,i |Ωi||Ωj|cos[(ωi − ωj )t − (ki − kj).(r + vt)]
    I_ktot c'est proportionel à |Ω|^2 ktot= Σ j,i (ki+kj)/2 |Ωi||Ωj|cos[(ωi − ωj )t − (ki − kj).(r + vt)] donc ktot = I_ktot/Itot

    We could use the total Rabi Frequency but when the laser is shaped it is not easy. So we prefer to use the local intensity (that is ~I/Gamma)
    So |Ω| is proportional to sqrt(I_tot_loc) and |Ωi| to sqrt(I_loc) with the same proportionality factors (for instance (Gamma_tot/2pi)/((Gamma_tot/2)^2 + delta^2) ) so it will be remove at the end.
    This is true because all lasers are similar (if not we will have no interferences!)
    ***/

    double omegai,omegaj,Itot_loc;
    Vecteur3D ki,kj,I_ktot;
    for (int n_classe_coh_las =  0; n_classe_coh_las < Nb_laser; n_classe_coh_las++) // On scan sur les classes de lasers cohérents. un laser j interfère avec un laser i<=j et le plus petit i interfère avec lui même!
    {
        int num_laser_coherent = laser[n_classe_coh_las].get_coherent_avec_laser_num(); // Numéro de la classe (en fait du premier laser de la classe de cohérence)
        if (n_classe_coh_las != num_laser_coherent)
            continue; // Ce laser n'est pas le premier de la classe. On a donc déjà calculé cette classe on passe au suivant!
        Itot_loc = 0.;
        I_ktot = Vecteur3D(0.,0.,0.);
        double Energy_transition_class_cm = abs(Internal_state_out.Energy_cm - Internal_state_in.Energy_cm) + laser[num_laser_coherent].wave_vector()*v/Conv_Ecm_delta ;// k.v effet Doppler (à resonnance)

        for (int i = n_classe_coh_las; i < Nb_laser; i++) // On scan sur les lasers de la classe (qui sont ordonnés donc n_classe_coh_las est le premier
        {
            if (laser[i].get_coherent_avec_laser_num() != num_laser_coherent)
                continue ; // Si pas cohérent on stoppe avec ce laser et on continue avec un autre
            double Itoti = sqrt_intensity_loc[i]*sqrt_intensity_loc[i];
            ki = laser[i].wave_vector();
            omegai = laser[i].get_omega();
            I_ktot += ki*Itoti;
            Itot_loc += Itoti;

            /****
            Calcul de l'intensité locale liée à l'interférence
            |Ω|^2= Σ L′,L |ΩL||ΩL'|cos[(ωL − ωL′ )t − (kL − kL′ ).(r + vt)] avec Ω^2 ~I
            Le moment transféré est calculé par
            |Ω|^2 ktot= Σ j,i (ki+kj)/2 |Ωi||Ωj|cos[(ωi − ωj )t − (ki − kj).(r + vt)]
            et Σ j,i = Σi + 2 Σ j>i
            Remark: In the case of lattice, we can have an absorption rate but no force. That mean that the momentum transfer at the absorption is zero.
            But we still have the "random" spontaneous emission existing. I think this is how it works in reality!
            It calculate also the dipole detuning in s^-1 which is just rate * Gamma/delta.
            In this case we do not necessarily to calculate the rate, for instance when calculate only dipolar shift, so is_rate_recorded would be false
            **/


            for (int j = i+1; j < Nb_laser; j++) // On scan les autres lasers j>i et on ne va garder que ceux cohérent avec i
            {
                if (laser[j].get_coherent_avec_laser_num() != num_laser_coherent)
                    continue ; // Si pas cohérent on stoppe avec ce laser et on continue avec un autre
                kj = laser[j].wave_vector();
                omegaj = laser[j].get_omega();
// We should put the phase given by the polarizations sometimes to take into account linear polarisation
                double sqrt_Ii_Ij =  sqrt_intensity_loc[i]*sqrt_intensity_loc[j]*cos((omegai-omegaj)*t - (ki-kj)*r);
                Itot_loc += 2.*sqrt_Ii_Ij;
                I_ktot += 2.*sqrt_Ii_Ij*(ki+kj)/2.; //
            }
        }
        double I_laser_tot = laser[num_laser_coherent].intensity_lab_axis(r) * laser[num_laser_coherent].intensity_t_nanosecond(t/ns)  * laser[num_laser_coherent].transmission_spectrum(abs(Internal_state_out.Energy_cm - Internal_state_in.Energy_cm)); // Intensité laser à la position de la molécule. WE DO NOT USE THE INTERFERENCE HERE
        double I_tot = I_laser_tot*Itot_loc/(sqrt_intensity_loc[num_laser_coherent]*sqrt_intensity_loc[num_laser_coherent]); // phenomenoligue way to treate I_tot. We take the same ratio for the total laser intensity that for the local one!
        if (I_laser_tot < SMALL_NUMBER)
            continue; // Pour une transition trop faible inutile de calculer le taux!
        double rate_exc = rate_excitation(reaction_list, rate,laser[n_classe_coh_las], Mol, n_mol, num_laser_coherent, I_ktot/Itot_loc, Itot_loc, I_tot, dipole_debye, delta[num_laser_coherent], Energy_transition_class_cm,
                                          fieldB, fieldE, is_bound_transition, Internal_state_out, is_rate_recorded); // calculate rate and reaction_list for this laser class and transition
        delta_pot_dipolaire += sgn(Internal_state_out.Energy_cm - Internal_state_in.Energy_cm) * rate_exc * delta[num_laser_coherent]/Gamma_spon_tot;
    }
}


// return the decay rate from Level i = sum_j<i Gamma_ij
double Gamma_Level_from_diagonalized_dipole(vector <Internal_state> &Level, MatrixXcd d[], const int i)
{
    double  Gamma=0.;

    for (int j=0; j<i; j++) // We scan only for levels that are below level i in energy
    {
        for(int n_polar = -1; n_polar <= 1; n_polar++)
        {
            complex<double> dipole_ij;
            dipole_ij = d[n_polar+1](i,j);
            double dip_square= norm(dipole_ij); // norm[z]=|z|^2
            double DE_cm = Level[i].Energy_cm - Level[j].Energy_cm;
            Gamma +=  dip_square * Spol_Debye_A_s * DE_cm*DE_cm*DE_cm;
        }
    }
    return Gamma;
}



// Copie tous les taux sauf celui de la molécule numero_mol
// Retourne le nombre de taux de ces molécules
// Could be done more efficiently using swap or other!
int copie_rates_molecules(vector <type_codage_react> &reaction_list, vector <double> &rate, const int numero_mol, const int nombre_old_rate)
{

    int n_mol; // Numéro de la molécule
    int nb_rate = 0;

    for (int n_reac = 0; n_reac < nombre_old_rate; n_reac++)
    {
        n_mol = reaction_list[n_reac].n_mol; // Numéro de la molécule affectée
        if (n_mol == numero_mol)
            continue; // Fait la boucle sauf si c'est la molécule numero_mol
        rate[nb_rate] = rate[n_reac]; // A mon avis il n'y a pas de risque car nb_rate <= n_reac
        reaction_list[nb_rate] = reaction_list[n_reac];
        nb_rate ++;
    }

    rate.erase(rate.begin()+nb_rate,rate.end()); // Remove the other rates rate[nb_rate], rate[nb_rate + 1] .... So we keep all rates except the one form the n_mol
    reaction_list.erase(reaction_list.begin()+nb_rate,reaction_list.end());
    return nb_rate;
}


// Calcul de tous les taux de toutes les molécules si numero_mol = aucune
// Sinon on ne recalcule que celui de la molécule numero_mol
int calcul_rates_molecules(vector <Internal_state> &Level, MC_algorithmes Algorithme_MC, vector <type_codage_react> &reaction_list, vector <double> &rate, const vector <Molecule> &Mol, const Field &fieldB, const Field &fieldE, const vector <Laser> &laser, const double t,
                           const int numero_mol, const int N_Mol, FitParams &params,
                           MatrixXcd &H, MatrixXcd &E0_cm, MatrixXcd &Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[])
{
    if (Algorithme_MC == Aucun_MC)
    {
        rate.clear();
        reaction_list.clear();
        return -1;
    }

    double delta_pot_dipolaire = 0.; // not used here but needed to call the functions

    if ((numero_mol != aucune)&&(Mol[numero_mol].bound_level == 1)) // Une molécule a été affectée donc on ne recalcule que ses taux, les autres sont recopiés
        //  FOr a state in the continuum there is no transition (back) it stay where it is. So we do not calculate anything
    {
        copie_rates_molecules(reaction_list, rate, numero_mol, rate.size());
        rates_molecule_spon(Level, reaction_list, rate, Mol[numero_mol], fieldB, fieldE, numero_mol, params, H, E0_cm, Zeeman_cm_B, d0, d); // Emission spontanée
        rates_molecule(Level, reaction_list, rate, Mol, numero_mol, fieldB, fieldE, laser, t, delta_pot_dipolaire,  params, H, E0_cm, Zeeman_cm_B, d0, d);        // Absorption ou Emission stimulée
    }
    else  //  Si on a (forcé) numero_mol = aucune on recalcule tout
    {
        rate.clear();
        reaction_list.clear();
        for (int n_mol = 0; n_mol < N_Mol; n_mol++)
        {
            if (Mol[n_mol].bound_level == 1) // 0=continuum, 1 = bound level //  FOr a state in the continum there is no transition (back) it stay where it is. So we do not calculate anything
            {
                rates_molecule_spon(Level, reaction_list, rate, Mol[n_mol], fieldB, fieldE, n_mol, params,  H, E0_cm, Zeeman_cm_B, d0, d); // Emission spontanée
                rates_molecule(Level, reaction_list, rate, Mol, n_mol, fieldB, fieldE, laser, t, delta_pot_dipolaire,  params,  H, E0_cm, Zeeman_cm_B, d0, d);  // Absorption ou Emission stimulée
            }
        }
    }

    return rate.size();
}


/***********************************************/
/********** réalisation d'une réaction *********/
/***********************************************/


// We scan all molecules and for each of them calculate one random
// Scan all rates and realize the reaction if the rate*random < dt=0.1/rate_max
// Typically  a tens of Molecules will be affected
int do_reaction_FastRoughMethod(const MC_algorithmes Algorithme_MC, const gsl_rng * r, const vector <double> &rate, const vector <type_codage_react> &reaction_list, vector <Molecule> &Mol,
                                const int n_reac, const vector <Laser> &laser, const double dt_KMC, ofstream & file_rate, int &number_photons, FitParams &params)
{
    int  nb_rate = rate.size();
    int nb_molecule_affected = -1; // number of the molecule that has undergone a reaction (none at the beginning of this loop)

    double random = 1.-gsl_rng_uniform(r); //  uniform random number u \in (0,1]

    for (int n_reac = 0; n_reac < nb_rate; n_reac++) // Scan all the rates (for all molecules)
    {
        int numero_mol = reaction_list[n_reac].n_mol; // Number for molecule to be affected
        if (numero_mol == nb_molecule_affected)
            continue; // we do only one rate per molecule to avoid problems

        double Delta_t = (- log(random))/rate[n_reac];    // Typical time for this reaction to occur
        if (Delta_t < dt_KMC) // The reaction occurs only if Gamma dt (random) < 0.1. Calculated like:  (random)/ Gamma <  0.1/rate_Max.
        {
            do_reaction(Algorithme_MC, r, rate, reaction_list, Mol, n_reac, laser, dt_KMC, file_rate, false, number_photons, params); // Do the reaction as if it was
            nb_molecule_affected = numero_mol;
            random = 1.-gsl_rng_uniform(r); // New random number for the other molecule
        }
    }
    return aucune; // In a sens "Aucune" Molecule has not been affected. So "all" (a lot) of Molecules have been affected and thus we will have to calculate all potential again
}

// Makes the reaction that has been chosen: reaction {n_mol; n_laser; quant_axis; pol_vector; k_eff_laser; final_internal_state;}.
int do_reaction(const MC_algorithmes Algorithme_MC, const gsl_rng * r, const vector <double> &rate, const vector <type_codage_react> &reaction_list, vector <Molecule> &Mol,
                const int n_reac, const vector <Laser> &laser, const double dt_KMC, ofstream & file_rate,  bool first_call, int &number_photons, FitParams &params )
{

    if (n_reac == -1)
        return -1; // No Reaction like in Aucun_MC
    if (Algorithme_MC == Fast_Rough_Method && first_call == true)
        do_reaction_FastRoughMethod(Algorithme_MC, r, rate, reaction_list, Mol, n_reac, laser, dt_KMC, file_rate, number_photons, params); // Do the reaction as if it was

    int numero_mol = reaction_list[n_reac].n_mol; // Numéro de la molécule affectée
    int numero_laser = reaction_list[n_reac].n_laser; // Numéro du laser
// cout << " laser numéro " << numero_laser << endl;

    Internal_state Internal_state_in = Mol[numero_mol]; //  état interne de la molecule initialement
    Internal_state Internal_state_out = reaction_list[n_reac].final_internal_state; //  état interne de la molecule après la réaction
    Vecteur3D k_photon_transfer_to_particle; // momentum (wave vector) of the photon transfert to the atoms
    Vecteur3D k_laser = reaction_list[n_reac].k_eff_laser; // it can be zero in the case of lattice for instance (no momentum transfer)
    // Is   Vecteur3D(k,0.,0.) in spontaneous emission and laser_wavector for absorption or emission.


    if (numero_laser != spon_emission) // Absorption or stimulated emission
    {
        if (Internal_state_out.Energy_cm > Internal_state_in.Energy_cm) // TODO (dc#8#): Could be a problem for real diagonalization but here all upper states are above the lower ones so this order is respected
            k_photon_transfer_to_particle = k_laser;   // For emission
        else
            k_photon_transfer_to_particle = -k_laser; // For absorption (low --> up) photon = photon laser but for stimulated emission (up --> low) k_photon = - k_photon laser
    }
    else  // Il y a de l'emission spontanée
    {
        int delta_M = (Internal_state_out.two_M - Internal_state_in.two_M)/2 ;
        complex<double> e_pol_dipole_transition[3] = {reaction_list[n_reac].pol_vector[0], reaction_list[n_reac].pol_vector[1], reaction_list[n_reac].pol_vector[2]}  ;
        Vecteur3D  quantization_axis = reaction_list[n_reac].quant_axis;
        Vecteur3D k_photon_unit_vector = get_unit_vector_spontaneous_emission(r, e_pol_dipole_transition, quantization_axis, delta_M, params); // Gives a random unit vector (in the lab frame) for the spontaneous emission
        // double p_trans = 100.*hPlanck*(Internal_state_in.Energy_cm - Internal_state_out.Energy_cm);   // impulsion de recul (i.e. celle du photon pour l'absorption, - celle pour l'émission)  // Impulsion p = hbar*k = 100 h sigma(cm^-1)
        k_photon_transfer_to_particle =   k_laser.mag()*k_photon_unit_vector;
    }
    Mol[numero_mol] = Internal_state_out; //  c'est l'état interne de la molecule  après réaction
    Vecteur3D delta_v = HBAR*k_photon_transfer_to_particle / Mol[numero_mol].get_mass(); // Emission spontanée aléatoire, p = hbar k
    Mol[numero_mol].inc_vel(delta_v); //  impulsion p = \hbar k = m v ;

    // MODIFY THE STATUS OF THE PARTICLES
    // May be also need to change the name ?
    if (Mol[numero_mol].exc ==  photoionized)
    {
        Mol[numero_mol].set_mass (Mol[numero_mol].get_mass() - ME);
        Mol[numero_mol].set_charge(-QE);
    }
    if (Mol[numero_mol].exc ==  photodetached)
    {
        Mol[numero_mol].set_mass(Mol[numero_mol].get_mass() - ME);
        Mol[numero_mol].set_charge(0.);
    }
    if (Mol[numero_mol].exc ==  annihilation)
    {
// Mol[numero_mol].set_mass(0.); // We keep the same mass for now
        Mol[numero_mol].set_charge(0.); // Change only the charge --> No problem for Lorentz force
    }

    number_photons ++;

    return numero_mol;
}

// Gives a random unit vector (in the lab frame) for the spontaneous emission for a transition delta_M=-1,0,+1 and a quantization axis
// Or (if diagonalized) for a polarization vector e_pol:
//  e_pol[q] = normalized dipole transition <i|d_(q)|j>; q=-1,0,1. So in the quantification axis e_pol = sum_q epol_q E^q
// the probability distribution linked with f(r)= (3/8π)[1-|r.e_pol|^2]
Vecteur3D get_unit_vector_spontaneous_emission(const gsl_rng * r,  complex<double> e_pol_dipole_transition[3], Vecteur3D quantization_axis, int delta_M, FitParams &params)
{
    Vecteur3D point;
    if ( (params.LocateParam("is_Levels_Lines_Diagonalized")->val) )
        /** There is  diagonalization so  vector e_pol


        So to produce the desired spontaneous emission patern for the random variable r. We use the
        Acceptance-rejection sampling method (Von Neumann):
        https://fr.wikipedia.org/wiki/Méthode_de_rejet (The wikipedia algorithm is better in the French version)


        such that r(θ, φ) has the probability distribution linked with f(r)= (3/8π)[1-|r.e_pol|^2]

                But because the differential element is dΩ = d(cos θ)dφ.


        We first choose randomly an r in an isotropy distribution:
        (Isotropy means the density is proportional to solid angle
         (so cos θ is uniform u=(2u1 − 1) and φ is uniform (2π u2) for u1 and u2 uniform in [0,1]).

        1) we choose u randomly in [-1,1] et φ uniform in [0,2π].

        So with the probability distribution g(u,φ)=1/(4 π) X_[-1,1](u)   X_[0,2 π](φ) (that is uniform =1/(4pi))

        f(u,φ) = 3/(8 π) [1-|r(u,φ).e_p|^2] =<  3/(8 π) X_[-1,1](u)   X_[0,2 π](φ) =  3/(8 π) * (4 π) g(u,φ) = 1.5 g(u,φ)

        2) we choose v randomly in [0,1] and keep the r(theta is given by the angle phi, theta) found if v < f(u,φ) / (1.5 g(u,φ)). If not we go back to 1.

        We use 3/(8 π) F= f and 1/(4 π) G = g so v < f(u,φ) / (1.5 g(u,φ)) = 3/(8 π) F/(1.5 * 1/(4 π) G ) = F/G

        **/
    {
        double cos_theta,phi,v; // random variables u = cos_theta
        double sin_theta;
        double norm2 =  norm(e_pol_dipole_transition[0])+ norm(e_pol_dipole_transition[1])+ norm(e_pol_dipole_transition[2]);// I use the C++ notation so norm[z]=|z|^2;
        complex<double> e_pol[3] = {e_pol_dipole_transition[0]/sqrt(norm2), e_pol_dipole_transition[1]/sqrt(norm2), e_pol_dipole_transition[2]/sqrt(norm2)} ;  // e_pol[q] = normalized dipole transition <i|d^(q)|j>; q=-1,0,1
        // e_pol = sum_q epol_q E^q in the quantification frame
        double F; // probability distribution
        do
        {
            cos_theta = gsl_ran_flat (r, -1., 1.); // cos (theta)
            phi= gsl_ran_flat (r, 0, 2.*pi);
            v = gsl_ran_flat (r, 0, 1);
            sin_theta = sqrt(1.-cos_theta*cos_theta);
            // f= is ( 3. /(8.*pi) )* F
            // F=   1. - pow((sin_theta)*(e_pol(2)-e_pol(0))*cos(phi)/sqrt(2.)+e_pol(1)*cos_theta,2) -  pow((sin_theta)*(e_pol(0)+e_pol(2))*sin(phi)/sqrt(2.),2) ; // F=1-|r(uφ,),.e_p|^2
            const complex<double> i(0., 1.);
            F=   1. - norm( (sin_theta)*cos(phi)*(e_pol[0]-e_pol[2])/sqrt(2.)+ e_pol[1]*cos_theta +  i*sin_theta*sin(phi)*(e_pol[0]+e_pol[2])/sqrt(2.) ) ; // F=1-|r(uφ,),.e_p|^2. recall norm[z]=|z|^2
            // I think the previous expression is wrong (contravariant vs covariant)
        }
        while (v > F);
        point= Vecteur3D(sin_theta*cos(phi),sin_theta*sin(phi),cos_theta);// Vecteur unitaire du vecteur d'émission spontanée dans le repère de l'axe de quantification.
    }
    else /** There is no diagonalization so
        the states are Pure and transition are between given M states so the calcul is analytical

            Selon D. Steck "Classical and Modern Optics" 17.70  ou "Quantum and Atom Optics" 1.47 ou De façon plus générale les formules sont données par Eq. (7.436)
            angular-distribution function for linear and circular polarization are
            fˆz(θ, φ) = 3/8π sin^2(θ) et f±(θ, φ) = 3/16π  (1 + cos^2(θ))
            De facon générale la formule est donnée par
            f_(e_p)(r)= (3/8π)[1-|r.e_p|^2]
            avec e_p le vecteur de polarisation e_p = z pour linéaire (pi) et e_(+/-) =  ∓(ˆx ± iˆy)/√2 en circulaire

            On retrouve aussi ( the probability density distribution for the projection of spontaneous emission along one axis
            (3/8p_recoil)[1 + p^2/p_recoil^2])

            Pauline Yzombard made  the calculation (cf Master thesis)! ON utilise ses résultats ici

            **/
    {
        double u = gsl_ran_flat (r, -1., 1.);  // tirage au hazard distribution uniforme sur [-1,1]; // nombre aléatoire pour ennsuite calculé la distribution pour theta par F^-1(u)
        double theta=0.; // Angle de l'émission spontanée

        if (delta_M !=0) //Pour transition sigma+/sigma-. Densité de proba :3/(16 PI)* 1+cos²(théta)
        {
            double D1=2.*u+sqrt(1.+4.*u*u);
            double g1=(-1.+pow(D1,2./3.))/(1.*pow(D1,1./3.));
            theta=acos(g1);
        }
        else   // transition Pi.  densité de proba :3/(8 PI)sin²(théta)
        {
            complex<double> D2;
            complex<double> g2;
            D2 =complex<double>(u,0.)+sqrt(complex<double>(-1.+u*u,0.));
            g2=(-1.*(complex<double>(1.,0.)+pow(D2,2./3.))+complex<double>(0,sqrt(3.))*((complex<double>(1.,0.)-pow(D2,2./3.))))/(2.*pow(D2,1./3.));
            theta=acos(g2.real());
        }

        double phi = gsl_ran_flat (r, 0., 2.*pi);
        point= Vecteur3D(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));// Vecteur unitaire du vecteur d'émission spontanée dans le repère de l'axe de quantification.
    }

    Vecteur3D  Angles_quantization_labo;
    Angles_quantization_labo = Euler_angles(quantization_axis); // We rotate the axis. because the quantization frame is not the lab frame

    return rotation_axis_lab(point, Angles_quantization_labo.x(), Angles_quantization_labo.y(), Angles_quantization_labo.z()); // k_spon_unit_vector
}


