/************************* File datacard ***************************

Can be used in the program by creating a variable such as
double Name_Parameter = params.LocateParam("Name_Parameter")->val

# 0 = false ; 1 =  true

***********************************************************************/
########## 	Particles  	###########################################################################
########## 	all in SI units (if not clearly expressed by the name) 	###############################
#
// Nb of different particles. For now only the first one is laser cooled (Level and Lines correspond to it)
@Nb_type_of_Mol	1

// Le choix du nom (BaF, Cs2, NH, Cs, CO, Li6Cs,Laminus, Li7Cs, Rb85Cs, Rb87Cs, Ps, C2minus, Ps_minus,P_bar,Osminus) ne donne que la masse mais pas le nom des fichiers qu'il faut changer ensuite"
#1st type of particle
// so Mol[0] to Mol[Nom_Mol[0]-1]
@Nom_Mol[0]	Ps
// It is the number of molecules that are laser cooled.
@N_Mol[0]  1
@Temp_ini_x[0] 1e-10
@Temp_ini_y[0] 1e-10
@Temp_ini_z[0] 1e-10


// Choice in position: fixed size (sigma_pos) or from density
// -1 fixed size and orders the positions at the start (according to an axis put the other axes randomly)
// 0 fixed size given by size (Gaussian)
// 1 a pot. magn. linear --> Laplace (the coefficient comes from F1)
// 2  pot. magn. quadratic for neutral particle --> Gaussian
// 3  pot. magn. quadratic for charged particle --> Gaussian
// 3 for pot. elec. quadratic (electric field linear) but for CHARGED particles (for instance in a Paul trap)--> Gaussian
// 4 perfect ordered gaussian in velocity (random in position)
// 5 effusive beam. Meaning as in 0 but we keep only the positive velocities

#
@Procedure_init_x[0]   0
@Procedure_init_y[0]   0
@Procedure_init_z[0]   0
// Taille (x,y,z) si on choisit taille fixe
@size_x[0]  0.1e-6
@size_y[0]	0.1e-6
@size_z[0]  0.1e-6
// Initial position added to the random one
@offset_x[0]	0
@offset_y[0]	0
@offset_z[0]	0
// Initial velocity added to the random one
@v0_x[0]	1e5
@v0_y[0]	0.
@v0_z[0]    0.
#
#2nd type of particle
#
// so Mol[Nom_Mol[0]] to Mol[Nom_Mol[0]+Nom_Mol[1]-1]
@Nom_Mol[1]	P_bar
@N_Mol[1]  10
@Temp_ini_x[1] 4
@Temp_ini_y[1] 4
@Temp_ini_z[1] 4
@Procedure_init_x[1]   0
@Procedure_init_y[1]   0
@Procedure_init_z[1]   0
@size_x[1]	1e-4
@size_y[1]	1e-4
@size_z[1]	6e-4
@offset_x[1]	0.
@offset_y[1]	0.
@offset_z[1]	0.
@v0_x[1]	0.
@v0_y[1]	0.
@v0_z[1]	0.
#
##########   Temps KMC, paramètres de sortie ####################
#
// For control parameter to determine the dynamical time step size in second  to check convergence
// typical is 0.001*waist/velocity (or 0.001*lambda/velocity for lattices)
// or near 0.001 cyclotron period  2pi m/(q B) for Leapfrog (in magnetic field case)
// 0.1 m/(q B) for Boris (10^-8 at 0.0001T for 3me mass; 2 e-8 for C2- 1Telsa). *
// A good test is to remove lasers and check Energy conservation

//for t< t_scaling_max
@dt_dyn_epsilon_param  1e-10
//for t> t_scaling_max
// fin du temps.
//@t_fin  10e-9
@t_fin  0.1e-9
// time interval between diagnostics (in cout) output
@dt_dia 5e-9
// time interval between output of snapshots (draw particles)
@dt_out 2e-9
#
###################### GRAPHICS and OUTPUT ###############################
#
@SIZE_affichage	7e-3
// Temps d'attende entre 2 affichages. Permet de ne pas avoir un affichage trop rapide
@t_wait_affichage   1e-1
// 0 = false ; 1 =  true
@Graphics 0
// rotate 90 degree along the vector (rot_axe_x, rot_axe_y, rot_axe_z)
// if no rotation we would have x (red arrow) toward the right, y (green arrow) up and z (blue arrow) toward the screen.
// DEFAULT IS (1,0,0) = x right, y  in screen, z up (gravity down)
// (0,1,0) = x out screen, y  up, z left
// (0,0,1) = x up, y  left, z out of screen
@rot_axe_x  1
@rot_axe_y  0
@rot_axe_z  0

// Scaling of  the drawing of the molecules for the vibration v_max and rotation J_max
@v_maX 8
@two_J_maX 10;


#
# OUTPUT
#
// gives the number of the manifold that we do not want to take into account for stats
// (can be dead level or photoionized (-1) one or  ..)
// To take all manifold into account just put a number that is not used such as -10
@num_manifold_not_studied   -3
// numéro du niveau étudié pour faire des stats.
// Be carefule: c'est le numéro (ordre dans la lecture du fichier en commençant par 0) du niveau pas de la manifold
// -1 pour toutes les molécules
@num_niveau_etudie  -1
#
######### 	Diagonalization	################################
#
// Are the energy levels diagonalized or simple calculed using the simple analytical formula used in the code (linear, quadratic or two level case)
//  0 = false (we use the standard simple formulas) ; 1 =  true (we diagonalized)
@is_Levels_Lines_Diagonalized   1
#
######### 	CHAMPS EXTERNES SI units (T, T/m, T/m^2 ...)	################################
#########   FOR MORE COMPLEX CASE  see is_Levels_Lines_Diagonalized    #########
#
// We cannot for now have both electric and magnetic field easily (just because we have only one parameter in the Level file for the field shift
// Thus we have to choose if this parameter is Zeeman or Stark shift
// But both fields will be used of the Lorentz force between charged particles
// 0: magnetic (the Delta and C parameter for Level file will be a Zeeman shift)
// 1: electric (Stark shift)
@type_of_field_for_internal_state_shift 0


// But both fields are used for charges particles (Lorentz) forces
// For instance in Penning trap --> There is a magnetic and electric. But the magnetic is fully treated (Zeeman + Lorentz force).
// But the electric is only for the Lorentz FOrce not for the Stark effect.

/*** How to define the fields: analtically (0), by coils(1), or by reading grids points on different geometries (2,3,4,..) ****/
// DEFAULT
// 0: 2nd order plus a nth order
// 1:  Field in  (anti-)Helmotlz coils coils (so usualy goes with a field that is magnetic)

// Then we have file grids + Electric 2nd order: EQUIPARTITIONED POINTS (at least per axis, the spacing can be different for each axis). Ordered by column (1st increasing then second then third etc ..)
// 2: Field map 3D from 2D cylindrical symmetry: 4 columns r,z, F_r(r,z); F_z(r,z).   // NOte if needed in the code, and for speed. the derivatiove car be calculated using efficiently Field::Calculate_Derivative_Matrix
// 3: Field map 3D from 2D cylindrical symmetry: F_r(r,z); F_z(r,z)+ derivative d/dr; d/dz and d^2/drdz
// 4: Field map 3D: 6x,y,z, Bx, By, Bz (TO BE DONE)
@type_field_read_E    0
@type_field_read_B    0


#
## MAGNETIC FIELD definition (if) coils ##
#
// Dans le cas 1: Nombre de bobines. Axe z positions  (N-1/2) * z_gap auquel on ajoute r si (anti-)Helmotz
@Nb_bobines 1
// Si  is_Helmholtz oui (1) les bobines sont doublées (on an ajoute une décalée de +r)
// Si -1 c'est anti-is_Helmholtz
// si autre (comme 0)  on n'ajoute pas de bobine on ne crée que les bobines prévues
@is_Helmholtz -1
// Ecart entre les bobines
@gap_bobines    25e-3
 // Courant dans les bobines --> Champ au centre µ0 I/ (2 r)  --> 0.63 mT for 1A and r=1mm
@courant_bobines  1
// rayon des bobines
@rayon_bobines  25e-3
#
## MAGNETIC FIELD definition analytical ##
#
// Champ magn selon x,y et z. se décompose par composante: Example selon Ox: B_x + grad_B_x x + grad_grad_B_x x^2 + Bn x^n
// NEVER put 0 but something small like 1e-10 to keep a quantization axis
@B_x	0.
@B_y	0.
@B_z	1e-10
@grad_B_x	0.
@grad_B_y	0.
@grad_B_z   0.
@grad_grad_B_x	0.
@grad_grad_B_y	0.
@grad_grad_B_z	0.
@n_value_B    3
@Bn_x	0.
@Bn_y	0.
@Bn_z	0
#
## ELECTRIC FIELD ##
// Electric Field along x,y and z. example Ox: E_x + grad_E_x x + grad_grad_E_x x^2 + En x^n
// NEVER put 0 but something small like 1e-10 to keep a quantization axis (if given by E)
@E_x	0.
@E_y	0.
@E_z	0.
@grad_E_x	0.
@grad_E_y	0.
@grad_E_z	0.
@grad_grad_E_x	0.
@grad_grad_E_y	0.
@grad_grad_E_z	0.
@n_value_E    3
@En_x	0.
@En_y	0.
@En_z	0
#
######### 	LASERS 	########################################
#
// Parametre multiplicatif de la puissance des tous les lasers
@scale_Power 1

// Paramètre additif de la fréquence de tous les lasers
@Offset_Detuning_cm  0

// Parametre multiplivatif de la largeur spectrale laser
@scale_Gamma 1

// Nb de laser utilisés (pas forcément le nombre ci-après qui peut être plus grand mais ne seront pas utilisés)
@Nb_laser 0

// We can swhich between lasers
@Is_Laser_Switched 0
// IF t is between 0 and T1 (modulo T1+T2) THEN the lasers On are the one with numbers between 0 and N_laser/2
// IF t is between T1, T1+T2 (modulo T1+T2) THEN the lasers On are between    N_laser/2+1 and N_laser.
@dt_switch_1 3e-8
@dt_switch_2 1e-8

# Premier laser. Laser n°1 (called number 0 in the C++ program)
@waist_pos_x[0]	0.
@waist_pos_y[0]	0.
@waist_pos_z[0]	0
@direction_x[0]	0.
@direction_y[0]	1.
@direction_z[0]	0.

// waist in SI units
@waist[0]	5e-3
@Energie_cm[0] 41148.3848

@Gamma_L_MHz[0] 1e-3
@Power[0]	200

// Vector laser polarization (in the laser propagation frame) See User Guide
// Polarization can be purely circular (sigma+ or sigma -). Example: sigma + --> Pol_circulaire_left_sp = 1 and @Pol_circulaire_right_sm =-1
// Can also be linear example eX = eX=(e-1-e+1)/sqrt(2). SoPol_circulaire_left_sp = -0.7071 and Pol_circulaire_right_sm[ = 0.7071;
// Then the angle_psi_degree is (for linear polarization) the angle (so 90° if we want eY polarisation)
// For linear polarization at 54.7356 degree it is  sp= -0.707107  sm= 0.707107 and angle 54.7356
// by the way this creates 1/3 sigma+, 1/3 sigma- and 1/3 pi polarization (for a Y laser beam and quantization axes along z)
@Pol_circulaire_left_sp[0]    1
@Pol_circulaire_right_sm[0]   0.
@polar_angle_degree[0]  0.

// The laser types are numerous (see laser.h file) and can be use on purpose (cf  intensity_Convolution_linewidth function) for field_ionization, black body type, femtosecond
//  spon_emission = -1,    CW = 0, femto =1,  multi_mode=2,  pulse=3,  faconne = 4 (Obsolete because of laser_spectrum now),  gaussien = 5,  lorentzien = 6,   comb = 7,  pseudo_BBR =8,  field_ionization =9
// the more standard is  gaussien = 5
// spectral shape and shaping can be done using Laser_spectrum file
@type_laser[0]  5

// If femtosecond + comb laser (type_laser = 7) we havee these parameters. If not they are useless
@nu_offset_MHz[0] 0
@nu_repetition_MHz[0] 80
@nu_individual_comb_line_MHz[0] 80

// Interférence de tous les lasers qui ont le numéro  coherent_avec_laser_num[0]
// si  coherent_avec_laser_num[0] = -1 ce laser est seul et n'interfère avec personne.
// Attention bien mettre dans l'ordre un laser j interfère TOUJOURS avec un laser i<=j. S'il y a interférence entre i et j alors le plus petit i interfère avec i aussi
@coherent_avec_laser_num[0]  -1

#Deuxieme laser. Laser n°2
@waist_pos_x[1]	0.
@waist_pos_y[1]	0.
@waist_pos_z[1]	0
@direction_x[1]	0.
@direction_y[1]	0.
@direction_z[1]	-1.

@waist[1]	10e-3
@Energie_cm[1] 41148.3848

@Gamma_L_MHz[1]	1e5
@Power[1]	1
@Pol_circulaire_left_sp[1]    0.7071067812
@Pol_circulaire_right_sm[1]   -0.7071067812
@polar_angle_degree[1]  45

@type_laser[1]  5

@coherent_avec_laser_num[1]  -1

#
##########   PARAMETRES DE SORTIE FICHIER ####################
#
// sortie de la datacard?
@is_DataCard_out 2
// # 0 = false ; 1 =  true ; 2 = 2 file (one with only datacard, one with data but without datacard)
@is_param_scanned_out 1
#


#
######### 	NOM_FICHIERS 	################################
#
// Fichier contenant les Niveaux (état, énergie, ...)

@nom_file_Levels    Data/Ps/Ps_Levels_23.dat
@nom_file_Lines     Data/Ps/Ps_Lines_23.dat

// Fichier contenant les spectres laser doit finir en .dat
// Les dichiers contenant les transitions lasers seront donc ce nom +"[numero laser]" où numéro=0 pour le premier laser
@nom_file_Laser_Spectrum    Data/Ps/Laser_Spectrum.dat
#
@nom_sortie_donnees	Data/donnee_Mol.dat
@nom_sortie_rate	Data/sortie_rate.dat
@nom_sortie_donnees_Data	Data/data_card.dat
#
@nom_fichier_random_gen Data/random_gen.txt
#
#
#############	Choix des algorithmes Monte Carlo,  N-corps; aléatoire 	 ##############
#
// Vérifier lesquels marchent avant -1, 0, are O.K. !!!
// Aucun_MC = -1 (does not even calculate the rates), Kinetic_Monte_Carlo = 0, Random_Selection_Method = 1, First_Reaction_Method = 2,
// Fast_Rough_Method = 3 is a new one evolving every ~0.1/rate_max time so typically N_Mol/10 evolve during one time step
@Choix_algorithme_Monte_Carlo	0

// Aucun_N_corps = -1 (mais photon recoil), Verlet_acc = 1 (sans force dipolaire), Verlet_pot (avec potentiel dipolaire) = 2,
// si = 6(c'est Verlet_pot_gradient_high_order avec potentiel dipolaire et calcul du gradient dans one_body à l'ordre supérieur)
// Yoshida6_acc = 3,       // 6 th order symplectic algorithm based on Verlet_acc
//    Yoshida6_pot = 4,        // 6 th order symplectic algorithm based on Verlet_pot
// Boris_Buneman  (with Magnetic field for charged particles) = 7
@Choix_algorithme_N_corps	1

// Choix du epsilon en position (en metre) pour calculer la dérivée du potentiel.
// Une variation de 1/100 du potentiel selon epsilon semble bien
// So 1e-6 for standard lasers or 1e-9 when interferences are present
// Ne pas hésiter à tester aussi epsilon <0. Et a vérifier avec le epsilon_param
@choix_epsilon  1e-8

// La même valeur permet d'avoir toujours la même séquence.
// This is EXTREMELLY useful for debuging because the program is always the same.
// We can use a negative value (but I never did it) to use a file to  have a new sequence a each run.
@Seed_Init_Random_Number_Generator  2
#
#
#############  lists of SCAN or time-varying VARIABLE VARIABLES #############

// The list starts with BEGIN_OF... and ends with END_OF...
// You have to put exactly the same name preceded by SCAN and followed by the following values (at least until nb_steps) separated by TAB:
// name_scanned, minvalue, maxvalue, nombrestep, bool_is_scanned, bool_is_time_dependent, tau_var
// Nb_steps = nb_intervall --> so 1 means that we will take 2 values min and max. nb_step=2 we take 3 values: min, half and max ...
 // val_t0 = minv + steps * (maxv - minv )/nbstep
 // tau = time rate of change --> val = val_t0 exp^(-t/tau).
// If true is given then the parameters will thus be changed in exp^(-t/tau). If we don't specify tau and it's time dependent the value will be Tau_Modif
// If there are several values, the table number is set [number].
#
# Tau_Modif is the default time if no other is marked
// Modification of parameters over time if no other value is specified.
@Tau_Modif 1e-3
#
#
// To find out if we scan randomly or orderly
@is_Scan_Random    false
# name  minv    maxv    nbstep  is_scanned  is_time tau
BEGIN_OF_FITPARAMS
@SCAN_scale_Power  5 40 4   false    false
@SCAN_Offset_Detuning_cm -2.5 -1 5 false false
@SCAN_scale_Gamma 0.3 0.6 3   false    false
@SCAN_Tau_Modif 0.5e-3 2e-3 2 false false
@SCAN_B_z   0   5  1000  false false
END_OF_FITPARAMS

