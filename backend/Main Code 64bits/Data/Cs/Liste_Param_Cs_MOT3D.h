/************************* File datacard ***************************

Can be used in the program by creating a variable such as
double Name_Parameter = params.LocateParam("Name_Parameter")->val

# 0 = false ; 1 =  true

***********************************************************************/
########## 	Particles  	###############################
#
// Nb of different particles. For now only the first one is laser cooled (Level and Lines correspond to it)
@Nb_type_of_Mol	1

// Le choix du nom (BaF, Cs2, NH, Cs, CO, Li6Cs, Li7Cs, Rb85Cs, Rb87Cs, Ps, C2minus, Ps_minus,P_bar) ne donne que la masse mais pas le nom des fichiers qu'il faut changer ensuite"
#1st type of particle
// so Mol[0] to Mol[Nom_Mol[0]-1]
@Nom_Mol[0]	Cs
// It is the number of molecules that are laser cooled.
@N_Mol[0]  100
@Temp_ini_x[0] 10e-6
@Temp_ini_y[0] 10e-6
@Temp_ini_z[0] 10e-6


// Choix en position: taille fixe (sigma_pos) ou à partir de la densité
//  -1 taille fixe et ordonne les positions au départ (selon un axe mettre les autres axes aléatoires)
//  0 taille fixe donnée par size (gaussien)
//  1 un pot. magn. linéaire --> Laplace
//  2 un pot. magn. quadratique --> Gaussien
#
@Procedure_init_x[0]   0
@Procedure_init_y[0]   0
@Procedure_init_z[0]   0
// Taille (x,y,z) si on choisit taille fixe
@size_x[0]  1000e-6
@size_y[0]	1000e-6
@size_z[0]  1000e-6
// Initial position added to the random one
@offset_x[0]	0.
@offset_y[0]	0.
@offset_z[0]	0.
// Initial velocity added to the random one
@v0_x[0]	0.
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
@size_x[1]	1e-3
@size_y[1]	1e-3
@size_z[1]	1e-3
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
// 0.1 m/(q B) for Boris (10^-8 at 0.0001T for 3me mass; 2 e-8 for C2- 1Telsa). A good test is to remove lasers and check Energy conservation

//for t< t_scaling_max
@dt_dyn_epsilon_param  1e-8
//for t> t_scaling_max
// fin du temps.
//@t_fin  10e-9
@t_fin  1000e-6
// time interval between diagnostics (in cout) output
@dt_dia 1e-5
//@dt_dia 60e-9
// time interval between output of snapshots (draw particles)
@dt_out 1e-6
//@dt_out 100e-6
#
###################### GRAPHICS and OUTPUT ###############################
#
@SIZE_affichage	500e-6
// Temps d'attende entre 2 affichages. Permet de ne pas avoir un affichage trop rapide
@t_wait_affichage   1e-10
// 0 = false ; 1 =  true
@Graphics 1
// rotate 90 degree along the vector (rot_axe_x, rot_axe_y, rot_axe_z)
// if no rotation we would have x (red arrow) toward the right, y (green arrow) up and z (blue arrow) toward the screen.
// DEFAULT IS (1,0,0) = x right, y  in screen, z up (gravity down)
// (0,1,0) = x out screen, y  up, z left
// (0,0,1) = x up, y  left, z out of screen
@rot_axe_x  1
@rot_axe_y  0
@rot_axe_z  0



##### VELOCITY  SCALING####

@velocity_scaling 0
// velocity scaling ON (1) or OFF (0)
@time_max_vel_scaling 3e-3
// temps pendant lequel on procède au time scaling. (temps max)
// Time for step in velocity scalling
@dt_scal 1e-7
@coupling_efficiency 0.01
// in percentage (max 1) 0.1=10%

// coupling parameter of the BERENDSEN THERMOSTAT Algorithm
// T(t+dt) ~  T(t)+coupling_efficiency(T0-T(t))



#
# OUTPUT
#
// gives the number of the manifold that we do not want to take into account for stats
// (can be dead level or photoionized (-1) one or  ..)
// To take all manifold into account just put a number that is not used such as -10
@num_manifold_not_studied   -10
// numéro du niveau étudié pour faire des stats.
// -1 pour toutes les molécules
@num_niveau_etudie  -1
#
######### 	CHAMPS EXTERNES SI units (T, T/m, T/m^2 ...)	################################
#
// We cannot for now have both electric and magnetic field easily (just because we have only one parameter in the Level file
// Thus we have to choose if this parameter is Zeeman (0) or Stark shift (10)
@type_of_field_for_internal_state_shift 0
// But both fields will be used of the Lorentz force between charged particles
// For instance in Penning trap --> There is a magnetic and electric. But the magnetic is fully treated (Zeeman + Lorentz force). But the electric is only for the Lorentz FOrce not for the Stark effect.
// Or for Paul trap with (not implemented) or without micro-motion.

@type_field_read_E    0
@type_field_read_B    0
// THIS is only for the "field_for_internal_state_shift" the other one will be by default at zero so in 2nd +nth order
// DEFAULT
// 0: 2nd order plus a nth order

// Helmotlz coils
// 1:  Field in Helmoltz coils (so usualy goes with @type_of_default_field 0)

// File grids + Electric 2nd order: EQUIPARTITIONED POINTS (at least per axis, the spacing can be different for each axis). Ordered by column (1st increasing then second then third etc ..)
// 2: Field map 3D from 2D cylindrical symmetry: 4 columns r,z, F_r(r,z); F_z(r,z)
// 3: Field map 3D from 2D cylindrical symmetry: F_r(r,z); F_z(r,z)+ derivative d/dr; d/dz and d^2/drdz
// 4: Field map 3D: 6x,y,z, Bx, By, Bz (TO BE DONE)

// For speed. WE SUGGEST TO CALCULATE THE DEIVATIVES IF TYPE 2 USING Field::Calculate_Derivative_Matrix

#
## MAGNETIC FIELD ##
#
#
#
@Nb_bobines 0
// Si oui (1) les bobines sont doublées (on an ajoute une décalée de +r)
@is_Helmholtz 1
// Ecart entre les bobines
@gap_bobines    15e-3
 // Courant dans les bobines --> Champ au centre µ0 I/ (2 r)
@courant_bobines  0
// rayon des bobines
@rayon_bobines  3.e-3
#
// Champ magn selon x,y et z. se décompose par composante: Example selon Ox: B_x + grad_B_x x + grad_grad_B_x x^2 + Bn x^n
@B_x	0.
@B_y	0.000000001
@B_z	0.
@grad_B_x	-0.05
@grad_B_y	-0.05
@grad_B_z	0.1
@grad_grad_B_x	0.
@grad_grad_B_y	0.
@grad_grad_B_z	0.
@n_value_B    3
@Bn_x	0.
@Bn_y	0.
@Bn_z	0.
#
## ELECTRIC FIELD ##
#
// Electric Field along x,y and z. example Ox: E_x + grad_E_x x + grad_grad_E_x x^2 + En x^n
@E_x	0.
@E_y	0.
@E_z	0.
@grad_E_x	0
@grad_E_y	0
// V(z) the voltage being only 0.1meV (10K) at -8cm (turning point) so - 0.0001 (z/8cm)^12
// E(z) = -12*0.0001/(0.08)^12 * z^11 ~1.7 e10 z^11
@grad_grad_E_x	0.
@grad_grad_E_y	0.
@grad_grad_E_z	0.
@n_value_E    3
@En_x	0.
@En_y	0.
@En_z	0.
#
######### 	LASERS 	########################################
#
// Parametre multiplicatif de la puissance des lasers
@scale_Power  1
// Paramètre additif de la fréquence de tous les lasers
// Si Offset_Detuning_cm est >0 le laser est plus bleu (*1K detunning*)
@Offset_Detuning_cm  0
// Parametre multiplivatif de la largeur spectrale laser
@scale_Gamma 1

// Nb de laser utilisés (pas forcément le nombre ci-après que peut être plus grand)
@Nb_laser 7
// Switch time laser. we switch between
//laser 0, 1, 2, .. Nb_laser/2-1 to
//laser Nb_laser/2, Nb_laser/2+1 ..., Nb_laser-1 after dt_switch_1 and back after dt_switch_2 and so on ans so forth
// 0 = false ; 1 =  true
@Is_Laser_Switched 0
@dt_switch_1 3e-8
@dt_switch_2 1e-8

# Premier laser. Laser n°1 (called number 0 in the C++ program)
@waist_pos_x[0]	0.
@waist_pos_y[0]	0.
@waist_pos_z[0]	0.
@direction_x[0]	0.
@direction_y[0]	0.
@direction_z[0]	1.
// @waist_x[0]	5e-3
// @waist_y[0]	5e-3
// Mettre si on veux un seul waist
@waist[0]	1e-2
// Energie_cm[1]	0 field 33.356005422904 29841.11126
// 0.1T: 29874.4672665725 , 33.3984724, 33.362644, 33.33569101
// 29840.606 NH
// 11630.258 BaF
// 15935.8455177462 LiCs
// 13709.648988 RbCs
// 11732.1813268101 Cs
// PS :
// T=150K , Gamma =2.5e5 MHz DeltaE=4.63
//Etransition sans decalage(B=0T)=41148.41996
//Etransition sans décalage (B=1T)=41154.453
// Pour Cs 11732.1813268101
@Energie_cm[0]  11732.1812
@Gamma_L_MHz[0]	1
@Power[0]	0.001
@Pol_circulaire_left_sp[0]    1
@Pol_circulaire_right_sm[0]     0
//  façonné ---> Energie_cm+2*Gamma > EnergiE_trans_cm > Energie_cm (OBSOLETTE)
 // gaussien = 5, lorentzien = 6. (En fait Gaussien ne marche que si le laser est à résonnance et large spectralement)
@type_laser[0]  6
@nu_offset_MHz[0] 0
@nu_repetition_MHz[0] 80
@nu_individual_comb_line_MHz[0] 80
// Interférence de tous les lasers qui ont le numéro  coherent_avec_laser_num[0]
// si  coherent_avec_laser_num[0] = -1 ce laser est seul et n'interfère avec personne.
// Attention bien mettre dans l'ordre un laser j interfère TOUJOURS avec un laser i<=j. S'il y a interférence entre i et j alors le plus petit i interfère avec i aussi
@coherent_avec_laser_num[0]  -1
// est t'il utilisé pour le (re)pompage (i.e. sa fréquence varie avec scale_temp_pompage? False  par défaut.
// Si oui on ajoute scale*Temp_initial à son énergie
@is_pompage[0] 0
@is_rep[0] 0
// CW = 0, femto =1, multi_mode=2, pulse=3, faconne = 4, gaussien = 5, lorentzien = 6

#Deuxieme laser. Laser n°2
@waist_pos_x[1]	0.
@waist_pos_y[1]	0.
@waist_pos_z[1]	0.
@direction_x[1]	0.
@direction_y[1]	0.
@direction_z[1]	-1.

@waist[1]	1e-2
@Energie_cm[1]  11732.1812
@Gamma_L_MHz[1]	1
@Power[1]	0.001
@Pol_circulaire_left_sp[1]    1
@Pol_circulaire_right_sm[1]    0

@type_laser[1]  6
@nu_offset_MHz[1] 0
@nu_repetition_MHz[1] 80
@nu_individual_comb_line_MHz[1] 80

@coherent_avec_laser_num[1]  -1

@is_pompage[1] 0
@is_rep[1] 0


#  Laser n°3
@waist_pos_x[2]	0.
@waist_pos_y[2]	0.
@waist_pos_z[2]	0.
@direction_x[2]	0.
@direction_y[2]	1.
@direction_z[2]	0.

@waist[2]	1e-2
@Energie_cm[2]  11732.1812
@Gamma_L_MHz[2]	1
@Power[2]	0.001
@Pol_circulaire_left_sp[2]    0
@Pol_circulaire_right_sm[2]    1

@type_laser[2]  6
@nu_offset_MHz[2] 0
@nu_repetition_MHz[2] 80
@nu_individual_comb_line_MHz[2] 80

@coherent_avec_laser_num[2]  -1

@is_pompage[2] 0
@is_rep[2] 0

# Laser n°4
@waist_pos_x[3]	0.
@waist_pos_y[3]	0.
@waist_pos_z[3]	0.
@direction_x[3]	0.
@direction_y[3]	-1.
@direction_z[3]	0

@waist[3]	1e-2
@Energie_cm[3]  11732.1812
@Gamma_L_MHz[3]	1
@Power[3]	0.001
@Pol_circulaire_left_sp[3]    0
@Pol_circulaire_right_sm[3]    1

@type_laser[3]  6
@nu_offset_MHz[3] 0
@nu_repetition_MHz[3] 80
@nu_individual_comb_line_MHz[3] 80

@coherent_avec_laser_num[3]  -1

@is_pompage[3] 0
@is_rep[3] 0


#Cinquième laser 5
@waist_pos_x[4]	0.
@waist_pos_y[4]	0.
@waist_pos_z[4]	0.
@direction_x[4]	-1.
@direction_y[4]	0.
@direction_z[4]	0.

@waist[4]	1e-2
@Energie_cm[4]  11732.1812
@Gamma_L_MHz[4]	1
@Power[4]	0.001
@Pol_circulaire_left_sp[4]    0
@Pol_circulaire_right_sm[4]    1

@type_laser[4]  6
@nu_offset_MHz[4] 0
@nu_repetition_MHz[4] 80
@nu_individual_comb_line_MHz[4] 80

@coherent_avec_laser_num[4]  -1

@is_pompage[4] 0
@is_rep[4] 0

#Sixieme Laser : n°6
@waist_pos_x[5]	0.
@waist_pos_y[5]	0.
@waist_pos_z[5]	0.
@direction_x[5]	1.
@direction_y[5]	0.
@direction_z[5]	0.

@waist[5]	1e-2
@Energie_cm[5]  11732.1812
@Gamma_L_MHz[5]	1
@Power[5]	0.001
@Pol_circulaire_left_sp[5]    0
@Pol_circulaire_right_sm[5]    1

@type_laser[5]  6
@nu_offset_MHz[5] 0
@nu_repetition_MHz[5] 80
@nu_individual_comb_line_MHz[5] 80

@coherent_avec_laser_num[5]  -1

@is_pompage[5] 0
@is_rep[5] 0

#Sixieme Laser : n°7 (repumping)
@waist_pos_x[6]	0.
@waist_pos_y[6]	0.
@waist_pos_z[6]	0.
@direction_x[6]	0.
@direction_y[6]	0.
@direction_z[6]	1.

@waist[6]	1e-2
// Cs repiumping 11732.479585
@Energie_cm[6]  11732.479585
@Gamma_L_MHz[6]	1
@Power[6]	0.001
@Pol_circulaire_left_sp[6]    0
@Pol_circulaire_right_sm[6]    1

@type_laser[6]  6
@nu_offset_MHz[6] 0
@nu_repetition_MHz[6] 80
@nu_individual_comb_line_MHz[6] 80

@coherent_avec_laser_num[6]  -1

@is_pompage[6] 0
@is_rep[6] 0
#
#
######### 	PARAMETRES POMPAGE OPTIQUE + SISYPHE 	####################
#
//
@num_niveau_first   1
 // numéro du niveau de moindre pente pour le Sisyphe
@num_niveau_second	6
// numéro du niveau excité
@num_niveau_exc	3
#
// Repompage ou pompage forcé oui ou non (1 true or 0 false)
@repompage_force    0
@Pompage_optique_force  0
#
@rate_repompe   1e7
@E0_coupure_repompe_cm 33.37
// scale en temp (pompage Sisyphe à scale_temp * k_B T)
@scale_temp_pompage 3
// scale en temp (repompage à scale_temp_rep * k_B T)
@scale_temp_rep 0.4
// Nb de répétition (0 = aucune) du processus Sisyphe. A t_repet on recommence, on garde les molécules mais on remet à t=0 les paramètres
@nb_repet    0
@t_repet    1000000
#
#
##########   PARAMETRES DE SORTIE FICHIER ####################
#
// sortie de la datacard?
@is_DataCard_out 2
// # 0 = false ; 1 =  true ; 2 = 2 file (one with only datacard, one with data but without datacard)
@is_param_scanned_out 1
#
#
######### 	FICHIERS FC 	################################
#
// Lecture du ficher Franck-COndon. Si oui cela crée un autre fichier "nom_file_Levels" et "nom_file_Lines" utilisant l'existant
// (supposé être entre v''=0 et v'=0) et multiplié par le fichier F.C. et avec les énergies et Bv donnés
// Si is_File_FC=true on crée ces fichiers avec noms  + "new_v_2Jmax"
@is_File_FC false
// Fichier contenant les facteurs FC entre vA (lignes)  et vX (colonnes)
@nom_file_FC_vAvX	Data/BaF/FC_vA_0_15Lignes_vX_0_14Colonnes.dat
// Fichier contenant les facteurs niveaux vibrationels, les énergies et les Bv (en cm^-1)
@nom_file_E_vA  Data/Ps/Ps_Levels_0B_photo.dat
@nom_file_E_vX  Data/Ps/Ps_Lines_0B_photo.dat
// Nb of vibrational levels in X state; 0,1,...,NXmax-1 in the FC file
@NXmax  15
// Nb of vibrational levels in A state; 0,1,...,NAmax-1
@NAmax  16
// Nb of vibrational levels used to calculate the new file for X state; 0,1,...,NX-1
// Used also to give (for the J,v output) the NX used in the current file if @is_File_FC false)
// This is also used for the drawing of the molecules if the graphic is on
@NX_out 20
// Nb of vibrational levels used to calculate A; 0,1,...,NA-1
@NA_out 20
// Max +1 of rotational levels à mettre dans le fichier de sortie. Si on veux 2J=0,1,2 mettre 3.
// Used also to give (for the J,v output) the 2JX_max used in the current file if @is_File_FC false)
// This is also used for the drawing of the molecules if the graphic is on
@N_Two_JX_out_max 20;
@N_Two_JA_out_max 20;
#
// true pour fixer la durée de vie à duree_vie sinon calculées à partir des forces de transition # 0 = false ; 1 =  true
// EN FAIT EST UTILSE JUSTE pour les largeur, car la durée de vie des transitions vient du dipole dans  rates_molecule_spon
@if_fixed_lifetime	0
// en force la durée de vie de l'état excité d'avoir une durée de vie de 15 ns
@duree_vie	15e9
#
######### 	NOM_FICHIERS 	################################
#
// Fichier contenant les Niveaux (état, énergie, ...)

@nom_file_Levels            Data/Cs/Cs_levels.dat

//@nom_file_Levels    Data/C2anion/c2anion_levels_2Jmax10_final_min.dat
//@nom_file_Levels    Data/C2anion/LEVEL_v0_v1_v2_v3_v4_minimal.dat
// @nom_file_Levels    Data/C2anion/LEVEL_v0_minimal.dat
@nom_file_Lines             Data/Cs/Cs_lines.dat
#
// Fichier contenant les spectres laser doit finir en .dat
// Fichier contenant les transitions lasers en fait ce nom +"[numero laser]" où numéro=0 pour le premier laser
@nom_file_Laser_Spectrum    Data/Cs/Laser_Spectrum.dat
#
// If needed File containing maps for Elec or Magn fields
// the type can be:
 // Field map 3D x,y,z: 6 columns,x,y,z, Fx,Fy, Fz
// Field map 3D but from 2D cylindrical symmetry: r,z: 4 columns r,z, Fr, Fz
// Field map 3 from 2D cylindrical symmetry: 10 columns, r,z, Fr, Fz, + derivative d/dr; d/dz and d^2/drdz
// THE FILE SHOULD BE ORDERED by increasing values of the point coordinates
// FOR 2D
// r0 z0 Fr Fz;
// r0 z0+deltaz0 Fr Fz ..
// ......
//  rfinal Zfinal Fr Fz
#
//@nom_file_Magn_Field_3D   Data/C2anion/Bfinal_Positif_10_columns_center_z0_add_r_negatif.dat
@nom_file_Magn_Field_3D     Data/C2anion/Bfinal_Positif_10_columns_center_z0_add_r_negatif.dat
@nom_file_Elec_Field_3D     Data/C2anion/ElectricField2D.dat
@nom_sortie_donnees	        Data/donnee_Mol.dat
@nom_sortie_pulse	        Data/sortie_pulse.dat
@nom_sortie_rate	        Data/sortie_rate.dat
@nom_sortie_donnees_Data	Data/data_card.dat
@nom_sortie_temp	        Data/sortie_temp.dat
@nom_sortie_scal            Data/sortie_scal.dat
#
@nom_fichier_random_gen     Data/random_gen.txt
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
@choix_epsilon  1e-6
// La même valeur permet d'avoir toujours la même séquence. Une valeur négative utilise le fichier pour renouveler la séquence. If 0 the standard seed from the original implementation
@Seed_Init_Random_Number_Generator  7
#
#
#############  listes des VARIABLES SCANNEES ou variées temporellement #############
#
// La liste commence par BEGIN_OF... et finie par END_OF...
// Il faut mettre exactement le même nom précédé de SCAN et suivit des valeurs suivantes (au moins jusqu'a nb_steps) séparée par TAB:
// name_scanned, minvalue, maxvalue, nombrestep, bool_is_scanned, bool_is_time_dependent, tau_var
// Nb_steps =ntervall --> 1 signifie que l'on va prendre 2 valeurs min et max. nb_step=2 on prend 3 valeurs: min, moitié et max  ...
 // val_t0 = minv + steps * (maxv - minv )/nbstep
 // tau =   taux de variation temporelle --> val = val_t0 exp^(-t/tau).
//  Les paramètres seront modifiés en exp^(-t/tau). Si on ne précise pas tau et que c'est time dependent la valeur sera Tau_Modif
//  Si il y a plusieurs valeur on met le numéro du tableau [numéro]
#
# Tau_Modif  est le temps par défaut si aucun autre n'est marqué
// Modification des paramètres au cours du temps si aucune autre valeur n'est précisée.
@Tau_Modif    0.17
#
#
// Pour savoir si on scan de façon aléatoire ou ordonnée
@is_Scan_Random    false
# name  minv    maxv    nbstep  is_scanned  is_time tau
BEGIN_OF_FITPARAMS
@SCAN_Energie_cm[0]  13000 13600   6   false    false
@SCAN_scale_Power  0 4 4   false    false
@SCAN_scale_Detuning_cm -10 0 5 false false
@SCAN_scale_Gamma 0.5 2  3 false false
END_OF_FITPARAMS

