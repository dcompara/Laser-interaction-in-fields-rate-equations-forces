/************************* Fichier datacard ***************************

Ce n'est pas un fichier C++ le .h est l� juste pour avoir un bon �diteur

On les utilise dans le programme par params.LocateParam("Nom_Parametre")->val

# 0 = false ; 1 =  true


***********************************************************************/
########## 	MOLECULES  	###############################
#
#
// Le choix du nom (BaF, Cs2, NH, Cs, CO, Li6Cs, Li7Cs, Rb85Cs, Rb87Cs, Ps, C2moins) ne donne que la masse mais pas le nom des fichiers qu'il faut changer ensuite
@Nom_Mol	Cs
@N_Mol  500
@Temp_ini 16e-3
#
#
##########   Temps KMC,affichage, param�tres de sortie ####################
#
// For control parameter to determine the dynamical time step size in second  to check convergence
// typical is 0.001*waist/velocity (or 0.001*lambda/velocity for lattices)
@dt_dyn_epsilon_param  1e-7
// fin du temps.
@t_fin  350e-6
// time interval between diagnostics (in cout) output
@dt_dia 349e-6
// time interval between output of snapshots (draw particles)
@dt_out 1e-6
#
@SIZE_affichage	3e-3
// Temps d'attende entre 2 affichages. Permet de ne pas avoir un affichage trop rapide
@t_wait_affichage   1e-10
// 0 = false ; 1 =  true
@Graphics 1
#
#
######### 	CHAMPS EXTERNES SI units (T, T/m, T/m^2 ...)	################################
#
// Type de champ
// 0: donn� au 2�me ordre
// 1: N-bobines d'axes Ox
// 2: carte de champ (A FAIRE)
@type_field	0
#
// Dans le cas 1: Nombre de bobines
@Nb_bobines 0
// Si oui (1) les bobines sont doubl�es (on an ajoute une d�cal�e de +r)
@is_Helmholtz 1
// Ecart entre les bobines
@gap_bobines    15e-3
 // Courant dans les bobines --> Champ au centre �0 I/ (2 r)
@courant_bobines  0
// rayon des bobines
@rayon_bobines  3.e-3
#
// Champ magn selon x,y et z. se d�compose par composante: Example selon Ox: B_x + grad_B_x x + grad_grad_B_x x^2
@B_x	0.
@B_y	0.
@B_z	0.001
@grad_B_x	0.
@grad_B_y	0.
@grad_B_z	0.
@grad_grad_B_x	0.
@grad_grad_B_y	0.
@grad_grad_B_z	0.
#
#
######### 	INITIALISATION SAMPLE	########################
#
// Choix en position: taille fixe (sigma_pos) ou � partir de la densit�
//  -1 taille fixe et ordonne les positions au d�part (selon un axe mettre les autres axes al�atoires)
//  0 taille fixe donn�e par size (gaussien)
//  1 un pot. magn. lin�aire --> Laplace
//  2 un pot. magn. quadratique --> Gaussien
#
@Procedure_init_x   0
@Procedure_init_y   0
@Procedure_init_z   0
// Taille (x,y,z) si on choisit taille fixe
@size_x	1e-3
@size_y	1e-3
@size_z	1e-3
// Initial position added to the random one
@offset_x	-0.05
@offset_y	0.
@offset_z	0.
// Initial velocity added to the random one
@v0_x	250.
@v0_y	0.
@v0_z	0.
#
#
######### 	LASERS 	########################################
#
// Parametre multiplicatif de la puissance des lasers
@scale_Power   1
// Param�tre additif de la fr�quence de laser
// Si scale_detuning est >0 le laser est plus bleu
@scale_Detuning_cm 0
// Parametre multiplivatif de la largeur spectrale laser
@scale_Gamma 1
// Nb de laser utilis�s (pas forc�ment le nombre ci-apr�s que peut �tre plus grand)
@Nb_laser	2
# Premier laser. Laser n�1 (called number 0 in the C++ program)
@waist_pos_x[0]	0.
@waist_pos_y[0]	0.
@waist_pos_z[0]	0.
@direction_x[0]	0.
@direction_y[0]	0.
@direction_z[0]	1.
// @waist_x[0]	5e-3
// @waist_y[0]	5e-3
// Mettre si on veux un seul waist
@waist[0]	1.8e-3
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
//Etransition sans d�calage (B=1T)=41154.453
// Pour Cs 11732.1813268101 = 11732.48796 - 0.3066331899
@Energie_cm[0]  11732.1828268101
@Gamma_L_MHz[0]	1
@Power[0]	0.070
@Pol_circulaire_gauche_sp[0]    1
@Pol_circulaire_droit_sm[0]     0
//  fa�onn� ---> Energie_cm+2*Gamma > EnergiE_trans_cm > Energie_cm (OBSOLETTE)
 // gaussien = 5, lorentzien = 6. (En fait Gaussien ne marche que si le laser est � r�sonnance et large spectralement)
@type_laser[0]  6
// Interf�rence de tous les lasers qui ont le num�ro  coherent_avec_laser_num[0]
// si  coherent_avec_laser_num[0] = -1 ce laser est seul et n'interf�re avec personne.
// Attention bien mettre dans l'ordre un laser j interf�re TOUJOURS avec un laser i<=j. S'il y a interf�rence entre i et j alors le plus petit i interf�re avec i aussi
@coherent_avec_laser_num[0]  0
// est t'il utilis� pour le (re)pompage (i.e. sa fr�quence varie avec scale_temp_pompage? False  par d�faut.
// Si oui on ajoute scale*Temp_initial � son �nergie
@is_pompage[0] 0
@is_rep[0] 0
// CW = 0, femto =1, multi_mode=2, pulse=3, faconne = 4, gaussien = 5, lorentzien = 6

#Deuxieme laser. Laser n�2
@waist_pos_x[1]	0.
@waist_pos_y[1]	0.
@waist_pos_z[1]	0.
@direction_x[1]	0.
@direction_y[1]	0.
@direction_z[1]	-1.

@waist[1]	1.8e-3
@Energie_cm[1]  11732.1828268101
@Gamma_L_MHz[1]	1
@Power[1]	0.070
@Pol_circulaire_gauche_sp[1]    0
@Pol_circulaire_droit_sm[1]    1

@type_laser[1]  6

@coherent_avec_laser_num[1]  0

@is_pompage[1] 0
@is_rep[1] 0


#  Laser n�3
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
@Pol_circulaire_gauche_sp[2]    0
@Pol_circulaire_droit_sm[2]    1

@type_laser[2]  6

@coherent_avec_laser_num[2]  -1

@is_pompage[2] 0
@is_rep[2] 0

# Laser n�4
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
@Pol_circulaire_gauche_sp[3]    0
@Pol_circulaire_droit_sm[3]    1

@type_laser[3]  6

@coherent_avec_laser_num[3]  -1

@is_pompage[3] 0
@is_rep[3] 0


#Cinqui�me laser 5
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
@Pol_circulaire_gauche_sp[4]    0
@Pol_circulaire_droit_sm[4]    1

@type_laser[4]  6

@coherent_avec_laser_num[4]  -1

@is_pompage[4] 0
@is_rep[4] 0

#Sixieme Laser : n�6
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
@Pol_circulaire_gauche_sp[5]    0
@Pol_circulaire_droit_sm[5]    1

@type_laser[5]  6

@coherent_avec_laser_num[5]  -1

@is_pompage[5] 0
@is_rep[5] 0

#Sixieme Laser : n�7 (repumping)
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
@Pol_circulaire_gauche_sp[6]    0
@Pol_circulaire_droit_sm[6]    1

@type_laser[6]  6

@coherent_avec_laser_num[6]  -1

@is_pompage[6] 0
@is_rep[6] 0
#
######### 	PARAMETRES POMPAGE OPTIQUE + SISYPHE 	####################
#
// num�ro du niveau �tudi� pour faire des stats.
// -1 pour toutes les mol�cules
@num_niveau_etudie	-1
 // num�ro du niveau de moindre pente pour le Sisyphe
@num_niveau_second	6
// num�ro du niveau excit�
@num_niveau_exc	3
#
// Repompage ou pompage forc� oui ou non (1 true or 0 false)
@repompage_force    0
@Pompage_optique_force  0
#
@rate_repompe   1e7
@E0_coupure_repompe_cm 33.37
// scale en temp (pompage Sisyphe � scale_temp * k_B T)
@scale_temp_pompage 3
// scale en temp (repompage � scale_temp_rep * k_B T)
@scale_temp_rep 0.4
// Nb de r�p�tition (0 = aucune) du processus Sisyphe. A t_repet on recommence, on garde les mol�cules mais on remet � t=0 les param�tres
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
// Lecture du ficher Franck-COndon. Si oui cela cr�e un autre fichier "nom_file_Levels" et "nom_file_Lines" utilisant l'existant
// (suppos� �tre entre v''=0 et v'=0) et multipli� par le fichier F.C. et avec les �nergies et Bv donn�s
// Si is_File_FC=true on cr�e ces fichiers avec noms  + "new_v_2Jmax"
@is_File_FC false
// Fichier contenant les facteurs FC entre vA (lignes)  et vX (colonnes)
@nom_file_FC_vAvX	Data/BaF/FC_vA_0_15Lignes_vX_0_14Colonnes.dat
// Fichier contenant les facteurs niveaux vibrationels, les �nergies et les Bv (en cm^-1)
@nom_file_E_vA  Data/BaF/A_v_E_Bv_th_v_0_15.dat
@nom_file_E_vX  Data/BaF/X_v_E_Bv_v_0_15.dat
// Nb of vibrational levels in X state; 0,1,...,NXmax-1 in the FC file
@NXmax  15
// Nb of vibrational levels in A state; 0,1,...,NAmax-1
@NAmax  16
// Nb of vibrational levels used to calculate the new file for X state; 0,1,...,NX-1
// Used also to give (for the J,v output) the NX used in the current file if @is_File_FC false)
// This is also used for the drawing of the molecules if the graphic is on
@NX_out 4
// Nb of vibrational levels used to calculate A; 0,1,...,NA-1
@NA_out 2
// Max +1 of rotational levels � mettre dans le fichier de sortie. Si on veux 2J=0,1,2 mettre 3.
// Used also to give (for the J,v output) the 2JX_max used in the current file if @is_File_FC false)
// This is also used for the drawing of the molecules if the graphic is on
@N_Two_JX_out_max 10;
@N_Two_JA_out_max 10;
#
// true pour fixer la dur�e de vie � duree_vie sinon calcul�es � partir des forces de transition # 0 = false ; 1 =  true
// EN FAIT EST UTILSE JUSTE pour les largeur, car la dur�e de vie des transitions vient du dipole dans  rates_molecule_spon
@if_fixed_lifetime	0
// en force la dur�e de vie de l'�tat excit� d'avoir une dur�e de vie de 15 ns
@duree_vie	15e9
#
######### 	NOM_FICHIERS 	################################
#
// Fichier contenant les Niveaux (�tat, �nergie, ...)
@nom_file_Levels    Data/Cs/Cs_two_levels.dat
@nom_file_Lines     Data/Cs/Cs_lines.dat
#
// Fichier contenant les spectres laser doit finir en .dat
// Fichier contenant les transitions lasers en fait ce nom +"[numero laser]" o� num�ro=0 pour le premier laser
@nom_file_Laser_Spectrum    Data/Ps/Laser_Spectrum.dat
#
@nom_sortie_donnees	Data/donnee_Mol.dat
@nom_sortie_pulse	Data/sortie_pulse.dat
@nom_sortie_rate	Data/sortie_rate.dat
@nom_sortie_donnees_Data	Data/data_card.dat
@nom_sortie_temp	Data/sortie_temp.dat
#
@nom_fichier_random_gen Data/random_gen.txt
#
#
#############	Choix des algorithmes Monte Carlo,  N-corps; al�atoire 	 ##############
#
// V�rifier lesquels marchent avant!!!
// Aucun_MC = -1, Kinetic_Monte_Carlo = 0, Random_Selection_Method = 1, First_Reaction_Method = 2
@Choix_algorithme_Monte_Carlo	0
// Aucun_N_corps = -1 (mais photon recoil), LeapFrog_Verlet_acc = 1 (sans force dipolaire), LeapFrog_Verlet_pot (avec potentiel dipolaire) = 2,
// si = 6(c'est LeapFrog_Verlet_pot_gradient_high_order avec potentiel dipolaire et calcul du gradient dans one_body � l'ordre sup�rieur)
@Choix_algorithme_N_corps	2
// Choix du epsilon en position (en metre) pour calculer la d�riv�e du potentiel.
// Une variation de 1/100 du potentiel selon epsilon semble bien
// So 1e-6 for standard lasers or 1e-9 when interferences are present
// Ne pas h�siter � tester aussi epsilon <0. Et a v�rifier avec le epsilon_param
@choix_epsilon  1e-8
// La m�me valeur permet d'avoir toujours la m�me s�quence. Une valeur n�gative utilise le fichier pour renouveler la s�quence. If 0 the standard seed from the original implementation
@Seed_Init_Random_Number_Generator  4
#
#
#############  listes des VARIABLES SCANNEES ou vari�es temporellement #############
#
// La liste commence par BEGIN_OF... et finie par END_OF...
// Il faut mettre exactement le m�me nom pr�c�d� de SCAN et suivit des valeurs suivantes (au moins jusqu'a nb_steps) s�par�e par TAB:
// name_scanned, minvalue, maxvalue, nombrestep, bool_is_scanned, bool_is_time_dependent, tau_var
// Nb_steps =ntervall --> 1 signifie que l'on va prendre 2 valeurs min et max. nb_step=2 on prend 3 valeurs: min, moiti� et max  ...
 // val_t0 = minv + steps * (maxv - minv )/nbstep
 // tau =   taux de variation temporelle --> val = val_t0 exp^(-t/tau).
//  Les param�tres seront modifi�s en exp^(-t/tau). Si on ne pr�cise pas tau et que c'est time dependent la valeur sera Tau_Modif
//  Si il y a plusieurs valeur on met le num�ro du tableau [num�ro]
#
# Tau_Modif  est le temps par d�faut si aucun autre n'est marqu�
// Modification des param�tres au cours du temps si aucune autre valeur n'est pr�cis�e.
@Tau_Modif    0.17
#
#
// Pour savoir si on scan de fa�on al�atoire ou ordonn�e
@is_Scan_Random    false
# name  minv    maxv    nbstep  is_scanned  is_time tau
BEGIN_OF_FITPARAMS
@SCAN_Energie_cm[0]  13000 13600   6   false    false
@SCAN_scale_Power  0 4 4   false    false
@SCAN_scale_Detuning_cm -10 0 5 false false
@SCAN_scale_Gamma 0.5 2  3 false false
END_OF_FITPARAMS

