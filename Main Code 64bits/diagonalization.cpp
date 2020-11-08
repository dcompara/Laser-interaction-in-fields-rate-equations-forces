#include "diagonalization.h"
#include "laser.h"                          // for the Euler_angles

#include  <iostream>                       // to include cout, cin





// diagonalized the Hamiltionian for the current molecule its field etc.. and give the eigenvectors and eigenvalues and dipoles (in Debye) update all Level[n].Energy_cm.
void Diagonalization(vector <Internal_state> &Level, const Molecule &my_mol, const Field &fieldB, const Field &fieldE,
                     FitParams &params,  SelfAdjointEigenSolver<MatrixXcd> &es, MatrixXcd &H, MatrixXcd &E0_cm, MatrixXcd &Zeeman_cm_B, MatrixXd d0[], MatrixXcd d[])
{
    const int nb_levels=32; //   Level.size();
    // I tried a dynamical size (or vector) but may be not enough and was not easily compatible with the matrix and speed. But should be tried again

    const int nb_incoherent_low_levels=5; //  These are the "dead" levels number 0, 1, ... nb_incoherent_low_levels-1  for annihilation so for incoherent dipole sum
    const int nb_incoherent_high_levels=7; //  These are the "continuum" levels number nb_levels-nb_incoherent_high_levels, nb_levels-2, nb_levels-1  for photoionization so for incoherent dipole sum

    /*** Hamiltonian matrix element: Zeeman + Stark (External field + Dynamical Stark) ***/

    Vecteur3D r,v,B,F,dip_ij;
    r = my_mol.get_pos();
    v = my_mol.get_vel();
    B = fieldB.get_Field(r); // I neglect the dynamical modification of B field for v << c cf https://en.wikipedia.org/wiki/Classical_electromagnetism_and_special_relativity
    F = fieldE.get_Field(r) + v.cross(B); // External field + Dynamical Stark, noted F to avoid confusion with E = Energy
    complex<double> d_dot_F;

    // Hamiltonian. We assume here an adiabatic following of the state and so that their quantization axis is the B axis
    for (int i=0; i<nb_levels; i++)
        for (int j=0; j<nb_levels; j++)
        {
            dip_ij = Vecteur3D(d0[0](i,j),d0[1](i,j),d0[2](i,j))*Debye/(100.*hPlanck*C); // (to have the dipole such that d.F in cm^-1)
            dipole_dot_Electric_Ffield(dip_ij, B, F, d_dot_F); // Hamiltonian Stark = - d.F = +e r.F
            H(i,j) = E0_cm(i,j) + B.mag()*Zeeman_cm_B(i,j) + d_dot_F; // H_ij = 0_<i | H | j>_0 = 0_<j | H | i>_0  (cf Eq (3) of Dermer PRA 40, 5526 (1989)
//           cout << " i,j " << i <<  " " << j << "   " << E0_cm[i][j] << " " << Zeeman_cm_B[i][j] << "  "  << Stark_cm_Bv[i][j] << endl;
        }

    /*** diagonalization: gives new eigen energies (in Level[n].Energy_cm) and eigen vectors  ***/
// TODO (Daniel#1#): For speed. We can think of computing the Eigenvectors and Eigenvalues using the small dt step and not recalculating all the matrix. The cost of the computation is about $ 9n^3 $ if the eigenvectors are required and $ 4n^3/3 $ if they are not required.
   es.compute(H);  // calculate the new eigenvectors |i> (i start from 0)  and new eingen_Energies
    // E_i = es.eigenvalues()(i) (in increasing order)
    // es.eigenvectors()(i,j) = 0<i | j>  gives (i=line, j = column index) the new (column) vector |j> in function of the old |i>_0

//
//    for( int j0 = 0; j0 < nb_levels; j0++ ) // Ligne
//    {
//        for( int j = 0; j < nb_levels; j++ )
//        {
//            // cout <<  es.eigenvectors()(j0,j).real() << " "; // 0<j0 | j>
//
//            cout <<  (round(100.* es.eigenvectors()(j0,j).real()))/100  << "+i" << (round(100.* es.eigenvectors()(j0,j).imag()))/100 << "   ";
//        }
//        cout << endl;
//    }
//    cout << endl << endl;





    for( int n = 0; n < nb_levels; n++ )
    {
        Level[n].Energy_cm = es.eigenvalues()(n);
    }


    /**** calcul of the new dipoles (in Debye) in d[q]

    evec = es.eigenvectors() verifie evec(j0,j) = 0<j0 | j>  gives (j0=line, j = column index) the new (column) vector |j> in function of the old |j0>_0
    d0[q+1]_i0 j0 = 0_<i0 | d_q | j0>_0 . So d[q+1]_ij = <i | d_q | j> = Sum i0,j0   <i |i0>0 0<i0| d_q | j0>0 0<j0|j> =  Sum i0,j0  evec^dag (i,i0)  d_q(i0,j0)  evec(j0,j)

    IN CONCLUSION: the new dipole are given by d[polar] = evec^dag.d0[polar].evec =  <i | d_q | j> with evec_j = |j> = sum_|j>_0   0_<j| j>.


    But for some transition (low to dead levels and high to continuum we use an incoherent sum. To do this we use the Hadamar product A°B = (A.cwiseProduct(B)) defined by (A°B)__ij = A_ij B_ij
    d_incoh[q+1]_ij^2 = <i | d_incoh_q | j> = Sum i0,j0   | <i |i0>0 0<i0| d_q | j0>0 0<j0|j> |^2 =  Sum i0,j0  |evec^dag (i,i0)  d_q(i0,j0)  evec(j0,j)|^2 = Sum i0,j0  evec^dag (i,i0) evec^dag(i,i0)^*   d_q(i0,j0) d_q(i0,j0)^*  evec(j0,j) evec(j0,j)^*
    d_incoh[q+1]_ij^2 = Sum i0,j0  evec^dag°evec^dag* (i,i0)    d_q°d_q* (i0,j0)  evec°evec* (j0,j)  =  (evec^dag°evec^dag* .   d_q°d_q* .  evec°evec*)_ij
    with for such state typically the |j> levels does not change so |j>=|j0>_0 and the last sum is Sum i0  |<i |i0>0 0<i0| d_q | j=j0>0 |^2
    But we keep in case of order changing in the levels files ...

    IN CONCLUSION dSQUARE_incoh = evec^dag°evec^dag* .   d_q°d_q* .  evec°evec*
    ***/

    for(int n_polar = -1; n_polar <= 1; n_polar++)
    {
        // We first calculate all dipoles as incoherent sum and then we replace the central part with the proper coherent dipoles. Remark: We do twice the work (so it is slower, but what is long is probably the diagonalization) but is is very simple to write

        // incoherent sum for all dipoles
        d[n_polar+1] = ( ( es.eigenvectors().adjoint()).cwiseProduct(es.eigenvectors().adjoint().conjugate()) ) * ( d0[n_polar+1].cwiseProduct(d0[n_polar+1].conjugate()) ) * ( (es.eigenvectors()).cwiseProduct(es.eigenvectors().conjugate()) );
        // this is not yet the matrix of the dipole but the one of dipole squared d_incoh[q+1]_ij^2
        for (int i=0; i<nb_levels; i++)
            for (int j=0; j<nb_levels; j++)
                d[n_polar+1](i,j) = sqrt(d[n_polar+1](i,j)); // To go back to d_incoh[q+1]_ij

        // We replace the central part with the proper block matrix of coherent dipole. Using matrix.block<p,q>(i,j) = Block of size (p,q), starting at (i,j)
        const int i = nb_incoherent_low_levels;
        const int p = nb_levels - nb_incoherent_high_levels - nb_incoherent_low_levels;  // size of coherent block
        d[n_polar+1].block<p,p>(i,i) = ( (es.eigenvectors().adjoint())*d0[n_polar+1]*(es.eigenvectors()) ).block<p,p>(i,i);
    }
}


//   d.E (for the Stark effect that is -d.E)
// where the dipole  vector d= sum_p d_p e^p is given in the local quantification axis (typically the magnetic field)
void dipole_dot_Electric_Ffield (const Vecteur3D& dipole, const Vecteur3D& axe_quant,  const Vecteur3D& E, complex<double> &d_dot_E)
{
    double dp,d0,dm;
    dm = dipole(0);
    d0 = dipole(1);
    dp = dipole(2);

    Vecteur3D Euler_angles_axe_quant = Euler_angles(axe_quant);
    Vecteur3D Euler_angles_axe_electric_field = Euler_angles(E);

// The link between the polar angles (theta, phi) and the Euler angle (alpha, beta, gamma) in ZXZ convention as we used them (Wikipedia)  are
    // alpha= phi +pi/2; beta = theta; gamma = -pi/2
    double phi_F = Euler_angles_axe_quant(0) - pi/2.; // polar angle for the quantization axis
    double theta_F = Euler_angles_axe_quant(1); // polar angle for the quantization axis
    double phi_k = Euler_angles_axe_electric_field(0) - pi/2.; // polar angle for the electric field axis (along the vector E)
    double theta_k = Euler_angles_axe_electric_field(1); // polar angle for the electric field axis (along the vector E)


    complex<double> twice_dip; // Twice the dipole
    const complex<double> i(0., 1.);

    /*** We try to optimize the calcul so to reduce at maximmum the numebr of operations ****/

    double sin_theta_k = sin(theta_k); // The notation k is used here because this is the notation for a field as if it was a laser propagating along k with a polarization along k by a laser (cf effectif_dipole_local() function)
    double sin_theta_F = sin(theta_F);
    double cos_theta_k = cos(theta_k);
    double cos_theta_F = cos(theta_F);
    double cos_phi_F_minus_phi_k = cos(phi_F-phi_k);
    double sin_phi_F_minus_phi_k = sin(phi_F-phi_k);

    double sqrt2 = sqrt(2.);

    double dp_minus_dm = dp-dm;
    double dp_plus_dm = dp+dm;

    twice_dip  = cos_theta_k*(2.*d0*cos_theta_F + sqrt2*dp_minus_dm*sin_theta_F)+
                 sin_theta_k*(-sqrt2*dp_minus_dm*cos_theta_F*cos_phi_F_minus_phi_k+
                              2.*d0*cos_phi_F_minus_phi_k*sin_theta_F - i*sqrt2*dp_plus_dm*sin_phi_F_minus_phi_k);

    d_dot_E = E.mag() * twice_dip/2.;
}

