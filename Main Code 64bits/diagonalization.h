#include <Eigen/Eigen>
#include  "Internal_state.h"
#include  "molecule.h"
#include  "params.h"

#include <complex>                      // Pour les calculs de d_dot_F


using namespace Eigen;
using namespace std;


// diagonalized the Hamiltionian for the current molecule its field etc.. and give the eigenvectors and eigenvalues and dipoles (in Debye) update all Level[n].Energy_cm
void Diagonalization(vector <Internal_state> &Level, const Molecule &my_mol, const Field &fieldB, const Field &fieldE,
                                      FitParams &params,  SelfAdjointEigenSolver<MatrixXcd> &es,  MatrixXcd d[]);


//   d.E (for the Stark effect that is -d.E)
// where the dipole  vector d= sum_p d_p e^p is given in the local quantification axis (typically the magnetic field)
void dipole_dot_Electric_Ffield (const Vecteur3D& dipole, const Vecteur3D& axe_quant,  const Vecteur3D& E, complex<double> &d_dot_E);
