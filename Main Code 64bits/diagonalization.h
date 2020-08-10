#include <Eigen/Eigen>
#include  "Internal_state.h"
#include  "molecule.h"
#include  "params.h"


using namespace Eigen;
using namespace std;

// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and  update all Level[n].Energy_cm
// Need to be improve with the proper v axis
void Diagonalization_Energy(vector <Internal_state> &Level, double B, double v, SelfAdjointEigenSolver<MatrixXd> &es);

// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and dipoles (in Debye) update all Level[n].Energy_cm
void Diagonalization_Energy_dipole(vector <Internal_state> &Level, double B, double v,  SelfAdjointEigenSolver<MatrixXd> &es,  MatrixXd d[]);

// diagonalized the Hamiltionian for the current molecule its field etc.. and give the eigenvectors and eigenvalues and dipoles (in Debye) update all Level[n].Energy_cm
void Diagonalization(vector <Internal_state> &Level, const Molecule &my_mol, const Field &fieldB, const Field &fieldE,
                                      FitParams &params,  SelfAdjointEigenSolver<MatrixXd> &es,  MatrixXd d[]);
