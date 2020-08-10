#include <Eigen/Eigen>
#include  "Internal_state.h"

using namespace Eigen;
using namespace std;

// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and  update all Level[n].Energy_cm
void Diagonalization_Energy(vector <Internal_state> &Level, double B, double v, SelfAdjointEigenSolver<MatrixXd> &es);

// diagonalized the Hamiltionian for B field and v velocity (orthogonal to B) and give the eigenvectors and eigenvalues and dipoles (in Debye) update all Level[n].Energy_cm
void Diagonalization_Energy_dipole(vector <Internal_state> &Level, double B, double v,  SelfAdjointEigenSolver<MatrixXd> &es,  MatrixXd d[]);

// A partir de la matrice des dipole qui contien en première ligne et dernière colonne les M
int Create_dipole_Lines_from_Matrices(const char *nom_file);
