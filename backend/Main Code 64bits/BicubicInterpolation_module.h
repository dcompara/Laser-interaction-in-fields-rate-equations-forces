<<<<<<< HEAD
/**
Program for interpolating data points on a two-dimensional regular grid.
http://en.wikipedia.org/wiki/Bicubic_interpolation

Originally written by:
Andrea Demetrio <andrea.demetrio@cern.ch>

Modified by:
Daniel Comparat (changes include replacing `double** matrix` with `std::vector<std::vector<double>> &matrix`).
Replaced dynamic allocation with `std::vector`.

In this program, the axes are referred to as x and y, but they are used as r and z.

The original program was based on a matrix `m[x][y]` (i.e., `matrix[m_xindex][m_yindex]`).
It was modified to work with `matrix_field[j]`, where `j` denotes a column corresponding to a given field (e.g., Br, Bz, Bx, By, or others).
Each column contains the values from `m[x][y]` using the index `m_index * m_ydimension + m_yindex`.
Thus, the old `matrix[m_x][m_y]` is transformed into `matrix_field[j][m_x * m_ydimension + m_y]`.

Additional changes:
- Declared some variables as `const`.
- Removed `m_xindex`, `m_yindex` from being class members because they are non-constant and more suitable as local variables.
- Local variables like `m_coefficients`, `m_nonzero_indices`, and `m_nonzero_counter` are no longer class members, making it easier to have functions marked `const` (even when these parameters are modified internally).
- All arrays were replaced with `std::vector` to avoid issues caused by frequent construction and destruction of dynamically allocated arrays.
**/


#ifndef BICUBIC_INTERPOLATION
#define BICUBIC_INTERPOLATION

#include<iostream>
#include<cmath>
#include<stdlib.h>
#include <vector>

using namespace std;

class Bicubic_Interpolator
{
public:
    Bicubic_Interpolator() {};

    Bicubic_Interpolator(int NewXDimension, int NewYDimension, double delta_x, double delta_y, double X0 = 0, double Y0 = 0);   // Create equidistant points with offset X0 -and Y0) and spacinf delta_x (or delta y).  Initialize the value of the matrix by the Xvalues[i]=X0 + i*delta_x and same for

    void Init(int NewXDimension, int NewYDimension, double delta_x, double delta_y, double X0 = 0, double Y0 = 0);


// Cannot be const because called other function like Get_index that modify the m_index
    bool GetValue(const vector<vector<double> > &matrix, const double &x, const double &y, double &ans, double &dervx, double &dervy, int column_number=0, int type_field_read = 0) const; // Gives  ans=f(x,y), df/dx(x,y), df/dy(x,y)
    double GetValue(const vector<vector<double> > &matrix, const double &x, const double &y, int column_number=0, int type_field_read = 0) const; //Return evaluation of f(x,y)

// same but faster if only f(x,y) is needed
    double GetValueInterpolated(const vector<vector<double> > &matrix, const double &x, const double &y, int column_number) const;

    bool EvaluateIndices( int &m_xindex, int &m_yindex) const;       //Return false if  indices OUT OF RANGE and thus put BOTH them at zero
    void ShowIndex(const double &x, const double &y) const;   //Print index of the point of the grid close to the X,Y point
    const char* GetInterpolationType(); //Get interpolator name

    /* PARTIAL DISCRETE DERIVATIVES using f'(x) = (f(x+eps)-f(x-eps) )/2 eps between 2 grid points */
    // same as bellow but with matrix
    double XDerivative(const vector<vector<double> > &matrix, int line_number, int column_number=0) const;
    double YDerivative(const vector<vector<double> > &matrix, int line_number, int column_number=0) const;
    double XYDerivative(const vector<vector<double> > &matrix, int line_number, int column_number=0) const;

protected:
    int ij_to_index(int i, int j) const; //Transform matrix index into vector index
    void PrepareMatrixIndicesArray(); // //Prepare nonzero matrix index from the 16*16 BicubicMatrix
    void GetMatrixIndices(int index, int& i, int& j) const;
    // void GetIndex(const double &x, const double &y);
    void GetIndex(const double &x, const double &y, int &m_xindex, int &m_yindex) const; // return the index near x,y, or (0,0) if one is out of the grid range

    void Prepare_Coefficients(const vector<vector<double> > &matrix,  int m_xindex, int m_yindex,   vector<double> &m_coefficients, int column_number=0, int type_field_read = 0) const ;

    /* PARTIAL DISCRETE DERIVATIVES using f'(x) = (f(x+eps)-f(x-eps) )/2 eps between 2 grid points */
    double PartialXDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int column_number=0) const;
    double PartialYDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int column_number=0) const;
    double PartialXYDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int column_number=0) const;



private: // m_ stands for Matrix_
    vector <double> m_xvalues, m_yvalues; // Values
    int m_xdimension, m_ydimension; // size of the Box (in number of indices)
    // int m_xindex, m_yindex;   // Index of the point in the box
    double m_delta_x, m_delta_y;  // spacing ion the box (Assume that points are equipartitioned)
    vector <int> m_nonzero_indices; // use only non zero of the matrix to be faster
    int m_nonzero_counter;
};

#endif
=======
/**
Program for interpolating data points on a two-dimensional regular grid.
http://en.wikipedia.org/wiki/Bicubic_interpolation

Originally written by:
Andrea Demetrio <andrea.demetrio@cern.ch>

Modified by:
Daniel Comparat (changes include replacing `double** matrix` with `std::vector<std::vector<double>> &matrix`).
Replaced dynamic allocation with `std::vector`.

In this program, the axes are referred to as x and y, but they are used as r and z.

The original program was based on a matrix `m[x][y]` (i.e., `matrix[m_xindex][m_yindex]`).
It was modified to work with `matrix_field[j]`, where `j` denotes a column corresponding to a given field (e.g., Br, Bz, Bx, By, or others).
Each column contains the values from `m[x][y]` using the index `m_index * m_ydimension + m_yindex`.
Thus, the old `matrix[m_x][m_y]` is transformed into `matrix_field[j][m_x * m_ydimension + m_y]`.

Additional changes:
- Declared some variables as `const`.
- Removed `m_xindex`, `m_yindex` from being class members because they are non-constant and more suitable as local variables.
- Local variables like `m_coefficients`, `m_nonzero_indices`, and `m_nonzero_counter` are no longer class members, making it easier to have functions marked `const` (even when these parameters are modified internally).
- All arrays were replaced with `std::vector` to avoid issues caused by frequent construction and destruction of dynamically allocated arrays.
**/


#ifndef BICUBIC_INTERPOLATION
#define BICUBIC_INTERPOLATION

#include<iostream>
#include<cmath>
#include<stdlib.h>
#include <vector>

using namespace std;

class Bicubic_Interpolator
{
public:
    Bicubic_Interpolator() {};

    Bicubic_Interpolator(int NewXDimension, int NewYDimension, double delta_x, double delta_y, double X0 = 0, double Y0 = 0);   // Create equidistant points with offset X0 -and Y0) and spacinf delta_x (or delta y).  Initialize the value of the matrix by the Xvalues[i]=X0 + i*delta_x and same for

    void Init(int NewXDimension, int NewYDimension, double delta_x, double delta_y, double X0 = 0, double Y0 = 0);


// Cannot be const because called other function like Get_index that modify the m_index
    bool GetValue(const vector<vector<double> > &matrix, const double &x, const double &y, double &ans, double &dervx, double &dervy, int column_number=0, int type_field_read = 0) const; // Gives  ans=f(x,y), df/dx(x,y), df/dy(x,y)
    double GetValue(const vector<vector<double> > &matrix, const double &x, const double &y, int column_number=0, int type_field_read = 0) const; //Return evaluation of f(x,y)

// same but faster if only f(x,y) is needed
    double GetValueInterpolated(const vector<vector<double> > &matrix, const double &x, const double &y, int column_number) const;

    bool EvaluateIndices( int &m_xindex, int &m_yindex) const;       //Return false if  indices OUT OF RANGE and thus put BOTH them at zero
    void ShowIndex(const double &x, const double &y) const;   //Print index of the point of the grid close to the X,Y point
    const char* GetInterpolationType(); //Get interpolator name

    /* PARTIAL DISCRETE DERIVATIVES using f'(x) = (f(x+eps)-f(x-eps) )/2 eps between 2 grid points */
    // same as bellow but with matrix
    double XDerivative(const vector<vector<double> > &matrix, int line_number, int column_number=0) const;
    double YDerivative(const vector<vector<double> > &matrix, int line_number, int column_number=0) const;
    double XYDerivative(const vector<vector<double> > &matrix, int line_number, int column_number=0) const;

protected:
    int ij_to_index(int i, int j) const; //Transform matrix index into vector index
    void PrepareMatrixIndicesArray(); // //Prepare nonzero matrix index from the 16*16 BicubicMatrix
    void GetMatrixIndices(int index, int& i, int& j) const;
    // void GetIndex(const double &x, const double &y);
    void GetIndex(const double &x, const double &y, int &m_xindex, int &m_yindex) const; // return the index near x,y, or (0,0) if one is out of the grid range

    void Prepare_Coefficients(const vector<vector<double> > &matrix,  int m_xindex, int m_yindex,   vector<double> &m_coefficients, int column_number=0, int type_field_read = 0) const ;

    /* PARTIAL DISCRETE DERIVATIVES using f'(x) = (f(x+eps)-f(x-eps) )/2 eps between 2 grid points */
    double PartialXDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int column_number=0) const;
    double PartialYDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int column_number=0) const;
    double PartialXYDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int column_number=0) const;



private: // m_ stands for Matrix_
    vector <double> m_xvalues, m_yvalues; // Values
    int m_xdimension, m_ydimension; // size of the Box (in number of indices)
    // int m_xindex, m_yindex;   // Index of the point in the box
    double m_delta_x, m_delta_y;  // spacing ion the box (Assume that points are equipartitioned)
    vector <int> m_nonzero_indices; // use only non zero of the matrix to be faster
    int m_nonzero_counter;
};

#endif
>>>>>>> 8a5eff06c2dff7ab9c152968cca57c83445d78f7
