/**
Program for interpolating data points on a two dimensional regular grid.
http://en.wikipedia.org/wiki/Tricubic_interpolation

 coming from
Andrea Demetrio <andrea.demetrio@cern.ch>
and commented by Daniel Comparat using the commented code by
 https://github.com/deepzot/likely/tree/master/likely
 and
 https://github.com/nbigaouette/libtricubic



 The headers required for the tricubic interpolation are
 "TricubicInterpolation_module.h" and "LekienMarsden_matrix.h".

 You have to instantiate an object of class LM_Interpolator and give it the required input values.


**/

#ifndef TRICUBIC_LM_INTERPOLATION
#define TRICUBIC_LM_INTERPOLATION

#include<iostream>
#include<cmath>
#include<stdlib.h>


/**
Map x,y,z to a point dx,dy,dz in the cube [0,n1) x [0,n2) x [0,n3)
**/

class LM_Interpolator
{
public:
    LM_Interpolator();
    LM_Interpolator(int NewXDimension, int NewYDimension, int NewZDimension, double* Xvalues, double* Yvalues, double* Zvalues);
    LM_Interpolator(int NewXDimension, int NewYDimension, int NewZDimension, double delta_x, double delta_y, double delta_z, double X0 = 0, double Y0 = 0, double Z0 = 0);
    void Init(int NewXDimension, int NewYDimension, int NewZDimension, double* Xvalues, double* Yvalues, double* Zvalues);
    void Init(int NewXDimension, int NewYDimension, int NewZDimension, double delta_x, double delta_y, double delta_z, double X0 = 0, double Y0 = 0, double Z0 = 0);
    ~LM_Interpolator();
    bool GetValue(double*** matrix, const double &x, const double &y, const double &z, double &ans, double &dervx, double &dervy, double &dervz);
    bool EvaluateIndices(const double &x, const double &y, const double &z);
    double GetValue(double*** matrix, const double &x, const double &y, const double &z);
    const char* GetInterpolationType();
protected:
    int ijk_to_index(int i, int j, int k);
    void PrepareLMIndicesArray();
    void GetLMIndices(int index, int& i, int& j);
    void GetIndex(const double &x, const double &y, const double &z);
    void Prepare_Coefficients(double*** matrix);
    double PartialXDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ);
    double PartialYDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ);
    double PartialZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ);
    double PartialXYDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ);
    double PartialXZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ);
    double PartialYZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ);
    double PartialXYZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ);
private:
    double *m_coefficients;
    double *m_xvalues, *m_yvalues, *m_zvalues;
    int m_xdimension, m_ydimension, m_zdimension;
    int m_xindex, m_yindex, m_zindex;
    double m_delta_x, m_delta_y, m_delta_z;
    int *m_nonzero_indices;
    int m_nonzero_counter;
};

#endif
