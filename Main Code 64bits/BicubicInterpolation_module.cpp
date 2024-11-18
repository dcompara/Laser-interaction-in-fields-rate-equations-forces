#include "BicubicInterpolation_module.h"
#include "LekienMarsdenMatrix.h"
#include <iostream>
#include  <fstream>                        // to read the data and put them in files
// #include <iomanip>


using namespace std;


// Create equidistant points with offset X0 -and Y0) and spaced delta_x (or delta y)
// Initialize the value of the matrix by the Xvalues[i]=X0 + i*delta_x and same for Y
void Bicubic_Interpolator::Init(int NewXDimension, int NewYDimension, double delta_x, double delta_y, double X0, double Y0)
{

    m_xdimension = NewXDimension;
    m_ydimension = NewYDimension;
    //Prepare X positions array

    // vector <double> m_xvalues(m_xdimension);   // m_xvalues = new double[m_xdimension];
    m_xvalues.clear();
    for(int i = 0; i< m_xdimension; i++)
    {
        m_xvalues.push_back(X0 + i*delta_x); // m_xvalues[i] = X0 + i*delta_x;
    }
    //Prepare Y positions array

    //  vector <double> m_yvalues(m_ydimension);   // m_yvalues = new double[m_ydimension];
    m_yvalues.clear();
    for(int i = 0; i< m_ydimension; i++)
    {
        m_yvalues.push_back(Y0 + i*delta_y);// m_yvalues[i] = Y0 + i*delta_y;
    }
    //Prepare array for coefficients
    //Assume that points are equipartitioned
    m_delta_x = delta_x;
    m_delta_y = delta_y;
    //Allocate nonzero indices
    PrepareMatrixIndicesArray();
};


//Constructor
Bicubic_Interpolator::Bicubic_Interpolator(int NewXDimension, int NewYDimension, double delta_x, double delta_y, double X0, double Y0)
{
    Init(NewXDimension, NewYDimension, delta_x, delta_y, X0, Y0);
};

//Prepare nonzero matrix indices from the 16*16 BicubicMatrix
// The indices are calculated using i*100+j -this is useful because we can easily recognize i and j even if the matrix is only 16*16 and i*16+j would have worked but would be less easy to read!
void Bicubic_Interpolator::PrepareMatrixIndicesArray()
{
    m_nonzero_counter = 0;
    for(int i=0; i< 16; i++)
    {
        for(int j=0; j<16; j++)
        {
            if(Global::BicubicMatrix[i][j] != 0) m_nonzero_counter++;
        }
    }
    int index = 0;
    // m_nonzero_indices = new int[m_nonzero_counter];
    m_nonzero_indices.clear();
    for(int i=0; i< 16; i++)
    {
        for(int j=0; j<16; j++)
        {
            if(Global::BicubicMatrix[i][j] != 0)
            {
                m_nonzero_indices.push_back(i*100+j); // m_nonzero_indices[index] = i*100+j;
                index++;
            }
        }
    }
}

//Transform nonzero indices in a duplet i,j
void Bicubic_Interpolator::GetMatrixIndices(int index, int& i, int& j) const
{
    i = index / 100;
    j = index - 100*i;
};

//Get index of the point (just before) of the grid close to the X,Y point
// return the index near x,y, or (0,0) if one is out of the grid range
void Bicubic_Interpolator::GetIndex(const double& X, const double& Y, int &m_xindex, int &m_yindex) const
{
    //cout<<X<<" "<<Y<<" "<<m_xvalues[0]<<" "<<m_yvalues[0]<<" "<<m_delta_x;
    m_xindex = static_cast<int>(floor((X-m_xvalues[0])/m_delta_x));
    m_yindex = static_cast<int>(floor((Y-m_yvalues[0])/m_delta_y));
    EvaluateIndices(m_xindex, m_yindex);

};


//Show indices of the point of the grid close to the X,Y point
void Bicubic_Interpolator::ShowIndex(const double& X, const double& Y) const
{
    cout<<"Xindex = "<<static_cast<int>(floor((X-m_xvalues[0])/m_delta_x))<<endl;
    cout<<"Yindex = "<<static_cast<int>(floor((Y-m_yvalues[0])/m_delta_y))<<endl;
};

//Return false if  indices OUT OF RANGE and thus put BOTH of them at zero

bool Bicubic_Interpolator::EvaluateIndices(int &m_xindex, int &m_yindex) const
{
    if(m_xindex >= m_xdimension || m_xindex < 0 || m_yindex >= m_ydimension || m_yindex < 0)
    {
        m_xindex = 0 ;
        m_yindex = 0 ;
        return false;
    };
    return true;
};

//Transform matrix indices into vector indices
int Bicubic_Interpolator::ij_to_index(int i, int j) const
{
    return(i+4*j);
}

//Get interpolator name
const char* Bicubic_Interpolator::GetInterpolationType()
{
    return "Bicubic interpolator";
}


/* PARTIAL DISCRETE DERIVATIVES using f'(x) = (f(x+eps)-f(x-eps) )/2 eps between 2 grid points */
// cf http://en.wikipedia.org/wiki/Bicubic_interpolation
// because the derivation are not known we calculate them
double Bicubic_Interpolator::PartialXDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int j) const
{
    if(IndexX == 0) return 0;
    if(IndexX >= m_xdimension - 1) return 0;
    double DeltaF = matrix[j][(IndexX + 1)*m_ydimension + IndexY] - matrix[j][(IndexX - 1)*m_ydimension + IndexY];
    return (DeltaF / (2*m_delta_x));
};

double Bicubic_Interpolator::PartialYDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int j) const
{
    if(IndexY == 0) return 0;
    if(IndexY >= m_ydimension - 1) return 0;
    double DeltaF = matrix[j][IndexX *m_ydimension + IndexY+1] - matrix[j][IndexX * m_ydimension + IndexY-1];
    return (DeltaF / (2*m_delta_y));
};


double Bicubic_Interpolator::PartialXYDerivative(const vector<vector<double> > &matrix, int IndexX, int IndexY, int j) const
{
    if(IndexY == 0 || IndexX == 0) return 0;
    if(IndexY >= m_ydimension - 1 || IndexX >= m_xdimension - 1) return 0;
    double DeltaF = matrix[j][(IndexX + 1)*m_ydimension + IndexY + 1] + matrix[j][(IndexX - 1)*m_ydimension + IndexY - 1] -  matrix[j][(IndexX + 1)*m_ydimension + IndexY - 1] - matrix[j][(IndexX - 1)*m_ydimension + IndexY + 1];
    return (DeltaF / (4*m_delta_x*m_delta_y));
};



// Same but directly write from the matrix m[j][i]. Recall matrix[m_x][m_y] is transformed in m[j][m_x*m_ydimension + my]. So i = m_x*m_ydimension + my

double Bicubic_Interpolator::XDerivative(const vector<vector<double> > &matrix, int i, int j) const
{
    int IndexX = i/m_ydimension;
    if(IndexX == 0) return 0;
    if(IndexX == m_xdimension - 1) return 0;
    double DeltaF = matrix[j][i+m_ydimension] - matrix[j][i-m_ydimension];
    return (DeltaF / (2.*m_delta_x));
};

double Bicubic_Interpolator::YDerivative(const vector<vector<double> > &matrix, int i, int j) const
{
    int IndexY = i- m_ydimension*(i/m_ydimension);
    if(IndexY == 0) return 0;
    if(IndexY >= m_ydimension - 1) return 0;
    double DeltaF = matrix[j][i+1] - matrix[j][i-1];
    return (DeltaF / (2.*m_delta_y));
};


double Bicubic_Interpolator::XYDerivative(const vector<vector<double> > &matrix, int i, int j) const
{
    int IndexX = i/m_ydimension;
    int IndexY = i- m_ydimension*(i/m_ydimension);
    if(IndexY == 0 || IndexX == 0) return 0;
    if(IndexY >= m_ydimension - 1 || IndexX >= m_xdimension - 1) return 0;
    double DeltaF = matrix[j][i+m_ydimension  + 1] + matrix[j][i-m_ydimension - 1] -  matrix[j][i+m_ydimension - 1] - matrix[j][i-m_ydimension+ 1];
    return (DeltaF / (4*m_delta_x*m_delta_y));
};







/* PROTECTED MAIN FUNCTIONS */

/**  cf discussion in
LekienMarsdenMatrix.h
http://en.wikipedia.org/wiki/Bicubic_interpolation

Grouping the unknown parameters a_{ij} in a vector a, values of f, df/dx, df/dy, d2f/dxdy at 4 corner point in a vector F
the interpolation
f(x,y) = \sum_{i=0}^3 \sum_{j=0}^3 a_{ij} x^i y^j.
 can be reformulated into a matrix for the linear equation: BicubicMatrix a = F.

 Be CAREFUL The order of point is not the same than in Wikipedia
 in wikipedia it is (0,0) (1,0)   (0,1)   (1,1)
 Here it "seems"  to be (0,0) (1,0)   (1,1)   (0,1)
 because this comes from Numerical Recipies
cf http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c3-6.pdf
**/


void Bicubic_Interpolator::Prepare_Coefficients(const vector<vector<double> > &matrix, int m_xindex,  int m_yindex,     vector<double> &m_coefficients, int j, int type_field_read) const
{
// assuming that m_xindex, m_yindex have been already calculated
    double f[4];
    double dfdx[4], dfdy[4];
    double d2fdxdy[4];
    double x[16];
    int index_to_points_x[4] = {0, 1, 1, 0};
    int index_to_points_y[4] = {0, 0, 1, 1};
//std::cout<<
// filling arrays
    std::cout.setf(ios::fixed, ios::floatfield);

    for(int i = 0; i< 4; i++)
    {
        //cout<<m_xindex<<" "<<m_yindex<<" "<<m_zindex<<endl;
        //cout<<"derivatives "<<i<<" ok"<<endl;
        // matrix[m_x][m_y] is transformed in matrix[j][m_x*m_ydimension + my].

        f[i] = matrix[j][(m_xindex+index_to_points_x[i])*m_ydimension + m_yindex+index_to_points_y[i]];
        if (type_field_read == 2) // We have to calculate the derivatives
        {
            dfdx[i] = PartialXDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i],j);
            dfdy[i] = PartialYDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i],j);
            d2fdxdy[i] = PartialXYDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i],j);
        }
        if (type_field_read == 3) // The derivatives are part of the matrix: r	z	Br	Bz	dBr/dr	dBz/dr	dBr/dz	dBz/dz	d^2Br/drdz	d^2Bz/drdz
        {
            dfdx[i] = matrix[j+2][(m_xindex+index_to_points_x[i])*m_ydimension + m_yindex+index_to_points_y[i]];
            dfdy[i] = matrix[j+4][(m_xindex+index_to_points_x[i])*m_ydimension + m_yindex+index_to_points_y[i]];
            d2fdxdy[i] = matrix[j+6][(m_xindex+index_to_points_x[i])*m_ydimension + m_yindex+index_to_points_y[i]];
        }
        // Fill X as well
        x[0+i]=f[i];
        x[4+i]=dfdx[i]*m_delta_x;
        x[8+i]=dfdy[i]*m_delta_y;
        x[12+i]=d2fdxdy[i]*m_delta_x*m_delta_y;
    }
    std::cout.unsetf(ios::floatfield);
// calculate coefficients by means of 2D (similar to Lekien-Marsden) matrix

    int row = 0, column = 0;
//    for(int i = 0; i < 16; i++)
//    {
//        m_coefficients[i]= 0.0;
//    }
    for(int index = 0; index < m_nonzero_counter; index++)
    {
        GetMatrixIndices(m_nonzero_indices[index], row, column);
        m_coefficients[row] += Global::BicubicMatrix[row][column]*x[column];
    }
};




/* PUBLIC GET VALUE */
// Gives  f(x,y), f_x(x,y), f_y(x,y)
bool Bicubic_Interpolator::GetValue(const vector<vector<double> > &matrix, const double &x, const double &y,double &ans, double &dervx, double &dervy, int column_number, int type_field_read) const
{
    ans = 0.; // answer ! that is f(x,y)
    dervx = 0.; // derivation along x that is  f_x(x,y)
    dervy = 0.;  // derivation along y that is  f_y(x,y)
    // calculate LL indices basing on the point
    int m_xindex, m_yindex;
    GetIndex(x,y,m_xindex, m_yindex);
    if(m_xindex >= m_xdimension || m_xindex < 0 || m_yindex >= m_ydimension || m_yindex < 0)
    {
        cout<<endl<<"WARNING"<<endl<<"Exception::Bounds broken: "<<x<<" "<<y<<endl;
        return false;
    }
    // prepare coefficients for the selected matrix
    vector<double> m_coefficients(16,0.);
    Prepare_Coefficients(matrix,m_xindex, m_yindex,m_coefficients,column_number,type_field_read);
    double xd = (x - m_xvalues[m_xindex])/m_delta_x;
    double yd = (y - m_yvalues[m_yindex])/m_delta_y;

    // obtain values
    // COMPACT PROCEDURE: 1 cycle only (saves 12 operation per interpolation) - see commented lines for the explicit version
    for (int i=3; i>=0; i--)
    {
        ans = xd*(ans)+((m_coefficients[ij_to_index(3,i)]*yd + m_coefficients[ij_to_index(2,i)])*yd + m_coefficients[ij_to_index(1,i)])*yd + m_coefficients[ij_to_index(0,i)];
        dervy = xd*(dervy)+(3.0*m_coefficients[ij_to_index(3,i)]*yd + 2.0*m_coefficients[ij_to_index(2,i)])*yd + m_coefficients[ij_to_index(1,i)];
        dervx = yd*(dervx)+(3.0*m_coefficients[ij_to_index(i,3)]*xd + 2.0*m_coefficients[ij_to_index(i,2)])*xd + m_coefficients[ij_to_index(i, 1)];
    }
    dervx /= m_delta_x;
    dervy /= m_delta_y;

    /*
    //EXPLICIT VERSION: clearer, helds the same result but proves slower
    double arrayX[5] = {0, 1, xd, xd*xd, xd*xd*xd};
    double arrayY[5] = {0, 1, yd, yd*yd, yd*yd*yd};
    for(int j = 0;j < 4;j++) {
    		for(int i = 0;i < 4;i++) {
    		double coeff = m_coefficients[ij_to_index(j,i)];
    		//cout<<coeff<<endl;
    		ans += coeff*arrayX[i+1]*arrayY[j+1];
    		dervx += (i)*coeff*arrayX[i]*arrayY[j+1]/m_delta_x;
    		dervy += (j)*coeff*arrayX[i+1]*arrayY[j]/m_delta_y;
    		}
    }
    */
    return true;
};


// same but faster if only f(x,y) is needed
double Bicubic_Interpolator::GetValue(const vector<vector<double> > &matrix, const double &x, const double &y, int column_number, int type_field_read) const
{
    double ans = 0; // answer ! that is f(x,y)
    // calculate LL indices basing on the point
    int m_xindex, m_yindex;
    GetIndex(x,y,m_xindex, m_yindex);
    // prepare coefficients for the selected matrix
    vector<double> m_coefficients(16,0.);
    Prepare_Coefficients(matrix,m_xindex, m_yindex,m_coefficients,column_number,type_field_read);
    double xd = (x - m_xvalues[m_xindex])/m_delta_x;
    double yd = (y - m_yvalues[m_yindex])/m_delta_y;

    // obtain values
    // COMPACT PROCEDURE: 1 cycle only (saves 12 operation per interpolation) - see commented lines for the explicit version
    for (int i=3; i>=0; i--)
    {
        ans = xd*(ans)+((m_coefficients[ij_to_index(3,i)]*yd + m_coefficients[ij_to_index(2,i)])*yd + m_coefficients[ij_to_index(1,i)])*yd + m_coefficients[ij_to_index(0,i)];
    }
    return ans;
};


// same but faster if only f(x,y) is needed
// Version following http://en.wikipedia.org/wiki/Bicubic_interpolation
double Bicubic_Interpolator::GetValueInterpolated(const vector<vector<double> > &matrix, const double &x, const double &y, int column_number) const
{
    double ans = 0; // answer ! that is f(x,y)
    // calculate LL indices basing on the point
    int m_xindex, m_yindex;
    GetIndex(x,y,m_xindex, m_yindex);
    double xd = (x - m_xvalues[m_xindex])/m_delta_x; // TO be betwwen
    double yd = (y - m_yvalues[m_yindex])/m_delta_y;

    // obtain values
    // COMPACT PROCEDURE: 1 cycle only (saves 12 operation per interpolation) - see commented lines for the explicit version
    for (int i=3; i>=0; i--)
    {
//        ans = xd*(ans)+((m_coefficients[ij_to_index(3,i)]*yd + m_coefficients[ij_to_index(2,i)])*yd + m_coefficients[ij_to_index(1,i)])*yd + m_coefficients[ij_to_index(0,i)];
    }
    return ans;
    cerr << " the Bicubic_Interpolator::GetValueInterpolated function is not working (to be done) " << endl;
};


