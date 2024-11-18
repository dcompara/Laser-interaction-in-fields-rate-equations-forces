#include "TricubicInterpolation_module.h"
#include "LekienMarsdenMatrix.h"

using namespace std;

//Costructor #0
LM_Interpolator::LM_Interpolator() {};


//Costructor #1
LM_Interpolator::LM_Interpolator(int NewXDimension, int NewYDimension, int NewZDimension, double* Xvalues, double* Yvalues, double* Zvalues)
{
    Init(NewXDimension, NewYDimension, NewZDimension, Xvalues, Yvalues, Zvalues);
};

//Costructor #2
LM_Interpolator::LM_Interpolator(int NewXDimension, int NewYDimension, int NewZDimension, double delta_x, double delta_y, double delta_z, double X0, double Y0, double Z0)
{
    Init(NewXDimension, NewYDimension, NewZDimension, delta_x, delta_y, delta_z, X0, Y0, Z0);
};


void LM_Interpolator::Init(int NewXDimension, int NewYDimension, int NewZDimension, double* Xvalues, double* Yvalues, double* Zvalues)
{
    m_xdimension = NewXDimension;
    m_ydimension = NewYDimension;
    m_zdimension = NewZDimension;
    // Prepare X positions array
    m_xvalues = new double[m_xdimension];
    for(int i = 0; i< m_xdimension; i++)
    {
        m_xvalues[i] = Xvalues[i];
    }
    // Prepare Y positions array
    m_yvalues = new double[m_ydimension];
    for(int i = 0; i< m_ydimension; i++)
    {
        m_yvalues[i] = Yvalues[i];
    }
    // Prepare Z positions array
    m_zvalues = new double[m_zdimension];
    for(int i = 0; i< m_zdimension; i++)
    {
        m_zvalues[i] = Zvalues[i];
    }
    // Prepare array for coefficients
    m_coefficients = new double[64];
    // Assume that points are equipartitioned
    m_delta_x = m_xvalues[1]-m_xvalues[0];
    m_delta_y = m_yvalues[1]-m_yvalues[0];
    m_delta_z = m_zvalues[1]-m_zvalues[0];
    //Allocate nonzero indices
    PrepareLMIndicesArray();
};

void LM_Interpolator::Init(int NewXDimension, int NewYDimension, int NewZDimension, double delta_x, double delta_y, double delta_z, double X0, double Y0, double Z0)
{
    m_xdimension = NewXDimension;
    m_ydimension = NewYDimension;
    m_zdimension = NewZDimension;
    //Prepare X positions array
    m_xvalues = new double[m_xdimension];
    for(int i = 0; i< m_xdimension; i++)
    {
        m_xvalues[i] = X0 + i*delta_x;
    }
    //Prepare Y positions array
    m_yvalues = new double[m_ydimension];
    for(int i = 0; i< m_ydimension; i++)
    {
        m_yvalues[i] = Y0 + i*delta_y;
    }
    //Prepare Z positions array
    m_zvalues = new double[m_zdimension];
    for(int i = 0; i< m_zdimension; i++)
    {
        m_zvalues[i] = Z0 + i*delta_z;
    }
    //Prepare array for coefficients
    m_coefficients = new double[64];
    //Assume that points are equipartitioned
    m_delta_x = delta_x;
    m_delta_y = delta_y;
    m_delta_z = delta_z;

    /*
    cout<<m_xdimension<<" "<<m_ydimension<<" "<<m_zdimension<<endl;
    cout<<m_delta_x<<" "<<m_delta_y<<" "<<m_delta_z<<endl;
    cout<<m_xvalues[0]<<" "<<m_yvalues[0]<<" "<<m_zvalues[0]<<endl;
    */
    //Allocate nonzero indices
    PrepareLMIndicesArray();
};


//Destructor
LM_Interpolator::~LM_Interpolator()
{
    delete [] m_xvalues;
    delete [] m_zvalues;
    delete [] m_yvalues;
    delete [] m_coefficients;
    delete [] m_nonzero_indices;
}

//Prepare LM matrix indices
void LM_Interpolator::PrepareLMIndicesArray()
{
    m_nonzero_counter = 0;
    for(int i=0; i< 64; i++)
    {
        for(int j=0; j<64; j++)
        {
            if(Global::LekienMarsdenMatrix[i][j] != 0) m_nonzero_counter++;
        }
    }
    int index = 0;
    m_nonzero_indices = new int[m_nonzero_counter];
    for(int i=0; i< 64; i++)
    {
        for(int j=0; j<64; j++)
        {
            if(Global::LekienMarsdenMatrix[i][j] != 0)
            {
                m_nonzero_indices[index] = i*100+j;
                index++;
            }
        }
    }
}

//Transform nonzero indices in a duplet i,j
void LM_Interpolator::GetLMIndices(int index, int& i, int& j)
{
    i = index / 100;
    j = index - 100*i;
};

//Get indices
void LM_Interpolator::GetIndex(const double& X, const double& Y, const double& Z)
{
    m_xindex = static_cast<int>(floor((X-m_xvalues[0])/m_delta_x));
    m_yindex = static_cast<int>(floor((Y-m_yvalues[0])/m_delta_y));
    m_zindex = static_cast<int>(floor((Z-m_zvalues[0])/m_delta_z));
    //cout<<m_xindex<<" "<<m_yindex<<" "<<m_zindex<<endl;
    //cout<<X-m_xvalues[0]<<" "<<m_delta_x<<" "<<(X-m_xvalues[0])/m_delta_x<<endl;
    //cout<<Y-m_yvalues[0]<<" "<<m_delta_y<<" "<<(Y-m_yvalues[0])/m_delta_y<<endl;
};

//Evaluate indices - return false if OUT OF RANGE
bool LM_Interpolator::EvaluateIndices(const double& X, const double& Y, const double& Z)
{
    GetIndex(X,Y,Z);
    if(m_xindex >= m_xdimension || m_xindex < 0 || m_yindex >= m_ydimension || m_yindex < 0 || m_zindex >= m_zdimension || m_zindex < 0)
    {
        //cout<<endl<<"WARNING"<<endl<<"Exception::Bounds broken: "<<X<<" "<<Y<<" "<<Z<<endl;
        return false;
    };
    return true;
};

//Transform matrix indices into vector indices
int LM_Interpolator::ijk_to_index(int i, int j, int k)
{
    return(i+4*j+16*k);
}

//Get interpolator name
const char* LM_Interpolator::GetInterpolationType()
{
    return "Lekien-Marsden tricubic interpolator";
}

/* PARTIAL DISCRETE DERIVATIVES */

double LM_Interpolator::PartialXDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ)
{
    if(IndexX == 0) return 0;
    double DeltaF = matrix[IndexX + 1][IndexY][IndexZ] - matrix[IndexX - 1][IndexY][IndexZ];
    return (DeltaF / (2*m_delta_x));
};

double LM_Interpolator::PartialYDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ)
{
    if(IndexY == 0) return 0;
    double DeltaF = matrix[IndexX][IndexY+1][IndexZ] - matrix[IndexX][IndexY-1][IndexZ];
    return (DeltaF / (2*m_delta_y));
};

double LM_Interpolator::PartialZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ)
{
    if(IndexZ == 0) return 0;
    double DeltaF = matrix[IndexX][IndexY][IndexZ+1] - matrix[IndexX][IndexY][IndexZ-1];
    return (DeltaF / (2*m_delta_z));
};

double LM_Interpolator::PartialXYDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ)
{
    if(IndexX == 0 || IndexY == 0) return 0;
    double DeltaF = matrix[IndexX + 1][IndexY + 1][IndexZ] + matrix[IndexX - 1][IndexY - 1][IndexZ] -  matrix[IndexX + 1][IndexY - 1][IndexZ] - matrix[IndexX - 1][IndexY + 1][IndexZ];
    return (DeltaF / (4*m_delta_x*m_delta_y));
};

double LM_Interpolator::PartialXZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ)
{
    if(IndexX == 0 || IndexZ == 0) return 0;
    double DeltaF = matrix[IndexX + 1][IndexY][IndexZ+1] + matrix[IndexX - 1][IndexY][IndexZ - 1] -  matrix[IndexX + 1][IndexY][IndexZ - 1] - matrix[IndexX - 1][IndexY][IndexZ + 1];
    return (DeltaF / (4*m_delta_x*m_delta_z));
};

double LM_Interpolator::PartialYZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ)
{
    if(IndexZ == 0 || IndexY == 0) return 0;
    double DeltaF = matrix[IndexX][IndexY + 1][IndexZ + 1] + matrix[IndexX][IndexY - 1][IndexZ - 1] -  matrix[IndexX][IndexY - 1][IndexZ + 1] - matrix[IndexX][IndexY + 1][IndexZ - 1];
    return (DeltaF / (4*m_delta_y*m_delta_z));
};

/*
double LM_Interpolator::PartialXYZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ){
	double DeltaF = matrix[IndexX + 1][IndexY + 1][IndexZ + 1] - matrix[IndexX - 1][IndexY - 1][IndexZ - 1];
	DeltaF += matrix[IndexX - 1][IndexY - 1][IndexZ + 1] + matrix[IndexX + 1][IndexY - 1][IndexZ - 1] + matrix[IndexX - 1][IndexY + 1][IndexZ - 1];
	DeltaF -= matrix[IndexX - 1][IndexY + 1][IndexZ + 1] + matrix[IndexX + 1][IndexY - 1][IndexZ + 1] + matrix[IndexX + 1][IndexY + 1][IndexZ - 1];
	return (DeltaF / (8*m_delta_x*m_delta_y*m_delta_z));
};
*/

double LM_Interpolator::PartialXYZDerivative(double ***matrix, int IndexX, int IndexY, int IndexZ)
{
    if(IndexX == 0 || IndexY == 0 || IndexZ == 0) return 0;
    double DeltaF = matrix[IndexX + 1][IndexY + 1][IndexZ + 1] - matrix[IndexX - 1][IndexY - 1][IndexZ - 1];
    DeltaF += matrix[IndexX][IndexY + 1][IndexZ - 1] + matrix[IndexX - 1][IndexY - 1][IndexZ + 1] + matrix[IndexX + 1][IndexY][IndexZ - 1];
    DeltaF -= matrix[IndexX][IndexY + 1][IndexZ + 1] + matrix[IndexX + 1][IndexY][IndexZ + 1] + matrix[IndexX + 1][IndexY + 1][IndexZ - 1];
    return (DeltaF / (8*m_delta_x*m_delta_y*m_delta_z));
};

/* PROTECTED MAIN FUNCTIONS */

void LM_Interpolator::Prepare_Coefficients(double*** matrix)
{
    //cout<<m_xindex<<" "<<m_yindex<<" "<<m_zindex<<endl;
    // assuming that m_xindex, m_yindex and m_zindex have been already calculated
    double f[8];
    double dfdx[8], dfdy[8], dfdz[8];
    double d2fdxdy[8], d2fdxdz[8], d2fdydz[8];
    double d3fdxdydz[8];
    double x[64];
    int index_to_points_x[8] = {0, 1, 0, 1, 0, 1, 0, 1};
    int index_to_points_y[8] = {0, 0, 1, 1, 0, 0, 1, 1};
    int index_to_points_z[8] = {0, 0, 0, 0, 1, 1, 1, 1};
    // filling arrays
    std::cout.setf(ios::fixed, ios::floatfield);

    //DEBUG
    //cout<<"Point: "<<m_xindex<<" "<<m_yindex<<" "<<m_zindex<<endl;

    for(int i = 0; i< 8; i++)
    {
        //cout<<m_xindex<<" "<<m_yindex<<" "<<m_zindex<<endl;
        //cout<<"derivatives "<<i<<" ok"<<endl;
        f[i] = matrix[m_xindex+index_to_points_x[i]][m_yindex+index_to_points_y[i]][m_zindex+index_to_points_z[i]];
        dfdx[i] = PartialXDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i], m_zindex+index_to_points_z[i]);
        dfdy[i] = PartialYDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i], m_zindex+index_to_points_z[i]);
        dfdz[i] = PartialZDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i], m_zindex+index_to_points_z[i]);
        d2fdxdy[i] = PartialXYDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i], m_zindex+index_to_points_z[i]);
        d2fdxdz[i] = PartialXZDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i], m_zindex+index_to_points_z[i]);
        d2fdydz[i] = PartialYZDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i], m_zindex+index_to_points_z[i]);
        d3fdxdydz[i] = PartialXYZDerivative(matrix, m_xindex+index_to_points_x[i], m_yindex+index_to_points_y[i], m_zindex+index_to_points_z[i]);

        //DEBUG
        /*
        cout<<"::"<<i+1<<"::"<<endl;
        cout<<f[i]<<" "<<matrix[m_xindex+index_to_points_x[i]-1][m_yindex+index_to_points_y[i]][m_zindex+index_to_points_z[i]]<<" "<<matrix[m_xindex+index_to_points_x[i]+1][m_yindex+index_to_points_y[i]][m_zindex+index_to_points_z[i]]<<endl;
        cout<<dfdx[i]<<" "<<dfdy[i]<<" "<<dfdz[i]<<endl;
        cout<<d2fdxdy[i]<<" "<<d2fdxdz[i]<<" "<<d2fdydz[i]<<endl;
        cout<<d3fdxdydz[i]<<endl<<endl;;
        */

        // Fill X as well
        x[0+i]=f[i];
        x[8+i]=dfdx[i]*m_delta_x;
        x[16+i]=dfdy[i]*m_delta_y;
        x[24+i]=dfdz[i]*m_delta_z;
        x[32+i]=d2fdxdy[i]*m_delta_x*m_delta_y;
        x[40+i]=d2fdxdz[i]*m_delta_x*m_delta_z;
        x[48+i]=d2fdydz[i]*m_delta_y*m_delta_z;
        x[56+i]=d3fdxdydz[i]*m_delta_x*m_delta_y*m_delta_z;
    }
    std::cout.unsetf(ios::floatfield);

    // calculate coefficients by means of Lekien-Marsden matrix
    /*
    for (int i=0;i<64;i++) {
    	m_coefficients[i]= 0.0;
    	for (int j=0;j<64;j++) {
      		m_coefficients[i] += Global::LekienMarsdenMatrix[i][j]*x[j];
    	}
    }
    */
    int row = 0, column = 0;
    for(int i = 0; i < 64; i++)
    {
        m_coefficients[i]= 0.0;
    }
    for(int index = 0; index < m_nonzero_counter; index++)
    {
        GetLMIndices(m_nonzero_indices[index], row, column);
        m_coefficients[row] += Global::LekienMarsdenMatrix[row][column]*x[column];
    }
};

/* PUBLIC GET VALUE */

bool LM_Interpolator::GetValue(double*** matrix, const double &x, const double &y, const double &z, double &ans, double &dervx, double &dervy, double &dervz)
{
    ans = 0;
    dervx = 0;
    dervy = 0;
    dervz = 0;
    //cout<<m_xindex<<" "<<m_yindex<<" "<<m_zindex<<endl;
    // calculate lower limit indices according to the point
    GetIndex(x,y,z);
    //look for exceptions
    if(m_xindex >= m_xdimension || m_xindex < 0 || m_yindex >= m_ydimension || m_yindex < 0 || m_zindex >= m_zdimension || m_zindex < 0)
    {
        cout<<endl<<"WARNING"<<endl<<"Exception::Bounds broken: "<<x<<" "<<y<<" "<<z<<endl;
        return false;
    }
    // prepare coefficients for the selected matrix
    Prepare_Coefficients(matrix);
    double xd = (x - m_xvalues[m_xindex])/m_delta_x;
    double yd = (y - m_yvalues[m_yindex])/m_delta_y;
    double zd = (z - m_zvalues[m_zindex])/m_delta_z;
    double arrayX[5] = {0, 1, xd, xd*xd, xd*xd*xd};
    double arrayY[5] = {0, 1, yd, yd*yd, yd*yd*yd};
    double arrayZ[5] = {0, 1, zd, zd*zd, zd*zd*zd};
    // obtain values
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                double coeff = m_coefficients[ijk_to_index(i,j,k)];
                ans += coeff*arrayX[i+1]*arrayY[j+1]*arrayZ[k+1];
                dervx += (i)*coeff*arrayX[i]*arrayY[j+1]*arrayZ[k+1]/m_delta_x;
                dervy += (j)*coeff*arrayX[i+1]*arrayY[j]*arrayZ[k+1]/m_delta_y;
                dervz += (k)*coeff*arrayX[i+1]*arrayY[j+1]*arrayZ[k]/m_delta_z;
            }
        }
    };
    return true;
    //cout<<ans<<" "<<dervz<<endl;
};


double LM_Interpolator::GetValue(double*** matrix, const double &x, const double &y, const double &z)
{
    double ans = 0;
    // calculate lower limit indices according to the point
    GetIndex(x,y,z);
    // prepare coefficients for the selected matrix
    Prepare_Coefficients(matrix);
    double xd = (x - m_xvalues[m_xindex])/m_delta_x;
    double yd = (y - m_yvalues[m_yindex])/m_delta_y;
    double zd = (z - m_zvalues[m_zindex])/m_delta_z;
    double arrayX[5] = {0, 1, xd, xd*xd, xd*xd*xd};
    double arrayY[5] = {0, 1, yd, yd*yd, yd*yd*yd};
    double arrayZ[5] = {0, 1, zd, zd*zd, zd*zd*zd};
    // obtain values
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            for(int k = 0; k < 4; k++)
            {
                double coeff = m_coefficients[ijk_to_index(i,j,k)];
                ans += coeff*arrayX[i+1]*arrayY[j+1]*arrayZ[k+1];
            }
        }
    };
    return ans;
};


