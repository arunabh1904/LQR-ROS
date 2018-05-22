#include "integrationmethod.h"
#include <iostream>
using std::cout;
using std::endl;


// The predefined methods
// IMPLICIT METHODS
const IntegrationMethod BackwordEuler(BACKWORDS_EULER);
const IntegrationMethod Lobbato_IIIA_4th( LOBBATO_IIIA_4TH );
const IntegrationMethod Lobbato_IIIA_2nd( LOBBATO_IIIA_2ND );
const IntegrationMethod Lobbato_IIIB_4th( LOBBATO_IIIB_4TH );
const IntegrationMethod Lobbato_IIIB_2nd( LOBBATO_IIIB_2ND );
const IntegrationMethod Lobbato_IIIC_4th( LOBBATO_IIIC_4TH );
const IntegrationMethod Lobbato_IIIC_2nd( LOBBATO_IIIC_2ND );
const IntegrationMethod Randau_IA_3rd( RANDAU_IA_3RD );
const IntegrationMethod Randau_IA_5th( RANDAU_IA_5TH );
const IntegrationMethod Randau_IIA_3rd( RANDAU_IIA_3RD );
const IntegrationMethod Randau_IIA_5th( RANDAU_IIA_5TH );
const IntegrationMethod Gauss_Legendre_2nd( GAUSS_LEGENDRE_2ND );
const IntegrationMethod Gauss_Legendre_4th( GAUSS_LEGENDRE_4TH );
const IntegrationMethod Gauss_Legendre_6th( GAUSS_LEGENDRE_6TH );
const IntegrationMethod Gauss_Legendre_8th( GAUSS_LEGENDRE_8TH );
const IntegrationMethod Gauss_Legendre_10th( GAUSS_LEGENDRE_10TH );

// EXPLICIT METHODS
const IntegrationMethod Bogacki_Shampine_3rd( BOGACKI_SHAMPINE_23 );
const IntegrationMethod Runge_Kutta_3rd( RUNGE_KUTTA_3RD );
const IntegrationMethod Runge_Kutta_4th( RUNGE_KUTTA_4TH );
const IntegrationMethod Dormand_Prince_45( DORMAND_PRINCE_45 ); /**< Not Recomended For Interpolation **/
const IntegrationMethod Fehlberg_45( FEHLBERG_45 );
const IntegrationMethod Cash_Karp_45( CASH_KARP_45 );           /**< Not Recomended For Interpolation **/
const IntegrationMethod Cooper_Verner_8th( COOPER_VERNER_8TH ); /**< Not Recomended For Interpolation **/
const IntegrationMethod Curtis_8th( CURTIS_8TH ); /**< Not Great For Interpolation **/
const IntegrationMethod Dormand_Prince_8th( DORMAN_PRINCE_8TH );
const IntegrationMethod Rk_10( RK_10 );/**< Not Recomended For Interpolation **/
const IntegrationMethod Rk_10_8_Curtis( RK_10_8_CURTIS );
//static const IntegrationMethod Rk_12_10(RK_12_10); /**< Not Recomended For Interpolation **/


// A CONSTAINT DRIVEN NODE NO INTEGRATION PERFORMED
const IntegrationMethod ConstraintBased( CONSTRAINT_BASED );


std::string IntegrationMethod::name() const
{
    /**< Returns the method name **/
    return _name;
}
unsigned int IntegrationMethod::numStages() const
{
    /**< Returns the number of internode slopes **/
    return _a.cols();
}

unsigned int IntegrationMethod::order() const
{
    /**< Returns the order of the method **/
    return _order;
}

double IntegrationMethod::a(unsigned int i, unsigned int j) const
{
    /**< Returnds the a(i,j) slope multipler index **/
    return _a(i,j);
}

double IntegrationMethod::b(unsigned int i) const
{
    /**< Returns the b(i) index **/
    return _b(i);
}

double IntegrationMethod::b_errEst(unsigned int i) const
{
    /**< Returns the b(i) index **/
    return _b2(i);
}

double IntegrationMethod::c(unsigned int i) const
{
    /**< Returns the c(i) index **/
    return _c(i);
}

unsigned int IntegrationMethod::c_duplication_index(unsigned int i) const
{
    return c_duplicate[i];
}

Eigen::MatrixXd IntegrationMethod::buildInterpolationMatrix(double h) const
{
    // number of slopes
    unsigned int N_k = _c.rows();

    // count number of inter-node slopes that need to be resolved into states too
    unsigned int extraRows = 0;
    for( unsigned int i=0; i < N_k; i++ )
        if( _c(i) != 0 && _c(i) != 1 )
            extraRows ++;


    unsigned int matSize = (N_k+2); // number of states (Xo,Xf,K0..KN)
    unsigned int numCoeff = order()+1;


    Eigen::MatrixXd interpMatrix(matSize+extraRows,numCoeff);
    interpMatrix.setZero(matSize+extraRows,numCoeff);


    // upper left corner is identity
    interpMatrix(0,0) = 1;

    // slopes = to polynomial slopes
    for( unsigned int rowNum =0; rowNum< N_k; rowNum++ )
    {
        interpMatrix(1+rowNum,1) = 1;
        for( unsigned int colNum =0; (2+colNum) < numCoeff; colNum++ )
        {
            interpMatrix(1+rowNum, 2+colNum) = (colNum+2)*pow(_c(rowNum)*h,colNum+1);
        }
    }

    // final value is at end
    for( unsigned int colNum = 0; colNum < numCoeff; colNum++ )
    {
        interpMatrix(matSize-1, colNum) = pow(h,colNum);
    }

    // add rows for non-endpoint c's
    extraRows = 0;
    for( unsigned int i=0; i < N_k; i++ )
        if( _c(i) !=0 && _c(i) != 1 )
        {
            for( unsigned int colNum = 0; colNum < numCoeff; colNum++ )
            {
                interpMatrix(matSize+extraRows, colNum) = pow(_c(i)*h,colNum);
            }
            extraRows++;
        }


    // build matrix that maps to slopes and states at non-endpoint c's
    Eigen::MatrixXd stateMatrix = Eigen::MatrixXd::Zero(matSize+extraRows,matSize);
    stateMatrix.block(0,0,matSize,matSize).setIdentity(); // Xo, K0..K1, Xf

    extraRows = 0;
    for( unsigned int cCount =0; cCount < N_k; cCount ++ )
    {
        if( _c(cCount) != 0 && _c(cCount)!= 1)
        {
            stateMatrix(extraRows+matSize,0) = 1;
            stateMatrix.block(extraRows+matSize, 1, 1, N_k ) = h*_a.block(cCount,0,1,N_k);
            extraRows ++;
        }
    }

    //std::cout << "-----------" << std::endl <<std::endl;
    //std::cout << interpMatrix << std::endl << std::endl;
    //std::cout << "-----------" << std::endl <<std::endl;
    //std::cout << stateMatrix << std::endl;
    //std::cout << (interpMatrix.transpose()*interpMatrix).determinant() << std::endl;
    //std::cout << "-----------" << std::endl <<std::endl;
    //std::cout << interpMatrix.colPivHouseholderQr().solve(stateMatrix) << std::endl;


    return interpMatrix.colPivHouseholderQr().solve(stateMatrix);
}

bool IntegrationMethod::isAdaptive() const
{
    return _adaptive_method;
}
bool IntegrationMethod::first_same_as_last() const
{
    return _first_same_as_last;
}


IntegrationMethod::integratorType IntegrationMethod::type() const
{
    return _type;
}

IntegrationMethod::IntegrationMethod( IntegrationMethodName method_type )
{
    switch(method_type)
    {
    case BACKWORDS_EULER:
    {
        _name = "Backwords Euler";

        _order = 1;
        _numStages = 1;
        _type = IMPLICIT;
        _first_same_as_last = true;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1.0;

        _b.setZero(_numStages,1);
        _b << 1.0;

        _c.setZero(_numStages,1);
        _c << 1.0;
        break;
    }
    case  LOBBATO_IIIA_4TH:
    {
        _name = "Lobbato IIIA 4th Order";

        _order = 4;
        _numStages = 3;
        _type = IMPLICIT;
        _first_same_as_last = true;
        _adaptive_method = false;

        _a.setZero(3,3);
        _a << 0.0,          0.0,		0.0,
                5.0/24.0,     1.0/3.0,   -1.0/24.0,
                1.0 / 6.0,    2.0 / 3.0,	1.0 / 6.0;

        _b.setZero(3,1);
        _b << 1.0 / 6.0,
                2.0 / 3.0,
                1.0 / 6.0;

        _c.setZero(3,1);
        _c << 0.0,
                1.0 / 2.0,
                1.0;
        break;
    }
    case LOBBATO_IIIA_2ND:
    {
        _name = "Lobbato IIIA 2nd Order";

        _order = 2;
        _numStages = 2;
        _type = IMPLICIT;
        _first_same_as_last = true;
        _adaptive_method = false;


        _a.setZero(_numStages,_numStages);
        _a << 0.0,          0.0,
                0.5,          0.5;


        _b.setZero(_numStages,1);
        _b << 1.0 / 2.0,
                1.0 / 2.0;

        _c.setZero(_numStages,1);
        _c << 0,1;
        break;
    }
    case LOBBATO_IIIB_2ND:
    {
        _name = "Lobbato IIIB 2nd Order";

        _order = 2;
        _numStages = 2;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 0.0,   0.0,
                0.5,   0.0;


        _b.setZero(_numStages,1);
        _b << 1.0 / 2.0,
                1.0 / 2.0;

        _c.setZero(_numStages,1);
        _c << 0,1;
        break;
    }
    case LOBBATO_IIIB_4TH:
    {
        _name = "Lobbato IIIB 4th Order";

        _order = 4;
        _numStages = 3;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1.0/6.0,  -1.0/6.0,   0.0,
                1.0/6.0,   1.0/3.0,   0.0,
                1.0/6.0,   5.0/6.0,   0.0;


        _b.setZero(_numStages,1);
        _b << 1.0 / 6.0,
                2.0 / 3.0,
                1.0 / 6.0;

        _c.setZero(_numStages,1);
        _c << 0.0,
                1.0/2.0,
                1.0;

        break;
    }
    case LOBBATO_IIIC_2ND:
    {
        _name = "Lobbato IIIC 2nd Order";

        _order = 2;
        _numStages = 2;
        _type = IMPLICIT;
        _first_same_as_last = true;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1.0/2.0,  -1.0/2.0,
                1.0/2.0,   1.0/2.0;


        _b.setZero(_numStages,1);
        _b << 1.0 / 2.0,
                1.0 / 2.0;

        _c.setZero(_numStages,1);
        _c << 0.0,
                1.0;

        break;
    }
    case LOBBATO_IIIC_4TH:
    {
        _name = "Lobbato IIIC 4th Order";

        _order = 4;
        _numStages = 3;
        _type = IMPLICIT;
        _first_same_as_last = true;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1.0/6.0,  -1.0/3.0,   1.0/6.0,
                1.0/6.0,   5.0/12.0, -1.0/12.0,
                1.0/6.0,   2.0/3.0,   1.0/6.0;


        _b.setZero(_numStages,1);
        _b << 1.0 / 6.0,
                2.0 / 3.0,
                1.0 / 6.0;

        _c.setZero(_numStages,1);
        _c << 0.0,
                1.0/2.0,
                1.0;

        break;
    }

    case RANDAU_IA_3RD:
    {
        _name = "Randau IA 3rd Order";

        _order = 3;
        _numStages = 2;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1.0/4.0,  -1.0/4.0,
                1.0/4.0,   5.0/12.0;

        _b.setZero(_numStages,1);
        _b << 1.0 / 4.0,
                3.0 / 4.0;

        _c.setZero(_numStages,1);
        _c << 0.0,
                2.0/3.0;

        break;
    }
    case RANDAU_IA_5TH:
    {
        _name = "Randau IA 5th Order";

        _order = 5;
        _numStages = 3;
        _type = IMPLICIT;
        _first_same_as_last = true;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 11.0/45.0 - 7.0*sqrt(6.0)/360.0,      37.0/225.0 - 169.0*sqrt(6.0)/1800.0,    -2.0/225.0 + sqrt(6)/75.0,
                37.0/225.0 + 169.0*sqrt(6.0)/1800.0,  11.0/45.0 + 7.0*sqrt(6.0)/360.0,        -2.0/225.0 - sqrt(6)/75.0,
                4.0/9.0 - sqrt(6)/36.0,               4.0/9.0 + sqrt(6)/36.0,                 1.0/9.0;


        _b.setZero(_numStages,1);
        _b << 4.0 / 9.0 - sqrt(6)/36.0,
                4.0 / 9.0 + sqrt(6)/36.0,
                1.0 / 9.0;

        _c.setZero(_numStages,1);
        _c << 2.0/5.0 - sqrt(6)/10.0,
                2.0/5.0 + sqrt(6)/10.0,
                1.0;

        break;
    }

    case RANDAU_IIA_3RD:
    {
        _name = "Randau IIA 3rd Order";

        _order = 3;
        _numStages = 2;
        _type = IMPLICIT;
        _first_same_as_last = true;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 5.0/12.0, -1.0/12.0,
                3.0/4.0,   1.0/4.0;


        _b.setZero(_numStages,1);
        _b << 3.0 / 4.0,
                1.0 / 4.0;

        _c.setZero(_numStages,1);
        _c << 1.0/3.0,
                1.0;

        break;
    }
    case RANDAU_IIA_5TH:
    {
        _name = "Randau IIA 5th Order";

        _order = 5;
        _numStages = 3;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1.0/9.0,  (-1.0 - sqrt(6))/18.0,                (-1.0 + sqrt(6))/18.0,
                1.0/9.0,  (11.0/45.0 + (7.0*sqrt(6.0))/360.0),  (11.0/45.0 - (43.0*sqrt(6.0))/360.0),
                1.0/9.0,  (11.0/45.0 + (43.0*sqrt(6.0))/360.0), (11.0/45.0 - (7.0*sqrt(6.0))/360.0);


        _b.setZero(_numStages,1);
        _b << 1.0 / 9.0,
                4.0 / 9.0 + sqrt(6)/36.0,
                4.0 / 9.0 - sqrt(6)/36.0;

        _c.setZero(_numStages,1);
        _c << 0.0,
                3.0/5.0 - sqrt(6)/10.0,
                3.0/5.0 + sqrt(6)/10.0;

        break;
    }
    case GAUSS_LEGENDRE_2ND:
    {
        /*
         * Title: Implicit Runge-Kutta Processes
         * Author(s): J. C. Butcher
         * Source: Mathematics of Computation, Vol. 18, No. 85 (Jan., 1964), pp. 50-64
         * Stable URL: http://www.jstor.org/stable/2003405
         */

        _name = "Gauss Legendre 2nd Order";

        _order = 2;
        _numStages = 1;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1/2.0;


        _b.setZero(_numStages,1);
        _b << 1;

        _c.setZero(_numStages,1);
        _c << 1/2.0;

        break;
    }
    case GAUSS_LEGENDRE_4TH:
    {
        /*
         * Title: Implicit Runge-Kutta Processes
         * Author(s): J. C. Butcher
         * Source: Mathematics of Computation, Vol. 18, No. 85 (Jan., 1964), pp. 50-64
         * Stable URL: http://www.jstor.org/stable/2003405
         */

        _name = "Gauss Legendre 4th Order";

        _order = 6;
        _numStages = 2;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;

        _a.setZero(_numStages,_numStages);
        _a << 1.0/4.0,                 1.0/4.0 - sqrt(3.0)/6.0,
                1.0/4.0 + sqrt(3.0)/6.0, 1.0/4.0;


        _b.setZero(_numStages,1);
        _b << 1.0 / 2.0,
                1.0 / 2.0;

        _c.setZero(_numStages,1);
        _c << 1.0 / 2.0 - sqrt(3)/6,
                1.0 / 2.0 + sqrt(3)/6;

        break;
    }
    case GAUSS_LEGENDRE_6TH:
    {
        /*
         * Title: Implicit Runge-Kutta Processes
         * Author(s): J. C. Butcher
         * Source: Mathematics of Computation, Vol. 18, No. 85 (Jan., 1964), pp. 50-64
         * Stable URL: http://www.jstor.org/stable/2003405
         */

        _name = "Gauss Legendre 6th Order";

        _order = 6;
        _numStages = 3;
        _type = IMPLICIT;
        _adaptive_method = false;
        _first_same_as_last = false;


        _a.setZero(_numStages,_numStages);
        _a << 5.0/36.0, 2.0/9.0-sqrt(15.0)/15.0, 5.0/36.0-sqrt(15.0)/30.0,
                5.0/36.0+sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0-sqrt(15.0)/24.0,
                5.0/36.0+sqrt(15.0)/30.0, 2.0/9.0+sqrt(15.0)/15.0, 5.0/36.0;


        _b.setZero(_numStages,1);
        _b << 5.0 / 18.0,
                4.0 / 9.0,
                5.0 / 18.0;

        _c.setZero(_numStages,1);
        _c << 1.0 / 2.0 - sqrt(15)/10.0,
                1.0 / 2.0,
                1.0 / 2.0 + sqrt(15)/10.0;

        break;
    }
    case GAUSS_LEGENDRE_8TH:
    {
        /*
         * Title: Implicit Runge-Kutta Processes
         * Author(s): J. C. Butcher
         * Source: Mathematics of Computation, Vol. 18, No. 85 (Jan., 1964), pp. 50-64
         * Stable URL: http://www.jstor.org/stable/2003405
         */

        _name = "Gauss Legendre 8th Order";

        _order = 8;
        _numStages = 4;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;


        double w1,w1p,w2,w2p,w3,w3p,w4,w4p,w5,w5p;
        w1  = 1.0/8.0-sqrt(30.0)/144.0;
        w1p = 1.0/8.0+sqrt(30.0)/144.0;
        w2  = 1.0/2.0*sqrt((15.0+2*sqrt(30.0))/35.0);
        w2p = 1.0/2.0*sqrt((15.0-2*sqrt(30.0))/35.0);
        w3  = w2 *(1.0/6.0+sqrt(30.0)/24.0);
        w3p = w2p*(1.0/6.0-sqrt(30.0)/24.0);
        w4  = w2 *(1.0/21.0+5.0*sqrt(30.0)/168.0);
        w4p = w2p*(1.0/21.0-5.0*sqrt(30.0)/168.0);
        w5  = w2  - 2*w3;
        w5p = w2p - 2*w3p;

        _a.setZero(_numStages,_numStages);
        _a << w1,       w1p-w3+w4p,     w1p-w3-w4p,     w1 - w5,
                w1-w3p+w4,    w1p,        w1p-w5p,        w1-w3p-w4,
                w1+w3p+w4, w1p+w5p,       w1p,            w1+w3p-w4,
                w1+w5,    w1p+w3+w4p,     w1p+w3-w4p,     w1;

        _b.setZero(_numStages,1);
        _b << 2.0*w1,
                2.0*w1p,
                2.0*w1p,
                2.0*w1;

        _c.setZero(_numStages,1);
        _c <<  0.5 - w2,
                0.5 - w2p,
                0.5 + w2p,
                0.5 + w2;

        break;
    }
    case GAUSS_LEGENDRE_10TH:
    {
        /*
         * Title: Implicit Runge-Kutta Processes
         * Author(s): J. C. Butcher
         * Source: Mathematics of Computation, Vol. 18, No. 85 (Jan., 1964), pp. 50-64
         * Stable URL: http://www.jstor.org/stable/2003405
         */

        _name = "Gauss Legendre 10th Order";

        _order = 10;
        _numStages = 5;
        _type = IMPLICIT;
        _first_same_as_last = false;
        _adaptive_method = false;


        double w1,w1p,w2,w2p,w3,w3p,w4,w4p,w5,w5p, w6, w6p, w7 ,w7p;
        w1  = (322.0-13.0*std::sqrt(70.0))/3600.0;
        w1p = (322.0+13.0*std::sqrt(70.0))/3600.0;
        w2  = 1.0/2.0*std::sqrt((35.0+2.0*sqrt(70.0))/63.0);
        w2p = 1.0/2.0*std::sqrt((35.0-2.0*sqrt(70.0))/63.0);
        w3  = w2 *((452.0+59.0*std::sqrt(70.0))/3240.0);
        w3p = w2p*((452.0-59.0*std::sqrt(70.0))/3240.0);
        w4  = w2 *((64.0+11.0*std::sqrt(70.0))/1080.0);
        w4p = w2p*((64.0-11.0*std::sqrt(70.0))/1080.0);
        w5  = 8.0*w2 *((23.0-std::sqrt(70.0))/405.0);
        w5p = 8.0*w2p*((23.0+std::sqrt(70.0))/405.0);
        w6  = w2 - 2.0*w3 -w5;
        w6p = w2p- 2.0*w3p-w5p;
        w7  = w2 *((308.0-23.0*std::sqrt(70.0))/960.0);
        w7p  = w2p*((308.0+23.0*std::sqrt(70.0))/960.0);

        _a.setZero(_numStages,_numStages);
        _a << w1,        w1p-w3+w4p,    32.0/225.0 - w5,  w1p-w3-w4p,     w1 - w6,
              w1-w3p+w4, w1p,           32.0/225.0 - w5p, w1p-w6p,        w1-w3p-w4,
              w1 + w7,   w1p+w7p,       32.0/225.0,       w1p-w7p,        w1-w7,
              w1+w3p+w4, w1p+w6p,       32.0/225.0 + w5p, w1p,            w1+w3p-w4,
              w1+w6,     w1p+w3+w4p,    32.0/225.0 + w5,  w1p+w3-w4p,      w1;

        _b.setZero(_numStages,1);
        _b << 2.0*w1,
                2.0*w1p,
                64.0/225.0,
                2.0*w1p,
                2.0*w1;

        _c.setZero(_numStages,1);
        _c <<  0.5 - w2,
                0.5 - w2p,
                0.5,
                0.5 + w2p,
                0.5 + w2;

        break;
    }
    case DORMAND_PRINCE_45:
    {
        // Coefficients for Dormand–Prince 4th and 5th order integration
        _name = "Dormand Prince 45";
        _order = 5;
        _numStages = 7;
        _type = EXPLICIT;
        _adaptive_method = true;
        _first_same_as_last = true;


        _a.setZero(_numStages,_numStages);
        _a << 0.0,               0.0,              0.0,            0.0,            0.0,                     0.0,    0.0,
                1.0/5.0,         0.0,              0.0,            0.0,            0.0,                     0.0,    0.0,
                3.0/40.0,        9.0/40.0,         0.0,            0.0,            0.0,                     0.0,    0.0,
                44.0/45.0,      -56.0/15.0,        32.0/9.0,       0.0,            0.0,                     0.0,    0.0,
                19372.0/6561.0, -25360.0/2187.0,   64448.0/6561.0, -212.0/729.0,   0.0,                     0.0,    0.0,
                9017.0/3168.0,  -355.0/33.0,       46732.0/5247.0, 49.0/176.0,     -5103.0/18656.0,         0.0,    0.0,
                35.0/384.0,      0.0,              500.0/1113.0,   125.0/192.0,    -2187.0/6784.0,    11.0/84.0,    0.0;


        _b.setZero(_numStages,1);
        _b << 35.0/384.0,
               0.0,
               500.0/1113.0,
               125.0/192.0,
              -2187.0/6784.0,
               11.0/84.0,
               0.0;

        _b2.setZero(_numStages,1);
        _b2 << 5179.0/57600.0,
               0.0,
               7571.0/16695.0,
               393.0/640.0,
              -92097.0/339200.0,
               187.0/2100.0,
               1.0/40.0;

        _c.setZero(_numStages,1);
        _c << 0.0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0;


        break;
    }
    case BOGACKI_SHAMPINE_23:
    {
        _name = "Bogacki Shampine";
        _order = 3;
        _numStages = 4;
        _type = EXPLICIT;


        _a.setZero(_numStages,_numStages);
        _a << 0.0,        0.0,     0.0,     0.0,
                1.0/2.0,  0.0,     0.0,     0.0,
                0.0,      3.0/4.0, 0.0,     0.0,
                2.0/9.0,  1.0/3.0, 4.0/9.0, 0.0;

        _b.setZero(_numStages,1);
        _b << 2.0/9.0,   1.0/3.0,    4.0/9.0, 0.0;

        _b2.setZero(_numStages,1);
        _b2 << 7.0/24.0, 	1.0/4.0, 	1.0/3.0, 1.0/8.0;

        _c.setZero(_numStages,1);
        _c << 0.0, 1.0/2.0, 3.0/4.0, 1.0;

        _first_same_as_last = true;
        _adaptive_method = true;
        break;
    }
    case FEHLBERG_45:
    {
        _name = "Fehlberg 45";
        _order = 5;
        _numStages = 6;
        _type = EXPLICIT;


        _a.setZero(_numStages,_numStages);
        _a << 0.0,          0.0,            0.0,            0.0,          0.0,       0.0,
                1.0/4.0,      0.0,            0.0,            0.0,          0.0,       0.0,
                3.0/32.0,     9.0/32.0,       0.0,            0.0,            0.0,     0.0,
                1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0,            0.0,     0.0,
                439.0/216.0,   -8.0,         3680.0/513.0,  -845.0/4140.0,    0.0,     0.0,
                -8.0/27.0,    2.0,          -3544.0/2565.0, 1859.0/4140.0, -11.0/40.0, 0.0;

        _b.setZero(_numStages,1);
        _b  << 16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0;
        _b2.setZero(_numStages,1);
        _b2 << 25.0/216.0, 0.0, 1408.0/2565.0,  2197.0/4140.0,   -1.0/5.0,  0.0;

        _c.setZero(_numStages,1);
        _c << 0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0;


        _first_same_as_last = false;
        _adaptive_method = true;
        break;
    }
    case CASH_KARP_45:
    {
        _name = "Cash Karp 45";
        _order = 5;
        _numStages = 6;
        _type = EXPLICIT;


        _a.setZero(_numStages,_numStages);
        _a << 0.0,            0.0,          0.0,           0.0,                 0.0,        0.0,
                1.0/5.0,        0.0,          0.0,           0.0,                 0.0,        0.0,
                3.0/40.0,       9.0/40.0,     0.0,           0.0,                 0.0,        0.0,
                3.0/10.0,      -9.0/10.0,     6.0/5.0,       0.0,                 0.0,        0.0,
                -11.0/54.0,      5.0/2.0,     -70.0/27.0,     35.0/27.0,           0.0,        0.0,
                1631.0/55296.0,  175.0/512.0,   575.0/13824.0, 44275.0/110592.0, 253.0/4096.0, 0.0;

        _b.setZero(_numStages,1);
        _b << 37.0/378.0,    0.0, 250.0/621.0,    125.0/594.0,    0.0,          512.0/1771.0;

        _b2.setZero(_numStages,1);
        _b2 << 2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0,  1.0/4.0;


        _c.setZero(_numStages,1);
        _c<< 0.0, 1.0/5.0, 3.0/10.0, 3.0/5.0, 1.0, 7.0/8.0;

        _first_same_as_last = false;
        _adaptive_method = true;
        break;
    }
    case RUNGE_KUTTA_4TH:
    {
        _name = "Runge Kutta 4th Order";
        _order = 4;
        _numStages = 4;
        _type = EXPLICIT;


        _a.setZero(_numStages,_numStages);
        _a << 0.0, 0.0, 0.0, 0.0,
                0.5, 0.0, 0.0, 0.0,
                0.0, 0.5, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0;

        _b.setZero(_numStages,1);
        _b << 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0;

        _c.setZero(_numStages,1);
        _c << 0.0, 0.5, 0.5, 1.0;

        _first_same_as_last = false;
        _adaptive_method = false;
        break;
    }
    case RUNGE_KUTTA_3RD:
    {
        _name = "Runge Kutta 3rd Order";
        _order = 3;
        _numStages = 3;
        _type = EXPLICIT;


        _a.setZero(_numStages,_numStages);
        _a << 0.0, 0.0, 0.0,
                0.5, 0.0, 0.0,
                -1.0, 2.0, 0.0;

        _b.setZero(_numStages,1);
        _b << 1.0/6.0, 2.0/3.0, 1.0/6.0;

        _c.setZero(_numStages,1);
        _c << 0.0, 0.5, 1.0;

        _first_same_as_last = false;
        _adaptive_method = false;
        break;
    }
    case COOPER_VERNER_8TH:
    {
        _name = "Cooper Verner 8th Order";
        _order = 8;
        _numStages = 11;
        _type = EXPLICIT;


        _a.setZero(_numStages,_numStages);
        _a <<    0.0,                      0.0,                        0.0,                            0.0,                            0.0,                                0.0,                              0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                0.5,                      0.0,                        0.0,                            0.0,                            0.0,                                0.0,                              0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                0.25,                     0.25,                       0.0,                            0.0,                            0.0,                                0.0,                              0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                1.0/7.0,                  (-7.0+3.0*sqrt(21.0))/98.0, (21.0-5.0*sqrt(21.0))/49.0,     0.0,                            0.0,                                0.0,                              0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                (11.0-sqrt(21.0))/84.0,   0.0,                        (18.0-4.0*sqrt(21.0))/63.0,     (21.0+sqrt(21.0))/252.0,        0.0,                                0.0,                              0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                (5.0-sqrt(21.0))/48.0,    0.0,                        (9.0-sqrt(21.0))/36.0,          (-231.0-14.0*sqrt(21.0))/360.0, (63.0+7.0*sqrt(21.0))/80.0,         0.0,                              0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                (10.0+sqrt(21.0))/42.0,   0.0,                        (-432.0-92.0*sqrt(21.0))/315.0, (633.0+145.0*sqrt(21.0))/90.0,  (-504.0-115.0*sqrt(21.0))/70.0,     (63.0+13.0*sqrt(21.0))/35.0,      0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                (1.0/14.0),               0.0,                        0.0,                            0.0,                            (14.0+3.0*sqrt(21.0))/126.0,        (13.0+3.0*sqrt(21.0))/63.0,       1.0/9.0,                          0.0,                            0.0,                            0.0,                        0.0,
                (1.0/32.0),               0.0,                        0.0,                            0.0,                            (91.0+21.0*sqrt(21.0))/576.0,       (11.0/72.0),                      (-385.0+75.0*sqrt(21.0))/1152.0,  (63.0-13.0*sqrt(21.0))/128.0,   0.0,                            0.0,                        0.0,
                (1.0/14.0),               0.0,                        0.0,                            0.0,                            1.0/9.0,                            (-733.0+147.0*sqrt(21.0))/2205.0, (515.0-111.0*sqrt(21.0))/504.0,   (-51.0+11.0*sqrt(21.0))/56.0,   (132.0-28.0*sqrt(21.0))/245.0,  0.0,                        0.0,
                0.0,                      0.0,                        0.0,                            0.0,                            (-42.0-7.0*sqrt(21.0))/18.0,        (-18.0-28.0*sqrt(21.0))/45.0,     (-273.0+53.0*sqrt(21.0))/72.0,    (301.0-53.0*sqrt(21.0))/72.0,   (28.0+28.0*sqrt(21.0))/45.0,    (49.0+7.0*sqrt(21.0))/18.0, 0.0;

        _b.setZero(_numStages,1);
        _b << 1.0/20.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                49.0/180.0,
                (16.0/45.0),
                49.0/180.0,
                1.0/20.0;

        _c.setZero(_numStages,1);
        _c << 0,
                0.5,
                0.5,
                (7.0-sqrt(21.0))/14.0,
                (7.0-sqrt(21.0))/14.0,
                0.5,
                (7.0+sqrt(21.0))/14.0,
                (7.0+sqrt(21.0))/14.0,
                0.5,
                (7.0-sqrt(21.0))/14.0,
                1.0;

        _first_same_as_last = false;
        _adaptive_method = false;
        break;
    }
    case CURTIS_8TH:
    {
        /*
         * Curtis, A.R. "An Eigth order Runge-Kutta Process with Elven Function Evaluations per Step", Numer. Math. 16, 268-277, 1970
         *
         */

        _name = "A. Curtis 8th Order";
        _order = 8;
        _numStages = 11;
        _type = EXPLICIT;

        // the following are arbitrary but set to the values recomended in the paper
        double c2, b9, b10;
        c2 = 0.183850407856;
        b9 = 1.0/5.0;
        b10 = 13.0/180.0;

        _a.setZero(_numStages,_numStages);
        _a <<   0.0,                                0.0,                        0.0,                            0.0,                                0.0,                                    0.0,                                                         0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                c2,                                 0.0,                        0.0,                            0.0,                                0.0,                                    0.0,                                                         0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                0.183850407856-0.016900486234/c2,   0.016900486234/c2,          0.0,                            0.0,                                0.0,                                    0.0,                                                         0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                0.068943902946,                     0.0,                        0.206831708839,                 0.0,                                0.0,                                    0.0,                                                         0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                0.689439029462,                     0.0,                       -2.585396360481,                 2.585396360481,                     0.0,                                    0.0,                                                         0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                0.082732683535,                     0.0,                        0.0,                            0.413663417677,                     0.330930734142,                         0.0,                                                         0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                0.097115513262,                     0.0,                        0.0,                            0.097308686537,                    -0.044005801467,                         0.022254766315,                                              0.0,                              0.0,                            0.0,                            0.0,                        0.0,
                -0.062320146988,                     0.0,                        0.0,                           -0.259948900289,                     0.249582954483,                        -0.113452633375,                                              0.686138726168,                   0.0,                            0.0,                            0.0,                        0.0,
                0.082732683535 + 0.045090916605/b9, 0.0,                        0.0,                            0.413663417677 + 0.021208416940/b9, 0.330930734142 - 0.218567298182/b9,     0.083972393539/b9,                                          -0.125928146333/b9,                0.194223717431/b9,              0.0,                            0.0,                        0.0,
                0.097115513262 + 0.003060721281/b10, 0.0,                       0.0,                            0.097308686537 + 0.001439603760/b10,-0.044005801467 - 0.014836105169/b10,   0.022254766315 + 0.005699952702/b10 - 0.046098506765*b9/b10,-0.008547862550/b10,               0.013183689976/b10,             0.046098506765*b9/b10,          0.0,                        0.0,
                -0.499040783934,                      0.0,                       0.0,                           -1.386394134872,                      1.331109090575,                       -0.167728029907 - 2.531493157611*b9,                         2.314646450731 - 16.546536707080*b10, -0.592592592593,             2.531493157611*b9,              16.546536707080*b10,        0.0;

        _b.setZero(_numStages,1);
        _b << 1.0/20.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.272222222222 - b9,
                0.272222222222 - b10,
                0.355555555556,
                b9,
                b10,
                1.0/20.0;

        _c.setZero(_numStages,1);
        _c << 0.0,
                c2,
                0.183850407856,
                0.275775611785,
                0.689439029462,
                0.827326835354,
                0.172673164646,
                0.5,
                0.827326835354,
                0.172673164646,
                1.0;

        _first_same_as_last = false;
        _adaptive_method = false;
        break;
    }
    case DORMAN_PRINCE_8TH:
    {
        /*
             * P.J. Prince & J.R. Dorman (1981)
             * High order embedded Runge-Kutta formulae. J.Comp. Appl. Math., Vol. 7. p.67-75.
             *
            */

        _name = "Dorman Prince 8th Order";
        _order = 8;
        _numStages = 13;
        _type = EXPLICIT;

        _a.setZero(_numStages,_numStages);
        _a <<   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
                1.0/18.0,                 0,        0,          0,                              0,                          0,                          0, 0, 0, 0, 0, 0, 0,
                1.0/48.0,                 1.0/16.0, 0,          0,                              0,                          0,                          0, 0, 0, 0, 0, 0, 0,
                1.0/32.0,                 0,        3.0/32.0,   0,                              0,                          0,                          0, 0, 0, 0, 0, 0, 0,
                5.0/16.0,                 0,       -75.0/64.0,  75.0/64.0,                      0,                          0,                          0, 0, 0, 0, 0, 0, 0,
                3.0/80.0,                 0,        0,          3.0/16.0,                       3.0/20.0,                   0,                          0, 0, 0, 0, 0, 0, 0,
                29443841.0/614563906.0,   0,        0,          77736538.0/692538347.0,         -28693883.0/1125000000.0,  23124283.0/1800000000.0,     0, 0, 0, 0, 0, 0, 0,
                16016141.0/946692911.0,   0,        0,          61564180.0/158732637.0,         22789713.0/633445777.0,    545815736.0/2771057229.0,   -180193667.0/1043307555.0,   0, 0, 0, 0, 0, 0,
                39632708.0/573591083.0,   0,        0,          -433636366.0/683701615.0,       -421739975.0/2616292301.0, 100302831.0/723423059.0,     790204164.0/839813087.0,    800635310.0/3783071287.0, 0, 0, 0, 0, 0,
                246121993.0/1340847787.0, 0,        0,          -37695042795.0/15268766246.0,   -309121744.0/1061227803.0, -12992083.0/490766935.0,     6005943493.0/2108947869.0,  393006217.0/1396673457.0,   123872331.0/1001029789.0, 0, 0, 0, 0,
                -1028468189.0/846180014.0, 0,        0,          8478235783.0/508512852.0,       1311729495.0/1432422823.0, -10304129995.0/1701304382.0, -48777925059.0/3047939560.0, 15336726248.0/1032824649.0,-45442868181.0/3398467696.0, 3065993473.0/597172653.0, 0, 0, 0,
                185892177.0/718116043.0,  0,        0,          -3185094517.0/667107341.0,      -477755414.0/1098053517.0, -703635378.0/230739211.0,    5731566787.0/1027545527.0,  5232866602.0/850066563.0,   -4093664535.0/808688257.0,  3962137247.0/1805957418.0, 65686358.0/487910083.0, 0, 0.0,
                403863854.0/491063109.0,  0,        0,          -5068492393.0/434740067.0,      -411421997.0/543043805.0,   652783627.0/914296604.0,    11173962825.0/925320556.0, -13158990841.0/6184727034.0, 3936647629.0/1978049680.0, -160528059.0/685178525.0,   248638103.0/1413531060.0, 0, 0.0;


        _b.setZero(_numStages,1);
        _b << 13451932.0/455176623.0,
                0,
                0,
                0,
                0,
                -808719846.0/976000145.0,
                1757004468.0/5645159321.0,
                656045339.0/265891186.0,
                -3867574721.0/1518517206.0,
                465885868.0/322736535.0,
                53011238.0/667516719.0,
                2.0/45.0,
                0.0;


        _b2.setZero(_numStages,1);
        _b2 << 14005451.0/335480064.0,
                0,
                0,
                0,
                0,
                -59238493.0/1068277825.0,
                181606767.0/758867731.0,
                561292985.0/797845732.0,
                -1041891430.0/1371343529.0,
                760417239.0/1151165299.0,
                118820643.0/751138087.0,
                -528747749.0/2220607170.0,
                1.0/4.0;

        _c.setZero(_numStages,1);
        _c << 0.0,
                1.0/18.0,
                1.0/12.0,
                1.0/8.0,
                5.0/16.0,
                3.0/8.0,
                59.0/400.0,
                93.0/200.0,
                5490023248.0/9719169821.0,
                13.0/20.0,
                1201146811.0/1299019798.0,
                1.0,
                1.0;

        _first_same_as_last = false;
        _adaptive_method = true;
        break;
    }

    case RK_10:
    {

        _name = "10th Order";
        _order = 10;
        _numStages = 17;
        _type = EXPLICIT;

        _a.setZero(_numStages,_numStages);
        _a.block<1,1>(1,0) <<   0.100000000000000000000000000000000000000000000000000000000000;
        _a.block<1,2>(2,0) <<   -0.915176561375291440520015019275342154318951387664369720564660,
                1.45453440217827322805250021715664459117622483736537873607016;
        _a.block<1,3>(3,0) <<  0.202259190301118170324681949205488413821477543637878380814562,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.606777570903354510974045847616465241464432630913635142443687;
        _a.block<1,4>(4,0) <<  0.184024714708643575149100693471120664216774047979591417844635,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.197966831227192369068141770510388793370637287463360401555746,
                -0.0729547847313632629185146671595558023015011608914382961421311;
        _a.block<1,5>(5,0) <<  0.0879007340206681337319777094132125475918886824944548534041378,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.410459702520260645318174895920453426088035325902848695210406,
                0.482713753678866489204726942976896106809132737721421333413261;
        _a.block<1,6>(6,0) <<  0.0859700504902460302188480225945808401411132615636600222593880,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.330885963040722183948884057658753173648240154838402033448632,
                0.489662957309450192844507011135898201178015478433790097210790,
                -0.0731856375070850736789057580558988816340355615025188195854775;
        _a.block<1,7>(7,0) <<  0.120930449125333720660378854927668953958938996999703678812621,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.260124675758295622809007617838335174368108756484693361887839,
                0.0325402621549091330158899334391231259332716675992700000776101,
                -0.0595780211817361001560122202563305121444953672762930724538856;
        _a.block<1,8>(8,0) <<  0.110854379580391483508936171010218441909425780168656559807038,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.0605761488255005587620924953655516875526344415354339234619466,
                0.321763705601778390100898799049878904081404368603077129251110,
                0.510485725608063031577759012285123416744672137031752354067590;
        _a.block<1,9>(9,0) <<  0.112054414752879004829715002761802363003717611158172229329393,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.144942775902865915672349828340980777181668499748506838876185,
                -0.333269719096256706589705211415746871709467423992115497968724,
                0.499269229556880061353316843969978567860276816592673201240332,
                0.509504608929686104236098690045386253986643232352989602185060;
        _a.block<1,10>(10,0) <<  0.113976783964185986138004186736901163890724752541486831640341,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.0768813364203356938586214289120895270821349023390922987406384,
                0.239527360324390649107711455271882373019741311201004119339563,
                0.397774662368094639047830462488952104564716416343454639902613,
                0.0107558956873607455550609147441477450257136782823280838547024,
                -0.327769124164018874147061087350233395378262992392394071906457;
        _a.block<1,11>(11,0) <<  0.0798314528280196046351426864486400322758737630423413945356284,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.0520329686800603076514949887612959068721311443881683526937298,
                -0.0576954146168548881732784355283433509066159287152968723021864,
                0.194781915712104164976306262147382871156142921354409364738090,
                0.145384923188325069727524825977071194859203467568236523866582,
                -0.0782942710351670777553986729725692447252077047239160551335016,
                -0.114503299361098912184303164290554670970133218405658122674674;
        _a.block<1,12>(12,0) <<  0.985115610164857280120041500306517278413646677314195559520529,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.330885963040722183948884057658753173648240154838402033448632,
                0.489662957309450192844507011135898201178015478433790097210790,
                -1.37896486574843567582112720930751902353904327148559471526397,
                -0.861164195027635666673916999665534573351026060987427093314412,
                5.78428813637537220022999785486578436006872789689499172601856,
                3.28807761985103566890460615937314805477268252903342356581925,
                -2.38633905093136384013422325215527866148401465975954104585807,
                -3.25479342483643918654589367587788726747711504674780680269911,
                -2.16343541686422982353954211300054820889678036420109999154887;
        _a.block<1,13>(13,0) <<  0.895080295771632891049613132336585138148156279241561345991710,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.197966831227192369068141770510388793370637287463360401555746,
                -0.0729547847313632629185146671595558023015011608914382961421311,
                0.0000000000000000000000000000000000000000000000000000000000000,
                -0.851236239662007619739049371445966793289359722875702227166105,
                0.398320112318533301719718614174373643336480918103773904231856,
                3.63937263181035606029412920047090044132027387893977804176229,
                1.54822877039830322365301663075174564919981736348973496313065,
                -2.12221714704053716026062427460427261025318461146260124401561,
                -1.58350398545326172713384349625753212757269188934434237975291,
                -1.71561608285936264922031819751349098912615880827551992973034,
                -0.0244036405750127452135415444412216875465593598370910566069132;
        _a.block<1,14>(14,0) << -0.915176561375291440520015019275342154318951387664369720564660,
                1.45453440217827322805250021715664459117622483736537873607016,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.777333643644968233538931228575302137803351053629547286334469,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.0910895662155176069593203555807484200111889091770101799647985,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0910895662155176069593203555807484200111889091770101799647985,
                0.777333643644968233538931228575302137803351053629547286334469;
        _a.block<1,15>(15,0) << 0.100000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.157178665799771163367058998273128921867183754126709419409654,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.157178665799771163367058998273128921867183754126709419409654;
        _a.block<1,16>(16,0) << 0.181781300700095283888472062582262379650443831463199521664945,
                0.675000000000000000000000000000000000000000000000000000000000,
                0.342758159847189839942220553413850871742338734703958919937260,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.259111214548322744512977076191767379267783684543182428778156,
                -0.358278966717952089048961276721979397739750634673268802484271,
                -1.04594895940883306095050068756409905131588123172378489286080,
                0.930327845415626983292300564432428777137601651182965794680397,
                1.77950959431708102446142106794824453926275743243327790536000,
                0.100000000000000000000000000000000000000000000000000000000000,
                -0.282547569539044081612477785222287276408489375976211189952877,
                -0.159327350119972549169261984373485859278031542127551931461821,
                -0.145515894647001510860991961081084111308650130578626404945571,
                -0.259111214548322744512977076191767379267783684543182428778156,
                -0.342758159847189839942220553413850871742338734703958919937260,
                -0.675000000000000000000000000000000000000000000000000000000000;


        _b.setZero(_numStages,1);
        _b << 0.181781300700095283888472062582262379650443831463199521664945,
                0.675000000000000000000000000000000000000000000000000000000000,
                0.342758159847189839942220553413850871742338734703958919937260,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.259111214548322744512977076191767379267783684543182428778156,
                -0.358278966717952089048961276721979397739750634673268802484271,
                -1.04594895940883306095050068756409905131588123172378489286080,
                0.930327845415626983292300564432428777137601651182965794680397,
                1.77950959431708102446142106794824453926275743243327790536000,
                0.100000000000000000000000000000000000000000000000000000000000,
                -0.282547569539044081612477785222287276408489375976211189952877,
                -0.159327350119972549169261984373485859278031542127551931461821,
                -0.145515894647001510860991961081084111308650130578626404945571,
                -0.259111214548322744512977076191767379267783684543182428778156,
                -0.342758159847189839942220553413850871742338734703958919937260,
                -0.675000000000000000000000000000000000000000000000000000000000,
                0.0;

        _b2.setZero(_numStages,1);
        _b2 << 0.0333333333333333333333333333333333333333333333333333333333333,
                0.0250000000000000000000000000000000000000000000000000000000000,
                0.0333333333333333333333333333333333333333333333333333333333333,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0500000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0400000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.189237478148923490158306404106012326238162346948625830327194,
                0.277429188517743176508360262560654340428504319718040836339472,
                0.277429188517743176508360262560654340428504319718040836339472,
                0.189237478148923490158306404106012326238162346948625830327194,
                -0.0400000000000000000000000000000000000000000000000000000000000,
                -0.0500000000000000000000000000000000000000000000000000000000000,
                -0.0333333333333333333333333333333333333333333333333333333333333,
                -0.0250000000000000000000000000000000000000000000000000000000000,
                0.0333333333333333333333333333333333333333333333333333333333333;

        _c.setZero(_numStages,1);
        _c << 0.000000000000000000000000000000000000000000000000000000000000,
                0.100000000000000000000000000000000000000000000000000000000000,
                0.539357840802981787532485197881302436857273449701009015505500,
                0.809036761204472681298727796821953655285910174551513523258250,
                0.309036761204472681298727796821953655285910174551513523258250,
                0.981074190219795268254879548310562080489056746118724882027805,
                0.833333333333333333333333333333333333333333333333333333333333,
                0.354017365856802376329264185948796742115824053807373968324184,
                0.882527661964732346425501486979669075182867844268052119663791,
                0.642615758240322548157075497020439535959501736363212695909875,
                0.357384241759677451842924502979560464040498263636787304090125,
                0.117472338035267653574498513020330924817132155731947880336209,
                0.833333333333333333333333333333333333333333333333333333333333,
                0.309036761204472681298727796821953655285910174551513523258250,
                0.539357840802981787532485197881302436857273449701009015505500,
                0.100000000000000000000000000000000000000000000000000000000000,
                1.00000000000000000000000000000000000000000000000000000000000;

        _first_same_as_last = true;
        _adaptive_method = true;
        break;

        //The estimate of the local truncation error is  (1/360)  h ( f(t1,x1)-f(t15,x15) )
    }
    case RK_10_8_CURTIS:
    {
        // http://arxiv.org/abs/1305.6165
        // http://msp.org/camcos/2014/9-2/p01.xhtml
        // A. R. Curtis. High-order explicit Runge-Kutta formulae, their uses, and limitations. IMA
        // Journal of Applied Mathematics, 16(1):3552, 1975.
        //
        //  Coefficients from: https://github.com/ketch/high_order_RK_RR/tree/master/code

        _name = "Curtis 10 and 8th Order";
        _order = 10;
        _numStages = 21;
        _type = EXPLICIT;

        _a.setZero(_numStages,_numStages);
        _a.block<1,1>(1,0) <<   0.1452518960316150517617548528770033320314511251329947060838468741983976455607179673401;

        _a.block<1,2>(2,0) <<   .7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1,
                .7262594801580752588087742643850166601572556256649735304192343709919882278035898367003e-1;

        _a.block<1,3>(3,0) <<  .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252e-1,
                0.0,
                .1634083830355669332319742094866287485353825157746190443443277334731973512558077132576;

        _a.block<1,4>(4,0) <<  .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252,
                0.0,
                -2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720,
                2.042604787944586665399677618582859356692281447182738054304096668414966890697596415720;

        _a.block<1,5>(5,0) <<  .6536335321422677329278968379465149941415300630984761773773109338927894050232308530303e-1,
                0.0,
                0.0,
                0.3268167660711338664639484189732574970707650315492380886886554669463947025116154265151,
                0.2614534128569070931711587351786059976566120252393904709509243735571157620092923412121;

        _a.block<1,6>(6,0) <<  0.8233707757482716585173454344310125296066814318521742241762319051772963627695955263034e-1,
                0.0,
                0.0,
                .2119171963202803561687843468555305553175658807629274312902985594840086570224567152664,
                -.3997343508054218311577932550061320162379840049816347807630118786107674477850206579628e-1,
                .2037865317596006197606259822674324543477946306275935376058802473310199901934015124941e-1;

        _a.block<1,7>(7,0) <<  .8595305779007343831562027786771081909543936570474230618236291858949564375587825985001e-1,
                0.0,
                0.0,
                0.0,
                0.0,
                .2911769478058850960337179621761553399856026049598393013981874594942289837064329700000,
                .3964475145147024104912442607655312127779123206780991480607158892217361663405931088001;

        _a.block<1,8>(8,0) <<  .8612093485606967549983047372292119178898514571588438099912534616486575243508895957628e-1,
                0.0,
                0.0,
                0.0,
                0.0,
                .1397464826824442089036313891001189801074425314582326737716288563521183595455090268480,
                .3951098495815674599900526056001284215294125840404176924334653987770478924197803010468,
                -.4079412703708563576307759281612056453162454270752418047326990081493640904820003348350e-1;

        _a.block<1,9>(9,0) <<  .7233144422337948077616348229119326315582930871089020733092900891206129381937795204778e-1,
                0,
                0,
                0,
                0,
                .2200276284689998102140972735735070061373242800181187459951219347361114857342828430157,
                .8789533425436734013369780264792573637952226487753296416823846876217040795688489371334e-1,
                -.4445383996260350863990674880611108986832860648196030000580004690002268108984238641730e-1,
                -.2183282289488754689095532966861839909872150913926337371522805434288481649401165594213;

        _a.block<1,10>(10,0) <<  .8947100936731114228785441966773836169071038390882857211057269158522704971585365845223e-1,
                0,
                0,
                0,
                0,
                .3946008170285561860741397654755022300929434262701385530048127140223687993778661654316,
                .3443011367963333487713764986067104675654371857504670290688086760696354596195596354011,
                -.7946682664292661290694938113119430997053815140863772328764150866582492425892231395780e-1,
                -.3915218947895966123834967996391962853380545808840091268064277812752553499114569444180,
                0;

        _a.block<1,11>(11,0) <<  .3210006877963209212945282736072241886741425314298532400216927262619488479186214523312e-1,
                0,
                0,
                0,
                0,
                0,
                0,
                -.1846375997512050141835163881753227910996323204749769226655464078048769505209525299752e-3,
                .1560894025313219860759149162557283383430181475726228517203663063649626288079337909898,
                .1934496857654560252749984220385188727138526287670744309970093278715606577140084022992,
                .2611612387636636496908928477536452288263163392010050661129958478089356710938164130987;

        _a.block<1,12>(12,0) <<  .4423749328524996327035388417792688154433173133294892285295756457561276315648477233732e-1,
                0,
                0,
                0,
                0,
                0,
                0,
                .4640774434539039636406222168781981616534115643208114455689698789119941732444857047798e-2,
                .4704660282615136532130927218172390570903230981414159347904277946537920001824903276586e-1,
                .8620749948011488160369445167416002799205317397013619044391270706339561700281526529703e-1,
                -.2607983024682138093233254079066687623148682426317395111719299641390118652802949600035e-1,
                -.3858020174396621532493277639159499581333235076531298977820093139813399390137768850940e-1;

        _a.block<1,13>(13,0) <<  .2318046717429411567006043539613275607940758021709332569729352990777336390158311630529e-1,
                0,
                0,
                0,
                0,
                0,
                0,
                .3197856784116367067302124322582100058864027838197120089129330601737324659881765852593,
                .5933233331841898686063939886797828376866051205773280426848164018120869674204443797948,
                -1.937519548878479314706815782408229952008442222624773168771865465659822020582450444783,
                .1803950557030502357344063195737827904476240180662764468232042537858892203518134072359,
                -.4554014298857220726863505256926549022316460712353658688873150702827663762861750674926,
                2.158764106255762807077594619172645539322916635447781333204724468181634037726021280742;

        _a.block<1,14>(14,0) << .2624364325798105891527733985858552391723553030719144065844544880498188553839263944447e-1,
                0,
                0,
                0,
                0,
                0,
                0,
                .4863139423867266106526843913609225996253073727381961544415263239431571586043622332760e-1,
                .4274382538346478867636942429421724367591866585774144180215122660980822123988151132213e-1,
                -.4862259869465547771298976981868643277396586803130813159599600102115609499827986711663,
                .1326047194917652331781527125743684254490968718259563958293167893998110899691451568372,
                -.9402962152946515651634831658142934852383791641671387741034606371378082209616938685225e-1,
                .6993864679941022534190304512277131176659196396138275832136258135631963192299339871223,
                -.1197020013028860976492784934312243036670658451195397948726104511062042521592125912599e-1;

        _a.block<1,15>(15,0) << .5568066641536216461090823068917803436066365804361903532125349474551476120813558125830e-1,
                0,
                0,
                0,
                0,
                0,
                0,
                -.4324853319508358432896036654421685136736530810118924113940744870078036705505610668088,
                -.9979726994172038714656907882931844552238093285811791155499130927685987422432191170216,
                2.707893755718926115778725270396739994070337972517006747100005607751792006959604868323,
                -1.024823023512132929313567156576969954855232272749038347671818195935585095295127839150,
                1.334565206642246959252239602313589265188981560552694580059808406200559397799055652161,
                -2.587748998830690939658228913150922979184368065866213469477796089200252812362701917187,
                0.8992773696348355846430438306111181223414632598285854300924423251352733205187087732678e-1,
                1.497578446211167333777988534023066333042434967475357134513165331964695787890042760189;

        _a.block<1,16>(16,0) << -.8434891199686377639125188391985671318383858641413517143104162188088468627447515172982e-3,
                0,
                0,
                0,
                0,
                0,
                0,
                .7602144218856081893754106886111596435015500427480120290148318740899211421773423234728,
                1.769083927820959377467464871522349066447068428702073590698445112684989184432409492025,
                -4.499239797622297101452915424261016593995695456495268863455643396071539024609271033574,
                1.490558190212043468817221563278239942209691100326719140478588601720867838040211450448,
                -2.552203480132132516997563217309689292804518121743365818482497611667126218719069737195,
                4.795167551528575994217413424533259845001657006088189480440731104737960266616292993321,
                -.9161854401769482236671414092387917470686251714192236693920061138984202381209109248553e-1,
                -1.525735678746850818217653470352135651821164556169070505816135230784807058389577753184,
                .7371445601564892133467497107205798584829803038168267854389817508169123996459113657504;

        _a.block<1,17>(17,0) << .1017366974111576638766809656369828971944080018220332809259398740674738807023371082700,
                0,
                0,
                0,
                0,
                0,
                0,
                -1.696217553209432810711666838709742166182992092906177246174096517233561845662947862824,
                -3.825235846211624254528740857512255693551264719132875740261231165548583482101116676418,
                9.754768979885866648856431516333641627109105703674164986615824197909762854575668793816,
                -2.520767789227152291196336314591227486393143379933686189126240710041836742414125694941,
                5.472417145227780046950992000565734793413395536531652419585004300790370984185945495978,
                -9.781098113458736121002383874108051372067873053264954833376114258940736444388841687929,
                .3189152692455334369024560213486753019540464785641163242047782111839399471147176681561,
                3.447227036527756718156475010324322155277035924051392880570525223655410460762027138915,
                -.6051983612219277832241707671295607127814820499715293613761402732652780120810041653591,
                .3334525350307787459202631378414806560287636505658634784117511174230383993073398823363;

        _a.block<1,18>(18,0) << -.1012987737478284424676828232882617689682012456457322189102956361570156443805900941944,
                0,
                0,
                0,
                0,
                -.2409389328948775401304659380663043147167897928467308244359962659633933617326533285822e-1,
                -.6679880790275182076676283582867036095782150170801495251932447614617249253864579543857,
                1.600262798493100648047998296908183265688507618079976446601985464092263571149154964705,
                3.706958893826695766827011000213884379914407774639901049574259778345288538246990591819,
                -8.581755560147929325446798534254342948628755672447282004336563881429983605741487870996,
                .5607314974300953986559644699099897253584501767603091982484141468619493221310582281877e-1,
                -4.547761497422899514520768375507009011918601407646237921467449197008085790456674001879,
                9.255775439941294621826928846245618922061242300726600002589630404152665447900428712156,
                -.3450876657451631707159097079770789925142348071643902737346329921538351794816584861003,
                0,
                0,
                0,
                0;

        _a.block<1,19>(19,0) << .3826909723812638609001259641818040193828105314579492422836388985468479567237561247336e-1,
                0,
                0,
                0,
                0,
                .7786978965202527814624406274393101840018332461648638653990700950184871893714491273096,
                .4859454140913448249612202172501868752761599132465501266008866131088163955018926230543,
                1.814925350154666364151014269029611427420766367555858499108920245656959783343309816408,
                4.551165245704657956889158854062833952834232753889932986749613143631480116805870313264,
                -7.173770670344544101351160462586215092596352548535380880420409450623251883641801862305,
                -.3943009017000923237232456850787591816773705728833192412204243696911216045268772747196,
                -6.036544185898100312430357626685382432626027303329497026597513524312479466987506315664,
                7.338904299721887701527380004651998686389416058019429466200740313593568240326087171554,
                -.4143158595971836110248598960027762194900538872022960061452263646470675916118824501965,
                0,
                0,
                0,
                0,
                -.3732349451502749258108621577582478607301443393311959731632798508493352335121760204375;

        _a.block<1,20>(20,0) << .2162339046022045866878628785550588026780578552494608097931198882276791962912244674840e-1,
                0,
                0,
                0,
                0,
                .4611834700744369218866370212060318930941322187829670117414118166503940620998117275429,
                .1940797759547798743610542713744618433967649025379792966207862125676964319674160574624,
                .7041001229739959807963554405302474570280838416767002383409508232534658577705201658489,
                2.877431096792763528910415905652149398266490601780194388811216042455337979365709745445,
                0,
                -.4332742088749107411735902392606181444105337491234912425673655059805456011404518143074,
                -2.234178753588834452567105459024473991729105867012210449973082203886376638514123583334,
                .2235678086885984010238782832657956960650576194069632574873732156942360146780276407657,
                .1293532338308457711442786069651741293532338308457711442786069651741293532338308457711,
                0,
                0,
                0,
                0,
                .1418136968194278394808045812385429206355105705182818920178205766092934777719870449624,
                -1.085699633131323582531514699802817081967439754938101617737029931360398856861850276906;


        _b.setZero(_numStages,1);
        _b << .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                .1387145942588715882541801312803271702142521598590204181697361204933422401935856968980,
                .1892374781489234901583064041060123262381623469486258303271944256799821862794952728707,
                .9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1,
                .2774291885177431765083602625606543404285043197180408363394722409866844803871713937960,
                .1387145942588715882541801312803271702142521598590204181697361204933422401935856968980,
                .9461873907446174507915320205300616311908117347431291516359721283999109313974763643533e-1,
                .3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1,
                0,
                0,
                0;

        _b2.setZero(_numStages,1);
        _b2 << .3339829895931337572271945815422988633728883413227543303554098429202731077409488318421e-1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                .5024509803921568627450980392156862745098039215686274509803921568627450980392156862745e-1,
                -.1423859191318858946753152353981644782061337055184060977838998119673893661279423564924,
                .2126013199429258434998789109063801828540550730541648287733608913970804891935883227446,
                .3254854965632843133622967470840062095221514741629108993207015882688341071771214986692,
                .3312629399585921325051759834368530020703933747412008281573498964803312629399585921325,
                .1887845809230650005639203350759631573314744764356665687950807917985316096487827551997,
                0,
                0,
                0,
                0,
                .6159811094287144604404508847679200761569154839698701533406779753675999506388462962070e-1,
                -.9440109660594088037957791636147830082275023098021053120999763935859315478744850552959e-1,
                .3341117040855897708234682470384970584684876341854831047975628586614323631403861184369e-1;

        _c.setZero(_numStages,1);
        _c << 0.0,
                .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401,
                .1452518960316150517617548528770033320314511251329947060838468741983976455607179673401,
                .2178778440474225776426322793155049980471766876994920591257703112975964683410769510101,
                .5446946101185564441065806982887624951179417192487301478144257782439911708526923775252,
                .6536335321422677329278968379465149941415300630984761773773109338927894050232308530303,
                .2746594919905254008808021630247618520892150865127407293922085868737635475402543533498,
                .7735775201106609448405825008093973718589542913426807556412662673054607938029043386501,
                .5801831400829957086304368756070480288942157185070105667309497004790955953521782539876,
                .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383,
                .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092,
                .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908,
                .1174723380352676535744985130203309248171321557319478803362088220814723414805867429383,
                .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617,
                .3573842417596774518429245029795604640404982636367873040901247917361510345429002009092,
                .6426157582403225481570754970204395359595017363632126959098752082638489654570997990908,
                .8825276619647323464255014869796690751828678442680521196637911779185276585194132570617,
                1.0,
                .3510848126232741617357001972386587771203155818540433925049309664694280078895463510848,
                .6157407407407407407407407407407407407407407407407407407407407407407407407407407407407,
                1.0;


        //for( int i=0; i<_numStages-1; i++ )
        //    cout << "_a.block<1,20>("<<i<<",0)" << endl<< _a.block<1,20>(i,0).transpose() << endl << endl;


        _first_same_as_last = false;
        _adaptive_method = true;
        break;
    }
    case RK_12_10:
    {

        _name = "RK 12 and 10 th Order";
        _order = 12;
        _numStages = 25;
        _type = EXPLICIT;


        _c.setZero(_numStages,1);
        _c << 0.000000000000000000000000000000000000000000000000000000000000,
                0.200000000000000000000000000000000000000000000000000000000000,
                0.555555555555555555555555555555555555555555555555555555555556,
                0.833333333333333333333333333333333333333333333333333333333333,
                0.333333333333333333333333333333333333333333333333333333333333,
                1.00000000000000000000000000000000000000000000000000000000000,
                0.671835709170513812712245661002797570438953420568682550710222,
                0.288724941110620201935458488967024976908118598341806976469674,
                0.562500000000000000000000000000000000000000000000000000000000,
                0.833333333333333333333333333333333333333333333333333333333333,
                0.947695431179199287562380162101836721649589325892740646458322,
                0.0548112876863802643887753674810754475842153612931128785028369,
                0.0848880518607165350639838930162674302064148175640019542045934,
                0.265575603264642893098114059045616835297201264164077621448665,
                0.500000000000000000000000000000000000000000000000000000000000,
                0.734424396735357106901885940954383164702798735835922378551335,
                0.915111948139283464936016106983732569793585182435998045795407,
                0.947695431179199287562380162101836721649589325892740646458322,
                0.833333333333333333333333333333333333333333333333333333333333,
                0.288724941110620201935458488967024976908118598341806976469674,
                0.671835709170513812712245661002797570438953420568682550710222,
                0.333333333333333333333333333333333333333333333333333333333333,
                0.555555555555555555555555555555555555555555555555555555555556,
                0.200000000000000000000000000000000000000000000000000000000000,
                1.00000000000000000000000000000000000000000000000000000000000;

        _b.setZero(_numStages,1);
        _b << 0.0238095238095238095238095238095238095238095238095238095238095,
                0.0234375000000000000000000000000000000000000000000000000000000,
                0.0312500000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0416666666666666666666666666666666666666666666666666666666667,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0500000000000000000000000000000000000000000000000000000000000,
                0.0500000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.100000000000000000000000000000000000000000000000000000000000,
                0.0714285714285714285714285714285714285714285714285714285714286,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.138413023680782974005350203145033146748813640089941234591267,
                0.215872690604931311708935511140681138965472074195773051123019,
                0.243809523809523809523809523809523809523809523809523809523810,
                0.215872690604931311708935511140681138965472074195773051123019,
                0.138413023680782974005350203145033146748813640089941234591267,
                -0.0714285714285714285714285714285714285714285714285714285714286,
                -0.100000000000000000000000000000000000000000000000000000000000,
                -0.0500000000000000000000000000000000000000000000000000000000000,
                -0.0500000000000000000000000000000000000000000000000000000000000,
                -0.0416666666666666666666666666666666666666666666666666666666667,
                -0.0312500000000000000000000000000000000000000000000000000000000,
                -0.0234375000000000000000000000000000000000000000000000000000000,
                0.0238095238095238095238095238095238095238095238095238095238095;

        _b2.setZero(_numStages,1);
        _b2 << 1.47178724881110408452949550989023611293535315518571691939396,
                0.787500000000000000000000000000000000000000000000000000000000,
                0.421296296296296296296296296296296296296296296296296296296296,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.291666666666666666666666666666666666666666666666666666666667,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.348600717628329563206854421629657569274689947367847465753757,
                0.229499544768994849582890233710555447073823569666506700662510,
                5.79046485790481979159831978177003471098279506036722411333192,
                0.418587511856506868874073759426596207226461447604248151080016,
                0.307039880222474002649653817490106690389251482313213999386651,
                -4.68700905350603332214256344683853248065574415794742040470287,
                3.13571665593802262152038152399873856554395436199962915429076,
                1.40134829710965720817510506275620441055845017313930508348898,
                -5.52931101439499023629010306005764336421276055777658156400910,
                -0.853138235508063349309546894974784906188927508039552519557498,
                0.103575780373610140411804607167772795518293914458500175573749,
                -0.140474416950600941142546901202132534870665923700034957196546,
                -0.418587511856506868874073759426596207226461447604248151080016,
                -0.229499544768994849582890233710555447073823569666506700662510,
                -0.348600717628329563206854421629657569274689947367847465753757,
                -0.291666666666666666666666666666666666666666666666666666666667,
                -0.421296296296296296296296296296296296296296296296296296296296,
                -0.787500000000000000000000000000000000000000000000000000000000,
                0.0;


        _a.setZero(_numStages,_numStages);

        _a.block<1,1>(1,0) << 0.200000000000000000000000000000000000000000000000000000000000;

        _a.block<1,2>(2,0) << -0.216049382716049382716049382716049382716049382716049382716049,
                0.771604938271604938271604938271604938271604938271604938271605;

        _a.block<1,3>(3,0) << 0.208333333333333333333333333333333333333333333333333333333333,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.625000000000000000000000000000000000000000000000000000000000;

        _a.block<1,4>(4,0) <<  0.193333333333333333333333333333333333333333333333333333333333,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.220000000000000000000000000000000000000000000000000000000000,
                -0.0800000000000000000000000000000000000000000000000000000000000;

        _a.block<1,5>(5,0) << 0.100000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.400000000000000000000000000000000000000000000000000000000000,
                0.500000000000000000000000000000000000000000000000000000000000;

        _a.block<1,6>(6,0) << 0.103364471650010477570395435690481791543342708330349879244197,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.124053094528946761061581889237115328211074784955180298044074,
                0.483171167561032899288836480451962508724109257517289177302380,
                -0.0387530245694763252085681443767620580395733302341368038804290;

        _a.block<1,7>(7,0) << 0.124038261431833324081904585980175168140024670698633612292480,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.217050632197958486317846256953159942875916353757734167684657,
                0.0137455792075966759812907801835048190594443990939408530842918,
                -0.0661095317267682844455831341498149531672668252085016565917546;

        _a.block<1,8>(8,0) << 0.0914774894856882983144991846980432197088832099976660100090486,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                5   -0.00544348523717469689965754944144838611346156873847009178068318,
                0.0680716801688453518578515120895103863112751730758794372203952,
                0.408394315582641046727306852653894780093303185664924644551239;

        _a.block<1,9>(9,0) << 0.0890013652502551018954509355423841780143232697403434118692699,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.00499528226645532360197793408420692800405891149406814091955810,
                0.397918238819828997341739603001347156083435060931424970826304,
                0.427930210752576611068192608300897981558240730580396406312359,
                -0.0865117637557827005740277475955029103267246394128995965941585;

        _a.block<1,10>(10,0) << 0.0695087624134907543112693906409809822706021061685544615255758,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.129146941900176461970759579482746551122871751501482634045487,
                1.53073638102311295076342566143214939031177504112433874313011,
                0.577874761129140052546751349454576715334892100418571882718036,
                -0.951294772321088980532340837388859453930924498799228648050949,
                -0.408276642965631951497484981519757463459627174520978426909934;

        _a.block<1,11>(11,0) << 0.0444861403295135866269453507092463581620165501018684152933313,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.00380476867056961731984232686574547203016331563626856065717964,
                0.0106955064029624200721262602809059154469206077644957399593972,
                0.0209616244499904333296674205928919920806734650660039898074652,
                -0.0233146023259321786648561431551978077665337818756053603898847,
                0.00263265981064536974369934736325334761174975280887405725010964,
                0.00315472768977025060103545855572111407955208306374459723959783;

        _a.block<1,12>(12,0) << 0.0194588815119755475588801096525317761242073762016273186231215,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0000678512949171812509306121653452367476194364781259165332321534,
                -0.0000429795859049273623271005330230162343568863387724883603675550,
                0.0000176358982260285155407485928953302139937553442829975734148981,
                0.0653866627415027051009595231385181033549511358787382098351924;

        _a.block<1,13>(13,0) << 0.206836835664277105916828174798272361078909196043446411598231,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0166796067104156472828045866664696450306326505094792505215514,
                -0.00879501563200710214457024178249986591130234990219959208704979,
                0.00346675455362463910824462315246379209427513654098596403637231,
                -0.861264460105717678161432562258351242030270498966891201799225,
                0.908651882074050281096239478469262145034957129939256789178785;

        _a.block<1,14>(14,0) << 0.0203926084654484010091511314676925686038504449562413004562382,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0869469392016685948675400555583947505833954460930940959577347,
                -0.0191649630410149842286436611791405053287170076602337673587681,
                0.00655629159493663287364871573244244516034828755253746024098838,
                0.0987476128127434780903798528674033899738924968006632201445462,
                0.00535364695524996055083260173615567408717110247274021056118319,
                0.301167864010967916837091303817051676920059229784957479998077;

        _a.block<1,15>(15,0) << 0.228410433917778099547115412893004398779136994596948545722283,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.498707400793025250635016567442511512138603770959682292383042,
                0.134841168335724478552596703792570104791700727205981058201689,
                -0.0387458244055834158439904226924029230935161059142806805674360,
                -1.27473257473474844240388430824908952380979292713250350199641,
                1.43916364462877165201184452437038081875299303577911839630524,
                -0.214007467967990254219503540827349569639028092344812795499026,
                0.958202417754430239892724139109781371059908874605153648768037;

        _a.block<1,16>(16,0) << 2.00222477655974203614249646012506747121440306225711721209798,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                2.06701809961524912091954656438138595825411859673341600679555,
                0.623978136086139541957471279831494466155292316167021080663140,
                -0.0462283685500311430283203554129062069391947101880112723185773,
                -8.84973288362649614860075246727118949286604835457092701094630,
                7.74257707850855976227437225791835589560188590785037197433615,
                -0.588358519250869210993353314127711745644125882130941202896436,
                -1.10683733362380649395704708016953056176195769617014899442903,
                -0.929529037579203999778397238291233214220788057511899747507074;

        _a.block<1,17>(17,0) << 3.13789533412073442934451608989888796808161259330322100268310,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.129146941900176461970759579482746551122871751501482634045487,
                1.53073638102311295076342566143214939031177504112433874313011,
                0.577874761129140052546751349454576715334892100418571882718036,
                5.42088263055126683050056840891857421941300558851862156403363,
                0.231546926034829304872663800877643660904880180835945693836936,
                0.0759292995578913560162301311785251873561801342333194895292058,
                -12.3729973380186513287414553402595806591349822617535905976253,
                9.85455883464769543935957209317369202080367765721777101906955,
                0.0859111431370436529579357709052367772889980495122329601159540,
                -5.65242752862643921117182090081762761180392602644189218673969,
                -1.94300935242819610883833776782364287728724899124166920477873,
                -0.128352601849404542018428714319344620742146491335612353559923;

        _a.block<1,18>(18,0) << 1.38360054432196014878538118298167716825163268489922519995564,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.00499528226645532360197793408420692800405891149406814091955810,
                0.397918238819828997341739603001347156083435060931424970826304,
                0.427930210752576611068192608300897981558240730580396406312359,
                -1.30299107424475770916551439123047573342071475998399645982146,
                0.661292278669377029097112528107513072734573412294008071500699,
                -0.144559774306954349765969393688703463900585822441545655530145,
                -6.96576034731798203467853867461083919356792248105919255460819,
                6.65808543235991748353408295542210450632193197576935120716437,
                -1.66997375108841486404695805725510845049807969199236227575796,
                2.06413702318035263832289040301832647130604651223986452170089,
                -0.674743962644306471862958129570837723192079875998405058648892,
                -0.00115618834794939500490703608435907610059605754935305582045729,
                -0.00544057908677007389319819914241631024660726585015012485938593;

        _a.block<1,19>(19,0) << 0.951236297048287669474637975894973552166903378983475425758226,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.217050632197958486317846256953159942875916353757734167684657,
                0.0137455792075966759812907801835048190594443990939408530842918,
                -0.0661095317267682844455831341498149531672668252085016565917546,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.152281696736414447136604697040747131921486432699422112099617,
                -0.337741018357599840802300793133998004354643424457539667670080,
                -0.0192825981633995781534949199286824400469353110630787982121133,
                -3.68259269696866809932409015535499603576312120746888880201882,
                3.16197870406982063541533528419683854018352080342887002331312,
                -0.370462522106885290716991856022051125477943482284080569177386,
                -0.0514974200365440434996434456698127984941168616474316871020314,
                -0.000829625532120152946787043541792848416659382675202720677536554,
                0.00000279801041419278598986586589070027583961355402640879503213503,
                0.0418603916412360287969841020776788461794119440689356178942252,
                0.279084255090877355915660874555379649966282167560126269290222;

        _a.block<1,20>(20,0) << 0.103364471650010477570395435690481791543342708330349879244197,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.124053094528946761061581889237115328211074784955180298044074,
                0.483171167561032899288836480451962508724109257517289177302380,
                -0.0387530245694763252085681443767620580395733302341368038804290,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.438313820361122420391059788940960176420682836652600698580091,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.218636633721676647685111485017151199362509373698288330593486,
                -0.0312334764394719229981634995206440349766174759626578122323015,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0312334764394719229981634995206440349766174759626578122323015,
                0.218636633721676647685111485017151199362509373698288330593486,
                0.438313820361122420391059788940960176420682836652600698580091;

        _a.block<1,21>(21,0) << 0.193333333333333333333333333333333333333333333333333333333333,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.220000000000000000000000000000000000000000000000000000000000,
                -0.0800000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.0984256130499315928152900286856048243348202521491288575952143,
                -0.196410889223054653446526504390100417677539095340135532418849,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.436457930493068729391826122587949137609670676712525034763317,
                0.0652613721675721098560370939805555698350543810708414716730270,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.0652613721675721098560370939805555698350543810708414716730270,
                -0.436457930493068729391826122587949137609670676712525034763317,
                0.196410889223054653446526504390100417677539095340135532418849,
                -0.0984256130499315928152900286856048243348202521491288575952143;

        _a.block<1,22>(22,0) << -0.216049382716049382716049382716049382716049382716049382716049,
                0.771604938271604938271604938271604938271604938271604938271605,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.666666666666666666666666666666666666666666666666666666666667,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.390696469295978451446999802258495981249099665294395945559163,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.390696469295978451446999802258495981249099665294395945559163,
                0.666666666666666666666666666666666666666666666666666666666667;

        _a.block<1,23>(23,0) << 0.200000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                -0.164609053497942386831275720164609053497942386831275720164609,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.164609053497942386831275720164609053497942386831275720164609;

        _a.block<1,24>(24,0) << 1.47178724881110408452949550989023611293535315518571691939396,
                0.787500000000000000000000000000000000000000000000000000000000,
                0.421296296296296296296296296296296296296296296296296296296296,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.291666666666666666666666666666666666666666666666666666666667,
                0.000000000000000000000000000000000000000000000000000000000000,
                0.348600717628329563206854421629657569274689947367847465753757,
                0.229499544768994849582890233710555447073823569666506700662510,
                5.79046485790481979159831978177003471098279506036722411333192,
                0.418587511856506868874073759426596207226461447604248151080016,
                0.307039880222474002649653817490106690389251482313213999386651,
                -4.68700905350603332214256344683853248065574415794742040470287,
                3.13571665593802262152038152399873856554395436199962915429076,
                1.40134829710965720817510506275620441055845017313930508348898,
                -5.52931101439499023629010306005764336421276055777658156400910,
                -0.853138235508063349309546894974784906188927508039552519557498,
                0.103575780373610140411804607167772795518293914458500175573749,
                -0.140474416950600941142546901202132534870665923700034957196546,
                -0.418587511856506868874073759426596207226461447604248151080016,
                -0.229499544768994849582890233710555447073823569666506700662510,
                -0.348600717628329563206854421629657569274689947367847465753757,
                -0.291666666666666666666666666666666666666666666666666666666667,
                -0.421296296296296296296296296296296296296296296296296296296296,
                -0.787500000000000000000000000000000000000000000000000000000000;

        _first_same_as_last = false;
        _adaptive_method = true;
        break;
    }
    case CONSTRAINT_BASED:
    {
        _name = "Constraint Based";
        _order = 0;
        _numStages = 0;
        _type = NONE;
        _adaptive_method = false;
        _first_same_as_last = false;

        _a.setZero(0,0);
        _b.setZero(0);
        _c.setZero(0);
        break;
    }

    }

    c_duplicate = std::vector<unsigned int>(_numStages,0);

    // calculation of duplicate indicies so interpolation can be better for some methods
    for( int i = _numStages-1; i>=0; i-- )
    {
        if( (int)c_duplicate[i] < i )
        {
            c_duplicate[i] = i;

        } else if( (int)c_duplicate[i] > i )
        {
            continue;
        }

        for( int j=0; j < i; j++ )
        {
            //            std::cout << "C(i) " << _c(i) << std::endl;
            //            std::cout << "C(j) " << _c(j) << std::endl;
            //            std::cout << "abs diff: " << std::abs( _c(i) - _c(j) ) << std::endl<<std::endl;
            if( std::abs( _c(i) - _c(j) ) < 1e-6 )
            {
                c_duplicate[j] = i;
            }
        }
    }


    //cout << "Method Enum: " << (int) method_type << "\tPtr Val: " <<   (void*) this << endl;


}





//    case DORMAN_PRINCE_12TH:
//    {
//        /*
//         * P.J. Prince & J.R. Dorman (1981)
//         * High order embedded Runge-Kutta formulae. J.Comp. Appl. Math., Vol. 7. p.67-75.
//         *
//        */

//        _name = "Dorman Prince 12th Order";
//        _order = 12;
//        _numStages = 17;
//        _type = EXPLICIT;

//        _a.setZero(_numStages,_numStages);
//        _a.block<1,1>(1,0) << 2.0e-4;

//        _a.block<1,2>(2,0) << 2.66666666666666666666666666667e-4,
//                5.33333333333333333333333333333e-4;

//        _a.block<1,3>(3,0) << 2.91666666666666666666666666667e-3,
//                -4.16666666666666666666666666667e-3,
//                6.25e-3;

//        _a.block<1,4>(4,0) << 1.64609053497942386831275720165e-3,
//                0.0e0,
//                5.48696844993141289437585733882e-3,
//                1.75582990397805212620027434842e-3;

//        _a.block<1,5>(5,0) << 1.9456e-3,
//                0.0e0,
//                7.15174603174603174603174603175e-3,
//                2.91271111111111111111111111111e-3,
//                7.89942857142857142857142857143e-4;

//        _a.block<1,6>(6,0) << 5.6640625e-4,
//                0.0e0,
//                8.80973048941798941798941798942e-4,
//                -4.36921296296296296296296296296e-4,
//                3.39006696428571428571428571429e-4,
//                -9.94646990740740740740740740741e-5;

//        _a.block<1,7>(7,0) << 3.08333333333333333333333333333e-3,
//                0.0e0,
//                0.0e0,
//                1.77777777777777777777777777778e-3,
//                2.7e-3,
//                1.57828282828282828282828282828e-3,
//                1.08606060606060606060606060606e-2;

//        _a.block<1,8>(8,0) << 3.65183937480112971375119150338e-3,
//                0.0e0,
//                3.96517171407234306617557289807e-3,
//                3.19725826293062822350093426091e-3,
//                8.22146730685543536968701883401e-3,
//                -1.31309269595723798362013884863e-3,
//                9.77158696806486781562609494147e-3,
//                3.75576906923283379487932641079e-3;

//        _a.block<1,9>(9,0) << 3.70724106871850081019565530521e-3,
//                0.0e0,
//                5.08204585455528598076108163479e-3,
//                1.17470800217541204473569104943e-3,
//                -2.11476299151269914996229766362e-2,
//                6.01046369810788081222573525136e-2,
//                2.01057347685061881846748708777e-2,
//                -2.83507501229335808430366774368e-2,
//                1.48795689185819327555905582479e-2;

//        _a.block<1,10>(10,0) << 3.51253765607334415311308293052e-2,
//                0.0e0,
//                -8.61574919513847910340576078545e-3,
//                -5.79144805100791652167632252471e-3,
//                1.94555482378261584239438810411e0,
//                -3.43512386745651359636787167574e0,
//                -1.09307011074752217583892572001e-1,
//                2.3496383118995166394320161088e0,
//                -7.56009408687022978027190729778e-1,
//                1.09528972221569264246502018618e-1;

//        _a.block<1,11>(11,0) << 2.05277925374824966509720571672e-2,
//                0.0e0,
//                -7.28644676448017991778247943149e-3,
//                -2.11535560796184024069259562549e-3,
//                9.27580796872352224256768033235e-1,
//                -1.65228248442573667907302673325e0,
//                -2.10795630056865698191914366913e-2,
//                1.20653643262078715447708832536e0,
//                -4.13714477001066141324662463645e-1,
//                9.07987398280965375956795739516e-2,
//                5.35555260053398504916870658215e-3;

//        _a.block<1,12>(12,0) << -1.43240788755455150458921091632e-1,
//                0.0e0,
//                1.25287037730918172778464480231e-2,
//                6.82601916396982712868112411737e-3,
//                -4.79955539557438726550216254291e0,
//                5.69862504395194143379169794156e0,
//                7.55343036952364522249444028716e-1,
//                -1.27554878582810837175400796542e-1,
//                -1.96059260511173843289133255423e0,
//                9.18560905663526240976234285341e-1,
//                -2.38800855052844310534827013402e-1,
//                1.59110813572342155138740170963e-1;

//        _a.block<1,13>(13,0) << 8.04501920552048948697230778134e-1,
//                0.0e0,
//                -1.66585270670112451778516268261e-2,
//                -2.1415834042629734811731437191e-2,
//                1.68272359289624658702009353564e1,
//                -1.11728353571760979267882984241e1,
//                -3.37715929722632374148856475521e0,
//                -1.52433266553608456461817682939e1,
//                1.71798357382154165620247684026e1,
//                -5.43771923982399464535413738556e0,
//                1.38786716183646557551256778839e0,
//                -5.92582773265281165347677029181e-1,
//                2.96038731712973527961592794552e-2;

//        _a.block<1,14>(14,0) << -9.13296766697358082096250482648e-1,
//                0.0e0,
//                2.41127257578051783924489946102e-3,
//                1.76581226938617419820698839226e-2,
//                -1.48516497797203838246128557088e1,
//                2.15897086700457560030782161561e0,
//                3.99791558311787990115282754337e0,
//                2.84341518002322318984542514988e1,
//                -2.52593643549415984378843352235e1,
//                7.7338785423622373655340014114e0,
//                -1.8913028948478674610382580129e0,
//                1.00148450702247178036685959248e0,
//                4.64119959910905190510518247052e-3,
//                1.12187550221489570339750499063e-2;

//        _a.block<1,15>(15,0) << -2.75196297205593938206065227039e-1,
//                0.0e0,
//                3.66118887791549201342293285553e-2,
//                9.7895196882315626246509967162e-3,
//                -1.2293062345886210304214726509e1,
//                1.42072264539379026942929665966e1,
//                1.58664769067895368322481964272e0,
//                2.45777353275959454390324346975e0,
//                -8.93519369440327190552259086374e0,
//                4.37367273161340694839327077512e0,
//                -1.83471817654494916304344410264e0,
//                1.15920852890614912078083198373e0,
//                -1.72902531653839221518003422953e-2,
//                1.93259779044607666727649875324e-2,
//                5.20444293755499311184926401526e-3;

//        _a.block<1,16>(16,0) << 1.30763918474040575879994562983e0,
//                0.0e0,
//                1.73641091897458418670879991296e-2,
//                -1.8544456454265795024362115588e-2,
//                1.48115220328677268968478356223e1,
//                9.38317630848247090787922177126e0,
//                -5.2284261999445422541474024553e0,
//                -4.89512805258476508040093482743e1,
//                3.82970960343379225625583875836e1,
//                -1.05873813369759797091619037505e1,
//                2.43323043762262763585119618787e0,
//                -1.04534060425754442848652456513e0,
//                7.17732095086725945198184857508e-2,
//                2.16221097080827826905505320027e-3,
//                7.00959575960251423699282781988e-3,
//                0.0e0;



//        _b.setZero(_numStages,1);
//        _b << 1.21278685171854149768890395495e-2,
//              0.0e0,
//              0.0e0,
//              0.0e0,
//              0.0e0,
//              0.0e0,
//              8.62974625156887444363792274411e-2,
//              2.52546958118714719432343449316e-1,
//              -1.97418679932682303358307954886e-1,
//              2.03186919078972590809261561009e-1,
//              -2.07758080777149166121933554691e-2,
//              1.09678048745020136250111237823e-1,
//              3.80651325264665057344878719105e-2,
//              1.16340688043242296440927709215e-2,
//              4.65802970402487868693615238455e-3,
//              0.0e0,
//              0.0e0;

//        _b2.setZero(_numStages,1);
//        _b2 << 1.70087019070069917527544646189e-2,
//               0.0e0,
//               0.0e0,
//               0.0e0,
//               0.0e0,
//               0.0e0,
//               7.22593359308314069488600038463e-2,
//               3.72026177326753045388210502067e-1,
//               -4.01821145009303521439340233863e-1,
//               3.35455068301351666696584034896e-1,
//               -1.31306501075331808430281840783e-1,
//               1.89431906616048652722659836455e-1,
//               2.68408020400290479053691655806e-2,
//               1.63056656059179238935180933102e-2,
//               3.79998835669659456166597387323e-3,
//               0.0e0,
//               0.0e0;

// // ///// B Prime
// //        _b.setZero(_numStages,1);
// //        _b << 1.21278685171854149768890395495e-2,
// //              0.0e0,
// //              0.0e0,
// //              0.0e0,
// //              0.0e0,
// //              0.0e0,
// //              9.08394342270407836172412920433e-2,
// //              3.15683697648393399290429311645e-1,
// //              -2.63224906576909737811077273181e-1,
// //              3.04780378618458886213892341513e-1,
// //              -4.15516161554298332243867109382e-2,
// //              2.46775609676295306562750285101e-1,
// //              1.52260530105866022937951487642e-1,
// //              8.14384816302696075086493964505e-2,
// //              8.50257119389081128008018326881e-2,
// //              -9.15518963007796287314100251351e-3,
// //              2.5e-2;


// //        _b2.setZero(_numStages,1);
// //        _b2 << 1.70087019070069917527544646189e-2,
// //              0.0e0,
// //              0.0e0,
// //              0.0e0,
// //              0.0e0,
// //              0.0e0,
// //              7.60624588745593757356421093119e-2,
// //              4.65032721658441306735263127583e-1,
// //              -5.35761526679071361919120311817e-1,
// //              5.03182602452027500044876052344e-1,
// //              -2.62613002150663616860563681567e-1,
// //              4.26221789886109468625984632024e-1,
// //              1.07363208160116191621476662322e-1,
// //              1.14139659241425467254626653171e-1,
// //              6.93633866500486770090602920091e-2,
// //              2.0e-2,
// //              0.0e0;

//        _c.setZero(_numStages,1);
//        _c << 0.0,
//                2.0e-2,
//                4.0e-2,
//                1.0e-1,
//                1.33333333333333333333333333333e-1,
//                1.6e-1,
//                5.0e-2,
//                2.0e-1,
//                2.5e-1,
//                3.33333333333333333333333333333e-1,
//                5.0e-1,
//                5.55555555555555555555555555556e-1,
//                7.5e-1,
//                8.57142857142857142857142857143e-1,
//                9.45216222272014340129957427739e-1,
//                1.0,
//                1.0;

//        _first_same_as_last = false;
//         _adaptive_method = true;
//        break;
//    }
