#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <eigen3/Eigen/Dense>

#include <vector>
#include <iostream>
#include <iomanip>
using std::cout;
using std::setprecision;
using std::endl;

#include <cmath>

//#include "Utilities/NumericalDifferentiation.h"
#include "PseudoInverse.h"

//#include "Utilities/Timer.h"

struct _nr_default{};

/**
 * @brief The jacobian_update_method enum Describes the Jacobian update method
 *  EVALUATE -- Calculates the jacobian either numerically or analitically each step
 *  APPROXIMATE_BROYDEN_1 -- Uses the Broyden "Good" method to perform a rank 1 update on the inverse jacobian
 *  APPROXIMATE_BROYDEN_2 -- Uses Broyden's other method to perform a rank 1 update on the inverse jacobian
 *
 */
enum jacobian_update_method {EVALUATE, APPROXIMATE_BROYDEN_1, APPROXIMATE_BROYDEN_2};

/**
 * @brief The search_method enum Describes the search method to use:
 * - NEWTON_RAPHSON: Standard least-squares iterative root finding.  Can overshoot and blow up, does not work well with problems with negative curvature.
 * - NEWTON_RAPHSON_QUADRATIC_CORRECT: Does a simplifed line search using a least-squres parabolic fit to the error along the NEWTON_RAPHSON Direction.  Good for problems with negative curvature.
 * - LEVENBERG_MARQUARDT:  Damped least squres with a coefficeint that changes depending on the history.  Good from problems that are hard to solve.  Slower than NEWTON_RAPHSON.
 * - LEVENBERG_MARQUARDT_ADAPTIVE: Similar to LEVENBERG_MARQUARDT but a little bit faster.  It has fewer re-calcuations.
 * - LEVENBERG_MARQUARDT_ANTICIPATIVE: Uses a curvature estimat to choose the damping parameter.  Sometimes works well. Sometimes works poorly.
 *
 * Over all:  LEVENBERG_MARQUARDT_ADAPTIVE is the best at finding a solution robustly. NEWTON_RAPHSON_QUADRATIC_CORRECT is the best if the error function is well defined and behavied.
 */
enum search_method {NEWTON_RAPHSON, NEWTON_RAPHSON_PLUS_BRENTS, NEWTON_RAPHSON_QUADRATIC_CORRECT, LEVENBERG_MARQUARDT, LEVENBERG_MARQUARDT_ADAPTIVE, LEVENBERG_MARQUARDT_ANTICIPATIVE};


/**
 * @brief The NewtonRaphson class Performs function minimization using a newton decent method
 * This class can be used for least squares fitting of linear or nonlinear problems or for root finding of vector functions.
 * If the calculation of the jacobian is expensive, because it is done numerically or is very involved, the APPROXIMATE_BROYDEN_1
 * or APPROXIMATE_BROYDEN_2 update methods can be used to increase the speed of the program.  This funciton requires a good initial
 * guess of the solution if the funciton could have more than one maximum or minimum.
 *
 * Example:
 * \code
 *  #include "NewtonRaphson.h"
 *  #include <iostream>
 *  #include <eigen3/Eigen/Core>
 *
 *  // Define the function to be fit A+B*t+C*sin(6*t)+D^3
 *  Eigen::VectorXd errFcn(const Eigen::VectorXd& state)
 *  {
 *      Eigen::VectorXd output;
 *      output.setZero(60,1);
 *
 *      for(int i=0;i<60;i++)
 *      {
 *          double t = 3.141592654 * i/30.0;
 *
 *          output(i) = (5        + 3*t        + 4*sin(6*t) + pow(6,3))
 *                    - (state(0) + state(1)*t + state(2)*sin(6*t) + pow(state(3),3) );
 *      }
 *
 *      return output;
 *  }
 *
 *  // Define the analitical jacobian
 *  Eigen::MatrixXd jacobFcn(const Eigen::VectorXd& state)
 *  {
 *      Eigen::MatrixXd output;
 *      output.setZero(60,4);
 *
 *      for(int i=0;i<60; i++)
 *      {
 *          double t = 3.141592654 * i/30.0;
 *
 *          output.block<1,4>(i,0) << -1, -t, -sin(6*t), -3*pow(state(3),2);
 *
 *      }
 *
 *      return output;
 *  }
 *
 *  int main()
 *  {
 *      Eigen::VectorXd IC(4,1);
 *      IC(0) = 4.2;
 *      IC(1) = 2.2;
 *      IC(2) = 4.6;
 *      IC(3) = 8.9;
 *
 *      // Solve with analitical jacobian evaluation at each step
 *      NewtonRaphson<> leastSquaresSolver(errFcn,jacobFcn,EVALUATE,-1,-1);
 *
 *      Eigen::VectorXd soln(4,1);
 *
 *      soln = leastSquaresSolver.solve(IC,1e-10,100,true,NEWTON_RAPHSON);
 *      cout << "NEWTON_RAPHSON, Analitical Jacobian. Soln: " << soln.transpose() << std::endl << "----------" << std::endl << std::endl;
 *
 *      soln = leastSquaresSolver.solve(IC,1e-10,100,true,NEWTON_RAPHSON_QUADRATIC_CORRECT);
 *      cout << "NEWTON_RAPHSON_QUADRATIC_CORRECT, Analitical Jacobian. Soln: " << soln.transpose() << std::endl << "----------" << std::endl << std::endl;
 *
 *      soln = leastSquaresSolver.solve(IC,1e-10,100,true,LEVENBERG_MARQUARDT);
 *      cout << "LEVENBERG_MARQUARDT, Analitical Jacobian. Soln: " << soln.transpose() << std::endl << "----------" << std::endl << std::endl;
 *
 *      soln = leastSquaresSolver.solve(IC,1e-10,100,true,LEVENBERG_MARQUARDT_ADAPTIVE);
 *      cout << "LEVENBERG_MARQUARDT_ADAPTIVE, Analitical Jacobian. Soln: " << soln.transpose() << std::endl << "----------" << std::endl << std::endl;
 *
 *      soln = leastSquaresSolver.solve(IC,1e-10,100,true,LEVENBERG_MARQUARDT_ANTICIPATIVE);
 *      cout << "LEVENBERG_MARQUARDT_ANTICIPATIVE, Analitical Jacobian. Soln: " << soln.transpose() << std::endl << "----------" << std::endl << std::endl;
 *
 *
 *      // Solve with numerical jacobian evaluation at each step
 *      leastSquaresSolver  = NewtonRaphson<>(errFcn,0,EVALUATE,-1,-1);
 *      soln = leastSquaresSolver.solve(IC,1e-6, NEWTON_RAPHSON);
 *      cout << "Numerical Jacobian Evaluated Every Time \n\tSoln: " << soln.transpose() << std::endl;
 *
 *
 *      // Solve with numerical jacobian evaluation only once with Broyden's Good method used to update the inverse jacobian
 *      leastSquaresSolver.changeUpdateMethod(APPROXIMATE_BROYDEN_1);
 *      soln = leastSquaresSolver.solve(IC,1e-6, NEWTON_RAPHSON);
 *      cout << "Numerical Jacobian Evaluated Once (Broyden's Method 1) \n\tSoln: " << soln.transpose()  << std::endl;
 *
 *
 *      // Solve with numerical jacobian evaluation only once with Broyden's other method used to update the inverse jacobian
 *      leastSquaresSolver.changeUpdateMethod(APPROXIMATE_BROYDEN_2);
 *      soln = leastSquaresSolver.solve(IC,1e-6, NEWTON_RAPHSON);
 *      cout << "Numerical Jacobian Evaluated Once (Broyden's Method 2) \n\tSoln: " << soln.transpose() << std::endl;
 *
 *      return 0;
 *  }
 *  \endcode
 */

template< class T = _nr_default >
class NewtonRaphson
{
protected:

    Eigen::VectorXd (*errorFcnPtr)(const Eigen::VectorXd&); /**< The error function pointer **/
    Eigen::MatrixXd (*errorJacobianFcnPtr)(const Eigen::VectorXd&);/**< The jacobian function pointer **/
    void (*errorFcnPtr_2)(Eigen::VectorXd&, Eigen::VectorXd&); /**< The error function pointer **/
    void (*errorJacobianFcnPtr_2)(Eigen::VectorXd&, Eigen::MatrixXd&);/**< The jacobian function pointer **/

    T* classObjectPointer;
    void (T::*class_errorFcnPtr)(Eigen::VectorXd&, Eigen::VectorXd&); /**< The error function pointer **/
    void (T::*class_errorJacobianFcnPtr)(Eigen::VectorXd&, Eigen::MatrixXd&);/**< The jacobian function pointer **/


    void returnError(Eigen::VectorXd& state, Eigen::VectorXd & error_out );
    void calculateNumericalJacobian();
    Eigen::VectorXd numericalDeritiveHelper( const Eigen::VectorXd & );


    jacobian_update_method update_method; /**< The desired update function **/
    search_method method;

    std::vector<unsigned int> HtoJ_index;
    std::vector<double> H_col_scale;
    Eigen::MatrixXd H; /**< Teh scaling Matrix **/
    Eigen::MatrixXd J; /**< The current jacobian **/
    Eigen::MatrixXd J_inv; /**< The current inverse jacobian **/
    //Eigen::ColPivHouseholderQR<Eigen::MatrixXd> J_QR;
    //Eigen::HouseholderQR<Eigen::MatrixXd> J_QR;
    Eigen::LDLT<Eigen::MatrixXd> J_QR;
    //Eigen::LLT<Eigen::MatrixXd> J_QR;
    //Eigen::ConjugateGradient<Eigen::MatrixXd,Eigen::Lower|Eigen::Upper> J_QR;
    bool updateJacobianDecomposition;


    double modelReductionThreshold;
    double delta_error_prediction;
    double sqrtEps;
    double Lambda;
    double LambdaChangeRange;
    unsigned int lambdaRangeCount;

    Eigen::VectorXd states; /**< The current state vector **/
    Eigen::VectorXd states_last; /**< The current state vector **/
    Eigen::VectorXd delta_States; /**< The last change in states **/
    Eigen::VectorXd delta_States_last;

    Eigen::VectorXd error_last; /**< The previous step's error vector **/
    Eigen::VectorXd error_this; /**< The current error vector **/
    Eigen::VectorXd delta_error; /**< The current change in error **/

    double rmsError_this;
    double rmsError_last;

    int updateJacobianEvery;
    unsigned int iteration;

    virtual void updateJacobian(int count = 0);
    virtual void calculateDeltaStates();
    virtual void updateError();
    virtual void removeNullSpace();
    //virtual void inforceConstraints();

    unsigned int countErrorEval;
    unsigned int countJacobEval;

    /**
     * @brief NewtonRaphson Constructor for inherited classes only.
     * @param update_method The method desired to update the jacobian. \see jacobian_update_method
     * @param updateJacobianEvery Specifies the frequency the jacobian shoudl be recalculated for the Broyden methods.
     * @parblock
     *  - -1: Specifies the default value for the jacobian update method.
     *  -  0: Specifies the jacobian sould be updated every iteration.
     *  -  N>0: Specifies the jacobian should be explicitly calculated every N steps.  This is
     *          benifficial for the Broyden jacobian update methods because it prevents the jacobian
     *          error from growing too large.
     *  .
     * @endparblock
     */
    NewtonRaphson(jacobian_update_method update_method, int updateJacobianEvery);

public:
    /**
     * @brief NewtonRaphson The constructor
     * @param errorFcnPtr   A pointer to the function that calculates the current error (in a least squares problem sense)
     * @param errorJacobianFcnPtr A pointer to the function that provides the jacobian mapping between the error and the current states.  If set to zero a numerical approximation will be used. \see CalculateJacobian
     * @param update_method The method desired to update the jacobian. \see jacobian_update_method
     * @param updateJacobianEvery Specifies the frequency the jacobian shoudl be recalculated for the Broyden methods.
     * @parblock
     *  - -1: Specifies the default value for the jacobian update method.
     *  -  0: Specifies the jacobian sould be updated every iteration.
     *  -  N>0: Specifies the jacobian should be explicitly calculated every N steps.  This is
     *          benifficial for the Broyden jacobian update methods because it prevents the jacobian
     *          error from growing too large.
     *  .
     * @endparblock
     *
     */
    NewtonRaphson(  Eigen::VectorXd (*errorFcnPtr)(const Eigen::VectorXd&), Eigen::MatrixXd (*errorJacobianFcnPtr)( const Eigen::VectorXd& ) = 0, jacobian_update_method update_method = EVALUATE, int updateJacobianEvery = -1, double modelReductionThreshold = -1 );

    /**
     * @brief NewtonRaphson The constructor
     * @param errorFcnPtr   A pointer to the function that calculates the current error (in a least squares problem sense)
     * @param errorJacobianFcnPtr A pointer to the function that provides the jacobian mapping between the error and the current states.  This cannot be a null pointer.
     * @param update_method The method desired to update the jacobian. \see jacobian_update_method
     * @param updateJacobianEvery Specifies the frequency the jacobian shoudl be recalculated for the Broyden methods.
     * @parblock
     *  - -1: Specifies the default value for the jacobian update method.
     *  -  0: Specifies the jacobian sould be updated every iteration.
     *  -  N>0: Specifies the jacobian should be explicitly calculated every N steps.  This is
     *          benifficial for the Broyden jacobian update methods because it prevents the jacobian
     *          error from growing too large.
     *  .
     * @endparblock
     */
    NewtonRaphson(  void (*errorFcnPtr)(Eigen::VectorXd&, Eigen::VectorXd&), void (*errorJacobianFcnPtr)( Eigen::VectorXd&, Eigen::MatrixXd& ), jacobian_update_method update_method = EVALUATE, int updateJacobianEvery = -1, double modelReductionThreshold = -1 );

    /**
     * @brief NewtonRaphson_C The constructor for a Newton-Raphson minimizer with class-based functions.
     * @param classPointer A pointer to the class object that owns the functions.
     * @param errorFcnPtr   A pointer to the function that calculates the current error (in a least squares problem sense)
     * @param errorJacobianFcnPtr A pointer to the function that provides the jacobian mapping between the error and the current states.  This cannot be a null pointer.
     * @param update_method The method desired to update the jacobian. \see jacobian_update_method
     * @param updateJacobianEvery Specifies the frequency the jacobian shoudl be recalculated for the Broyden methods.
     * @parblock
     *  - -1: Specifies the default value for the jacobian update method.
     *  -  0: Specifies the jacobian sould be updated every iteration.
     *  -  N>0: Specifies the jacobian should be explicitly calculated every N steps.  This is
     *          benifficial for the Broyden jacobian update methods because it prevents the jacobian
     *          error from growing too large.
     *  .
     * @endparblock
     */
    NewtonRaphson( T* classPointer, void (T::*errorFcnPtr)( Eigen::VectorXd&, Eigen::VectorXd&), void (T::*errorJacobianFcnPtr)( Eigen::VectorXd&, Eigen::MatrixXd& ), jacobian_update_method update_method = EVALUATE, int updateJacobianEvery = -1, double modelReductionThreshold = -1 );


    /**
     * @brief solve Finds the solution to the problem
     * @param initial_guess The inital states guess
     * @param tolerance The termination tolerance
     * @return  The solution
     */
    Eigen::VectorXd solve(const Eigen::VectorXd& initial_guess, double tolerance = 1.0e-5, unsigned int maxAttempts = 1000, bool displayProgress = false, search_method method = NEWTON_RAPHSON, double initialLambdaValue = 1.0);

    /**
     * @brief solution Returns the current solution.
     * @return The solution
     */
    Eigen::VectorXd solution() const;

    /**
     * @brief solutionNormError Returns the current norm of the error vector (the thing were trying to minimize)
     * @return norm error (in a least squares sense)
     */
    double solutionRMSError() const;

    /**
     * @brief changeUpdateMethod
     * @param newMethod
     */
    void changeUpdateMethod(jacobian_update_method newMethod, int updateJacobianEvery = -1);

    void setModelReductionThreshold(double thresh = -1);

    //     ////////////// This functionality has been removed /////////////////
    //    /**
    //     * @brief changeJacobianFunction
    //     */
    //    void changeJacobianFunction( Eigen::MatrixXd (*errorJacobianFcnPtr)( const Eigen::VectorXd& ) );
    //
    //    /**
    //     * @brief changeErrorFunction
    //     */
    //    void changeErrorFunction( Eigen::VectorXd (*errorFcnPtr)(const Eigen::VectorXd&) );
    // ////////////////////////////////////////////////////////////////////////
};


template< class T >
NewtonRaphson<T>::NewtonRaphson( T* classPointer, void (T::*errorFcnPtr)(Eigen::VectorXd&, Eigen::VectorXd&), void (T::*errorJacobianFcnPtr)( Eigen::VectorXd&, Eigen::MatrixXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ )
{

    this->classObjectPointer = classPointer;
    this->class_errorFcnPtr = errorFcnPtr;
    this->class_errorJacobianFcnPtr = errorJacobianFcnPtr;
    //this->class_stateConstraintFcnPtr = stateConstraintFcnPtr;

    this->modelReductionThreshold = modelReductionThreshold_;

    iteration = 0;

    H.resize(0,0); /**< The scaling Matrix **/
    J.resize(0,0); /**< The current jacobian **/
    J_inv.resize(0,0); /**< The current inverse jacobian **/
    states.resize(0); /**< The current state vector **/
    delta_States.resize(0); /**< The last change in states **/
    error_last.resize(0); /**< The previous step's error vector **/
    error_this.resize(0); /**< The current error vector **/
    delta_error.resize(0); /**< The current change in error **/

    this->errorFcnPtr = 0;
    this->errorJacobianFcnPtr = 0;
    this->errorFcnPtr_2 = 0;
    this->errorJacobianFcnPtr_2 = 0;
    this->update_method = update_method;

    this->modelReductionThreshold = -1;

    this->updateJacobianEvery = updateJacobianEvery;

    sqrtEps = std::sqrt(std::numeric_limits<double>::epsilon());

}

template<>
NewtonRaphson<_nr_default>::NewtonRaphson(jacobian_update_method update_method, int updateJacobianEvery );

template<class T>
NewtonRaphson<T>::NewtonRaphson(jacobian_update_method update_method, int updateJacobianEvery )
{
    cout << "You are using the wrong NewtonRapson constructor.  Use the class specific one!" << endl;
    assert(false);
}

template<>
NewtonRaphson<_nr_default>::NewtonRaphson(  Eigen::VectorXd (*errorFcnPtr)(const Eigen::VectorXd&), Eigen::MatrixXd (*errorJacobianFcnPtr)( const Eigen::VectorXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ );

template<class T>
NewtonRaphson<T>::NewtonRaphson(  Eigen::VectorXd (*errorFcnPtr)(const Eigen::VectorXd&), Eigen::MatrixXd (*errorJacobianFcnPtr)( const Eigen::VectorXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ )
{
    cout << "You are using the wrong NewtonRapson constructor.  Use the class specific one!" << endl;
    assert(false);
}

template<>
NewtonRaphson<_nr_default>::NewtonRaphson(  void (*errorFcnPtr)(Eigen::VectorXd&, Eigen::VectorXd&), void (*errorJacobianFcnPtr)( Eigen::VectorXd&, Eigen::MatrixXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ );

template<class T>
NewtonRaphson<T>::NewtonRaphson(  void (*errorFcnPtr)(Eigen::VectorXd&, Eigen::VectorXd&), void (*errorJacobianFcnPtr)( Eigen::VectorXd&, Eigen::MatrixXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ )
{
    cout << "You are using the wrong NewtonRapson constructor.  Use the class specific one!" << endl;
    assert(false);
}

template<class T>
void NewtonRaphson<T>::calculateDeltaStates( )
{
    // Timing::stopWatch tic;

    int H_rows = H.rows();
    int H_cols = HtoJ_index.size();

    assert( H_rows > 0 && H_cols > 0 );

    Eigen::VectorXd ds_tmp(Eigen::VectorXd::Zero(H_cols));

    if( modelReductionThreshold <= 0 )
    {

        if( H_rows >= H_cols )
        {
            if( updateJacobianDecomposition )
            {
                if( method == NEWTON_RAPHSON || method == NEWTON_RAPHSON_QUADRATIC_CORRECT )
                    J_QR.compute( H.block(0,0,H_rows,H_cols).transpose()*H.block(0,0,H_rows,H_cols) );
                else
                {
                    Eigen::MatrixXd HtH(H.block(0,0,H_rows,H_cols).transpose()*H.block(0,0,H_rows,H_cols));
                    //Eigen::MatrixXd lambdaDiag((Lambda*HtH.diagonal()).asDiagonal());
                    for( int i=0;i<HtH.rows();i++)
                        HtH(i,i) *= 1+Lambda;

                    J_QR.compute( HtH );//Lambda*Eigen::MatrixXd::Identity(H_cols,H_cols) );
                }
                //cout << H.block(0,0,H_rows,H_cols).transpose()*H.block(0,0,H_rows,H_cols) << endl<<endl;
                updateJacobianDecomposition = false;
            }
            ds_tmp = -J_QR.solve(H.block(0,0,H_rows,H_cols).transpose()*error_this);
        }
        else
        {
            if ( updateJacobianDecomposition )
            {
                if( method == NEWTON_RAPHSON || method == NEWTON_RAPHSON_QUADRATIC_CORRECT)
                    J_QR.compute( H.block(0,0,H_rows,H_cols)*H.block(0,0,H_rows,H_cols).transpose() );
                else
                    J_QR.compute( H.block(0,0,H_rows,H_cols)*H.block(0,0,H_rows,H_cols).transpose() + Lambda*Eigen::MatrixXd::Identity(H_rows,H_rows) );

                updateJacobianDecomposition = false;
            }
            ds_tmp = -H.block(0,0,H_rows,H_cols).transpose()*J_QR.solve(error_this);
        }

    } else
    {
        if( method == NEWTON_RAPHSON || method == NEWTON_RAPHSON_QUADRATIC_CORRECT )
            ds_tmp = -Math::pseudoInverse(H.block(0,0,H_rows,H_cols), error_this, modelReductionThreshold, 0);
        else
            ds_tmp = -Math::pseudoInverse(H.block(0,0,H_rows,H_cols), error_this, modelReductionThreshold, Lambda);

    }

    delta_States.setZero(delta_States.rows(),delta_States.cols());
    for( unsigned int i=0; i<(unsigned int)H_cols; i++ )
    {
        delta_States(HtoJ_index[i]) = ds_tmp(i)/H_col_scale[i];
    }

    // rank one updates if requried
    switch( update_method )
    {
    case EVALUATE:
    {
        break;
    }
    case APPROXIMATE_BROYDEN_1:
    {
        delta_States += -1.0*J_inv * error_this;
        break;
    }
    case APPROXIMATE_BROYDEN_2:
    {
        delta_States += -1.0*J_inv * error_this;
        break;
    }
    }

    //std::cout << "Time for Inverse: " << tic.stop() << std::endl;
}


template< class T >
Eigen::VectorXd NewtonRaphson<T>::numericalDeritiveHelper( const Eigen::VectorXd & state)
{
    Eigen::VectorXd retValue(error_this.rows(), error_this.cols());
    Eigen::VectorXd tmp = state;

    if( errorFcnPtr )
        retValue = (*errorFcnPtr)(tmp);
    else if( errorFcnPtr_2 )
        (*errorFcnPtr_2)(tmp, retValue);
    else
        (classObjectPointer->*class_errorFcnPtr)(tmp, retValue);

    //(classObjectPointer->*class_errorFcnPtr)(tmp,retValue);

    return retValue;
}

template< class T >
void NewtonRaphson<T>::updateJacobian( int external_count )
{
    // update inverse jacobian
    jacobian_update_method method = update_method;

    if( external_count == 0 || update_method == EVALUATE || updateJacobianEvery == 0 || external_count%updateJacobianEvery == 0 || J_inv.rows()==0 || J_inv.cols()==0)
    {
        method = EVALUATE; // recalculate the jacobian
    }

    switch( method )
    {
    case EVALUATE:
    {
        //tic.start();
        countJacobEval ++;
        if( errorJacobianFcnPtr )
        {
            J = (*errorJacobianFcnPtr)(states);
        }
        else if( errorJacobianFcnPtr_2 )
        {
            (*errorJacobianFcnPtr_2)(states, J);
        }
        else if( class_errorJacobianFcnPtr && classObjectPointer )
        {
            (classObjectPointer->*class_errorJacobianFcnPtr)(states,J);
        }
        else
        {
            calculateNumericalJacobian();
            //J = CalculateJacobian< NewtonRaphson<T> >( *this, &NewtonRaphson<T>::numericalDeritiveHelper, states );
        }
        //std::cout << "Time to Build Jacobian: " << tic.stop() << std::endl;

//        Eigen::JacobiSVD<Eigen::MatrixXd> svd(J,Eigen::ComputeFullU | Eigen::ComputeFullV);
//        std::cout << "J is: " << J.rows() << " x " << J.cols() << std::endl << std::endl;
//        std::cout << svd.singularValues() << std::endl << std::endl;
//        std::cout << svd.matrixV().rightCols(4) << std::endl << std::endl;
//        std::cout << J.rightCols(4) << std::endl << std::endl;

        // Remove Zero columns
        HtoJ_index.clear();
        H_col_scale.clear();
        H.setZero(J.rows(),J.cols());
        unsigned int count = 0;
        for( unsigned int i=0; i<J.cols(); i++ )
        {
            double J_col_norm = J.col(i).norm();
            //cout << "Col: " << i << " Norm: " << J_col_norm << (J_col_norm > std::numeric_limits<double>::epsilon()*10.0?"\tKeep":"\tLoose")<< endl;
            if( J_col_norm > std::numeric_limits<double>::epsilon()*10.0 /*!=0*/ )
            {
                H.col(count) = J.col(i)/J_col_norm;
                H_col_scale.push_back(J_col_norm);
                HtoJ_index.push_back(i);
                count ++;
            }
        }

        if( update_method == APPROXIMATE_BROYDEN_1 )
            J_inv = sqrtEps*J_inv.setIdentity(J.cols(),J.rows());
        else
            J_inv.setZero(J.cols(),J.rows());

        updateJacobianDecomposition = true;

        //cout << "Size J: " << J.rows() << "x" << J.cols() << endl;
        //cout << "Size H: " << H.rows() << "x" << count << endl;

        //Eigen::MatrixXd HtH = H.block(0,0,H.rows(),count).transpose()*H.block(0,0,H.rows(),count);
        //cout << "Size H'H: " << HtH.rows() << "x" << HtH.cols() << endl;
        /*
        if( modelReductionThreshold <= 0 && method == NEWTON_RAPHSON )
        {
            if( H.rows() >= count )
                J_QR.compute(H.block(0,0,H.rows(),count).transpose()*H.block(0,0,H.rows(),count));
            else
                J_QR.compute(H.block(0,0,H.rows(),count)*H.block(0,0,H.rows(),count).transpose());
        }
        */

        break;
    }
    case APPROXIMATE_BROYDEN_1:
        //https://en.wikipedia.org/wiki/Broyden%27s_method
        if( external_count != 0 )
        {
            Eigen::VectorXd j_inv_delta_f = J_inv*(delta_error);
            J_inv += (delta_States - j_inv_delta_f)/(delta_States.dot(j_inv_delta_f))*(delta_States.transpose()*J_inv);
        }
        break;

    case APPROXIMATE_BROYDEN_2:
        //https://en.wikipedia.org/wiki/Broyden%27s_method
        Eigen::VectorXd j_inv_delta_f = J_inv*delta_error;
        J_inv += (delta_States - j_inv_delta_f)/(pow(delta_error.norm(),2.0))*(delta_error.transpose());
        break;
    }
}

template< class T >
void NewtonRaphson<T>::updateError()
{
    error_last = error_this;

    returnError( states, error_this );
}

template< class T >
void NewtonRaphson<T>::returnError(Eigen::VectorXd& state, Eigen::VectorXd& error_out )
{
    countErrorEval ++;

    if( errorFcnPtr )
        error_out = (*errorFcnPtr)(state);
    else if( errorFcnPtr_2 )
        (*errorFcnPtr_2)(state, error_out);
    else
        (classObjectPointer->*class_errorFcnPtr)(state, error_out);

//    cout << error_out << endl;
}

template< class T >
void NewtonRaphson<T>::calculateNumericalJacobian()
{
    Eigen::VectorXd errShifted(Eigen::VectorXd::Zero(error_this.rows(), 1));
    Eigen::VectorXd stateTmp(states);
    assert( error_this.cols() == 1 && states.cols() == 1);

    if(J.rows() != error_this.rows() || J.cols() != states.rows() )
        J.setZero(error_this.rows(), states.rows() );

    for( unsigned int i=0; i< states.rows(); i++ )
    {
        double delta = std::max(sqrtEps, std::abs(stateTmp(i))*sqrtEps);
        stateTmp(i) += delta;
        returnError(stateTmp, errShifted);
        J.col(i) = (errShifted-error_this)/delta;
        stateTmp(i) = states(i);
    }
}

template< class T >
Eigen::VectorXd NewtonRaphson<T>::solve(const Eigen::VectorXd& initial_guess, double tolerance, unsigned int maxAttempts, bool displayProgress, search_method useMethod, double initialLambdaValue)
{
    assert( ("Initial Lambda value must be greator than zero", useMethod ==NEWTON_RAPHSON || useMethod ==NEWTON_RAPHSON_QUADRATIC_CORRECT || initialLambdaValue > 0 ));
    method = useMethod;
    countErrorEval = 0;
    countJacobEval = 0;
    /*
    if( initial_guess.rows() == states.rows() )
        iteration &= iteration && (initial_guess-states).norm()/states.norm() < .01;
    else
        iteration = 0;

    if(iteration == 0 )
    {
        states = initial_guess;
        updateError();

        states_last = states;
        delta_States.setZero(initial_guess.rows(),initial_guess.cols());
        delta_error.setZero(error_this.rows());
        error_last = error_this;

        rmsError_this = solutionRMSError();
        rmsError_last = rmsError_this;

        delta_States_last.setZero(delta_States.rows(),1);

    } else
    {
        states = initial_guess;
        updateError();
        delta_States = states-states_last;
        delta_error = error_this-error_last;

        rmsError_last = rmsError_this;
        rmsError_this = solutionRMSError();

        delta_States_last.setZero(delta_States.rows(),1);

    }
*/
    iteration = 0;
    states = initial_guess;
    states_last = states;
    delta_States.setZero(initial_guess.rows(),initial_guess.cols());
    delta_States_last.setZero(initial_guess.rows(),initial_guess.cols());

    updateError();
    delta_error.setZero(error_this.rows(),error_this.cols());
    error_last = error_this;
    rmsError_this = solutionRMSError();
    rmsError_last = rmsError_this;


    // parameter initialization for adaptive LM
    switch( method )
    {
    case NEWTON_RAPHSON:
        Lambda = 0;
        break;
    case NEWTON_RAPHSON_PLUS_BRENTS:
        Lambda = 0;
        break;
    case NEWTON_RAPHSON_QUADRATIC_CORRECT:
        Lambda = 0;
        break;
    case LEVENBERG_MARQUARDT:
        Lambda = initialLambdaValue;
        break;
    case LEVENBERG_MARQUARDT_ANTICIPATIVE:
        Lambda = sqrtEps;
        break;
    case LEVENBERG_MARQUARDT_ADAPTIVE:
        Lambda = initialLambdaValue;//1.0;
        break;
    default:
        assert(("Unlknown solving method specified",false));
        break;
    }

    LambdaChangeRange = 1 + sqrtEps;
    lambdaRangeCount = 0;

    double deltaRmsError = rmsError_this = rmsError_last;

    bool oscilitory = false;
    bool jacobianUpdateRequired = true;
    bool converged;
    double predictionFactor = 1.0;
    double percentStateChange = 0.0;
    double percentErrorChange = 0.0;
    int errorIncreaseCounter = 0;

    if( displayProgress )
    {
        std::cout << std::setprecision(3) << "Iteration: " << iteration << "   RMS Error: " << rmsError_this << std::endl;// << "   Delta RMS Error: " <<  deltaRmsError  << "\tError Change: " << percentErrorChange << "%" << "\tState Change: " << percentStateChange << "%";
    }

    do {
        converged = true;

        if( jacobianUpdateRequired )
        {
            updateJacobian(iteration);
        }

        /*if( !initialPass && (delta_States+delta_States_last).norm() < tolerance * 10.0  )
        {
            //cout << "d states: " << (delta_States+delta_States_last).norm() << "\t\t" << (delta_States+delta_States_last).transpose() << endl;

            delta_States = (delta_States+delta_States_last)/2.0;
            oscilitory = true;
        }*/

        delta_States_last = delta_States;

        calculateDeltaStates();

        states += delta_States;
        jacobianUpdateRequired = true;


        removeNullSpace();

        Eigen::VectorXd delta_error_last = delta_error;
        Eigen::VectorXd error_last_last = error_last;
        updateError();

        if(method == NEWTON_RAPHSON_QUADRATIC_CORRECT )
        {
            // errorAct = a*f^2+b*f+c
            // c = rms_error_last
            // b = delta_error_prediction
            // a = deltaRMSError-delta_error_prediction
            // min => -b/a;

            //double tmp1 = solutionRMSError();

            double c = error_last.squaredNorm();

            delta_error = J*delta_States;
            double b = 2.0*error_last.dot(delta_error);

            double e_this_sqNorm = error_this.squaredNorm();
            double a = e_this_sqNorm-c-b;

            double factor = -b/(2.0*a);

            if( a == 0 || (0.9 <= factor && factor <= 1.1) )
                factor = 1;
            else if(factor < 0 || factor > 1.1 )
            {
                double deltaError = e_this_sqNorm - c;

                std::vector<double> e,ac,bc,cc;
                e.push_back(c);
                e.push_back(e_this_sqNorm);
                ac.push_back(0);
                ac.push_back(1);
                bc.push_back(0);
                bc.push_back(1);
                cc.push_back(1);
                cc.push_back(1);

                if( factor < 0 )
                {
                    factor = 1;

                    while (  deltaError < 0 )
                    {
                        factor *= 3;
                        states = states_last + factor*delta_States;
                        updateError();
                        double e_this_sqNorm2 = error_this.squaredNorm();

                        e.push_back(e_this_sqNorm2);
                        ac.push_back(factor*factor);
                        bc.push_back(factor);
                        cc.push_back(1.0);

                        deltaError = e_this_sqNorm2-e_this_sqNorm;
                        e_this_sqNorm = e_this_sqNorm2;
                    }
                }
                else // factor is greator than 1
                {
                    states = states_last + factor*delta_States;
                    updateError();
                    double e_this_sqNorm2 = error_this.squaredNorm();

                    e.push_back(e_this_sqNorm2);
                    ac.push_back(factor*factor);
                    bc.push_back(factor);
                    cc.push_back(1.0);
                }

                Eigen::MatrixXd A(e.size(),3);
                Eigen::VectorXd E(e.size(),1);
                for( int i=0; i< e.size(); i++ )
                {
                    A.block(i,0,1,3) << ac[i], bc[i], cc[i];
                    E(i) = e[i];
                }
                Eigen::Vector3d coeff = A.householderQr().solve(E);

                factor = -coeff(1)/(2.0*coeff(0));

            }

            if( factor < sqrtEps )
                factor = sqrtEps;

            if( factor != 1 )
            {
                delta_States *= factor;
                states = states_last + delta_States;
                updateError();
                error_last = error_last_last;
            }
            //double tmp2 = solutionRMSError();

            //cout << "\tFactor: " << factor << endl << endl;

        }

        double rmsError_last_last = rmsError_last;
        rmsError_last = rmsError_this;
        rmsError_this = solutionRMSError();

        delta_error = error_this-error_last;

        deltaRmsError = rmsError_this-rmsError_last;

        delta_error_prediction = std::sqrt((error_last+J*delta_States).squaredNorm()/error_this.rows())-rmsError_last;
        predictionFactor = deltaRmsError/( delta_error_prediction );

        delta_States = states-states_last; // just incase BC/Constraint inforcement changed something

        percentErrorChange = deltaRmsError / rmsError_this * 100.0;
        if( std::isnan(percentErrorChange) )
            percentErrorChange = deltaRmsError;

        double normStates = states.norm();
        percentStateChange = delta_States.norm()/normStates * 100.0;
        if( std::isnan(percentStateChange) )
            percentStateChange = sqrt(delta_States.squaredNorm()/delta_States.rows());


        if( deltaRmsError <= 0  || method == NEWTON_RAPHSON || method == NEWTON_RAPHSON_QUADRATIC_CORRECT || method == NEWTON_RAPHSON_PLUS_BRENTS)
        {
            errorIncreaseCounter = 1;

            if( method == LEVENBERG_MARQUARDT )
            {
                Lambda /= 2.0;
            }
            else if( method == LEVENBERG_MARQUARDT_ANTICIPATIVE )
            {
                double tmp = std::sqrt(1.0+4.0*std::abs(delta_error_prediction-deltaRmsError)/delta_States.squaredNorm()*rmsError_this)-1.0;
                Lambda = std::max( tmp, sqrtEps );
            }
            else if( method == LEVENBERG_MARQUARDT_ADAPTIVE )
            {

                // UPDATE Lambda based on convergance criteria

                // handel boundry conditions that occure at initialization
                if( std::isinf(predictionFactor) || std::isnan(predictionFactor)  )
                {
                    predictionFactor = 1.0-3.0/4.0*LambdaChangeRange;
                }

                if( predictionFactor > 1+LambdaChangeRange/4.0 )
                {
                    double e_this_sqNorm = error_this.squaredNorm();
                    double deltaError = e_this_sqNorm-error_last.squaredNorm();

                    std::vector<double> e,ac,bc,cc;
                    e.push_back(e_this_sqNorm-deltaError);
                    e.push_back(e_this_sqNorm);
                    ac.push_back(0);
                    ac.push_back(1);
                    bc.push_back(0);
                    bc.push_back(1);
                    cc.push_back(1);
                    cc.push_back(1);

                    double factor = 1.0;
                    while (  deltaError < 0 )
                    {
                        factor *= 3.0;
                        states = states_last + factor*delta_States;
                        updateError();
                        double e_this_sqNorm2 = error_this.squaredNorm();

                        e.push_back(e_this_sqNorm2);
                        ac.push_back(factor*factor);
                        bc.push_back(factor);
                        cc.push_back(1.0);

                        deltaError = e_this_sqNorm2-e_this_sqNorm;
                        e_this_sqNorm = e_this_sqNorm2;
                    }

                    Eigen::MatrixXd A(e.size(),3);
                    Eigen::VectorXd E(e.size(),1);
                    for( int i=0; i< e.size(); i++ )
                    {
                        A.block(i,0,1,3) << ac[i], bc[i], cc[i];
                        E(i) = e[i];
                    }
                    Eigen::Vector3d coeff = A.householderQr().solve(E);

                    factor = -coeff(1)/(2.0*coeff(0));

                    delta_States *= factor;
                    states = states_last + delta_States;
                    updateError();
                    delta_States = states-states_last; // just incase BC/Constraint inforcement changed something
                    jacobianUpdateRequired = true;

                    error_last = error_last_last;
                    rmsError_this = solutionRMSError();

                    delta_error = error_this-error_last;

                    deltaRmsError = rmsError_this-rmsError_last;

                    lambdaRangeCount ++;

                    percentErrorChange = deltaRmsError / rmsError_this * 100.0;
                    if( std::isnan(percentErrorChange) )
                        percentErrorChange = deltaRmsError;

                    double normStates = states.norm();
                    percentStateChange = delta_States.norm()/normStates * 100.0;
                    if( std::isnan(percentStateChange) )
                        percentStateChange = sqrt(delta_States.squaredNorm()/delta_States.rows());

                }
                else if( std::abs(1.0-predictionFactor) > LambdaChangeRange )
                {
                    // Error change was too small compared to a linear change. Getting close to an over-step.  Be more gradient decent.
                    double factor =  2;//std::abs(1.0-predictionFactor)/LambdaChangeRange;
                    Lambda = std::max( Lambda*factor, std::numeric_limits<double>::epsilon() );

                }
                else if( std::abs(1.0-predictionFactor) < LambdaChangeRange/2.0 )
                {
                    // Error change was too close to linear.  Get more aggresive!
                    double factor  =  1/2.0;//std::sqrt( std::abs(1.0-predictionFactor)/(LambdaChangeRange/2.0) );
                    Lambda = std::max( Lambda*factor, std::numeric_limits<double>::epsilon() );

                    lambdaRangeCount ++;

                }

                if( lambdaRangeCount > 3 )
                {
                    LambdaChangeRange = std::min(1.0, LambdaChangeRange*2);
                    lambdaRangeCount = 0;
                }

            }

            states_last = states;
        }
        else
        {
            errorIncreaseCounter ++;
            // change made things worse!

            if( method == LEVENBERG_MARQUARDT )
            {
                Lambda *= 2.0;
            }
            else if( method == LEVENBERG_MARQUARDT_ANTICIPATIVE )
            {
                //double secondDirEst = std::abs(delta_error_prediction-deltaRmsError)/delta_States.squaredNorm();

                double tmp = std::sqrt(1.0+4.0*std::abs(delta_error_prediction-deltaRmsError)/delta_States.squaredNorm()*rmsError_this)-1.0;
                Lambda = std::max( std::max(tmp, Lambda*2), sqrtEps );
            }
            else if( method == LEVENBERG_MARQUARDT_ADAPTIVE )
            {
                Lambda *= 2.0*errorIncreaseCounter;
                //Lambda *= std::abs(deltaRmsError-delta_error_prediction)/Lambda+2.0;//2.0;
                LambdaChangeRange = std::max(LambdaChangeRange*0.75, sqrtEps);
            }

            states = states_last;
            delta_States = delta_States_last;

            error_this = error_last;
            error_last = error_last_last;
            delta_error = delta_error_last;

            rmsError_this = rmsError_last;
            rmsError_last = rmsError_last_last;

            deltaRmsError = rmsError_this-rmsError_last;

            delta_States_last.setZero(delta_States_last.rows(),delta_States_last.cols());
            lambdaRangeCount = 0;

            jacobianUpdateRequired = false;
            updateJacobianDecomposition = true;

            converged = false;
        }


        /*if( !initialPass && (error_this+error_last).norm() < tolerance*10.0 )
        {
            // cout << "e Error: " << (error_this+error_last).norm() << "\t\t" << (error_this+error_last).transpose() << endl;

            //delta_States = (delta_States+delta_States_last)/2.0;
            oscilitory = true;
        }*/




       // cout << std::numeric_limits<double>::epsilon()*normStates << endl;

        iteration ++;

        assert( !std::isinf(Lambda) && !std::isnan(Lambda));
        assert( !std::isinf(rmsError_this) && !std::isnan(rmsError_this));


        converged = converged && (
                    oscilitory ||
                    iteration >= maxAttempts ||
                    rmsError_this < tolerance ||
                    (std::abs(deltaRmsError)*(1+Lambda) < tolerance && delta_error.norm()!=0) ||
                    (std::abs(percentErrorChange)*(1+Lambda)/100.0 < tolerance && percentErrorChange >= std::numeric_limits<double>::epsilon()*rmsError_this ) ||
                    (std::abs(percentStateChange)*(1+Lambda)/100.0 < tolerance && percentStateChange >= std::numeric_limits<double>::epsilon()*normStates ) );

        if( iteration == 1 && maxAttempts > 1 )
            converged = false;

        if( Lambda > 10000 )
            converged = true;


        if( displayProgress )
        {
            std::cout <</* "Conv: " << (converged?"YES   ":"NO   ") << */std::setprecision(3) << "Iteration: " << iteration << "   RMS Error: " << rmsError_this << "   Delta RMS Error: " <<  deltaRmsError  << "   Error Change: " << percentErrorChange << "%" << "   State Change: " << percentStateChange << "%";
            if( method != NEWTON_RAPHSON && method != NEWTON_RAPHSON_QUADRATIC_CORRECT)
                std::cout << std::setprecision(3) << "\t" << "Lambda: " << Lambda << "\tP. Factor: " << predictionFactor << "\tP. Factor Range: " <<  LambdaChangeRange;

            std::cout << std::endl;
        }


    }while( !converged );

    if( displayProgress )
    {
        std::cout << std::setprecision(3) << "Iteration: " << iteration << "\tRMS Error: " << rmsError_this << "\tDelta RMS Error: " <<  deltaRmsError  << "\tPercent Error Change: " << percentErrorChange << "%" << "\tPercent State Change: " << percentStateChange << "%\tTolerance: " << tolerance <<  "\tError Evaluations: " << countErrorEval << "  Jacobian Evaluations: " << countJacobEval << "\tDONE!"<<std::endl<<std::endl;
    }


    return states;
}

template< class T >
Eigen::VectorXd NewtonRaphson<T>::solution() const
{
    return states;
}
template< class T >
double NewtonRaphson<T>::solutionRMSError() const
{
    return std::sqrt(error_this.squaredNorm()/error_this.rows());
}

template< class T >
void NewtonRaphson<T>::changeUpdateMethod(jacobian_update_method newMethod, int updateJacobianEvery )
{
    this->updateJacobianEvery = updateJacobianEvery;
    this->update_method = newMethod;
}

template< class T >
void NewtonRaphson<T>::removeNullSpace()
{    
    if( H.rows() < (unsigned int)HtoJ_index.size() && H.rows() == J.rows() && update_method == EVALUATE)
    {
        Eigen::VectorXd ds_tmp(Eigen::VectorXd::Zero(HtoJ_index.size()));
        /*for( unsigned int i=0; i<HtoJ_index.size(); i++ )
        {
            ds_tmp(i) = states(HtoJ_index[i])*H_col_scale[i];;
        }

        if( modelReductionThreshold > 0 || update_method == LEVENBERG_MARQUARDT_ADAPTIVE )
        {
            J_QR.compute(H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose());
            updateJacobianDecomposition = true;
            //ds_tmp = H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*Math::pseudoInverse(H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose(), H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp, modelReductionThreshold);
        }
        //        cout << H.block(0,0,H.rows(),HtoJ_index.size()) << endl << "____" <<endl;
        //        cout << H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose() << endl << "____" <<endl;
        //        cout << (H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose()).inverse() << endl << "____" <<endl;
        //        cout << H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp << endl<< "----" << endl;

        //        cout << (H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose()).inverse()*H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp << endl << "____" <<endl;

        //        cout <<  H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*(H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose()).inverse()*H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp << endl << "____" <<endl;

        //        cout <<  J.transpose()*(J*J.transpose()).inverse()*J*ds_tmp2 << endl << "____" <<endl;

        ds_tmp = H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*J_QR.solve(H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp);

        for( unsigned int i=0; i<HtoJ_index.size(); i++ )
        {
            states(HtoJ_index[i]) = ds_tmp(i)/ H_col_scale[i];
        }*/

        for( unsigned int i=0; i<HtoJ_index.size(); i++ )
        {
            ds_tmp(i) = states(HtoJ_index[i])*H_col_scale[i];
        }

        //        cout << "\n---------------\n"
        //             << "Inital States:\n"<< states.transpose() << "\tLength: " << states.norm()<< endl << "____" <<endl
        //             << "Correct Nullspace Removed States:\n" <<  (J.transpose()*(J*J.transpose()).inverse()*J*states).transpose() << endl << "____" <<endl;

        ds_tmp = J.transpose()*(J*J.transpose()).ldlt().solve(H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp);

        for( unsigned int i=0; i<HtoJ_index.size(); i++ )
        {
            states(HtoJ_index[i]) = ds_tmp(i);
        }
        //        cout << "Calculated Nullspace Removed States:\n" << states.transpose() << "\tLength: " << states.norm() <<endl <<endl;
    }

}

template< class T >
void NewtonRaphson<T>::setModelReductionThreshold(double thresh )
{
    modelReductionThreshold = thresh;
}


//template< class T>
//void NewtonRaphson_C<T>::inforceConstraints()
//{
//    if( class_stateConstraintFcnPtr )
//    {
//        Eigen::VectorXd statesTmp = states;

//        (classObjectPointer->*class_stateConstraintFcnPtr)(states);

//        delta_States += states - statesTmp;
//    }
//}


#endif // NEWTONRAPHSON_H
