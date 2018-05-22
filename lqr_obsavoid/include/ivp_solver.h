#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H

#include <iostream>
#include <sstream>

#include <assert.h>
#include <Eigen/Dense>
#include <vector>

#include <cmath> // for isnan functionality

#include "integrationmethod.h"
#include "NewtonRaphson.h"

#include "utilities.h"

//#include "Utilities/Timer.h"



struct IVP_Solver_Default {};

template< class T = IVP_Solver_Default >
class IVP_Solver
{

public:
    IVP_Solver(const IntegrationMethod& integrationMethod, Eigen::MatrixXd (*ode_function_ptr)( double time, const Eigen::MatrixXd& currentState ), Eigen::MatrixXd (*ode_Jacobian_function_prt)( double time, const Eigen::MatrixXd& currentState ) = 0, void (*state_constraint_function_prt)( double time, Eigen::MatrixXd& currentState ) = 0, void (*state_constraint_jacobian_prt)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian )=0 );
    IVP_Solver(const IntegrationMethod& integrationMethod, T& owningObject, Eigen::MatrixXd (T::*ode_function_ptr)( double time, const Eigen::MatrixXd& currentState ), Eigen::MatrixXd (T::*ode_Jacobian_function_prt)( double time, const Eigen::MatrixXd& currentState ) = 0, void (T::*state_constraint_function_prt)( double time, Eigen::MatrixXd& currentState ) = 0, void (T::*state_constraint_jacobian_prt)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian ) = 0 );
    IVP_Solver( const IVP_Solver& other);
    ~IVP_Solver();

    bool solve(const Eigen::MatrixXd& initialCondition, double to, double tf, double tolerance = 1e-5, unsigned int maxNumSteps = 100 );
    bool isSolved() const;

    std::string getSolverName() const;

    Eigen::MatrixXd getSolution() const; /**< Returns the solution at the final time **/
    Eigen::MatrixXd getSolution(double time) const; /**< Returns the solution at an intermediate time based on the interpolation polynomials **/

    Eigen::MatrixXd getSolutionIntegral() const; /**< Returns the integral of the solution from the initial time to the final time based on the interpolation polynomials **/
    Eigen::MatrixXd getSolutionIntegral(double to, double tf); /**< Returns the integration of the solution from t0 to tf based on the interpolation polynomials **/

    int getNumberOfSteps() const;

    IVP_Solver& operator=(const IVP_Solver& other);


protected:
    const IntegrationMethod* methodPtr;

    struct NodeCoeff
    {
        double t_o;
        double t_f;
        std::vector< typename Eigen::MatrixXd > interpolationCoeff;
    };

    Eigen::MatrixXd solution_initial_state;
    Eigen::MatrixXd solution_final_state;
    std::vector< typename IVP_Solver<T>::NodeCoeff* > solution_node_list;

    // methods for implicit node solving
    NewtonRaphson< IVP_Solver<T> > NR_solver;
    void NR_errorFunction( Eigen::VectorXd& collocatoinStateCurrent, Eigen::VectorXd& errorOut ); /**< Packs the implicit error vector **/
    void NR_errorJacobian( Eigen::VectorXd& collocatoinStateCurrent, Eigen::MatrixXd& jacobOut ); /**< Packs the implicit error jacobian **/
    double packNode( double start, double finish, const Eigen::VectorXd& collocationState, IVP_Solver<T>::NodeCoeff* nodePtr ) const; /**< returns estimated error **/

    // methods for explicit node solving
    double integrateNode( double start, double finish, const Eigen::MatrixXd& IC, Eigen::MatrixXd& state_final, IVP_Solver<T>::NodeCoeff* outputNode ); /**< returns estimated error **/

    T* owningObjectPtr;
    Eigen::MatrixXd (T::*ode_function_clptr)( double time, const Eigen::MatrixXd& currentState );
    Eigen::MatrixXd (T::*ode_jacobian_function_clptr)( double time, const Eigen::MatrixXd& currentState );
    void (T::*state_constraint_function_clptr)( double time, Eigen::MatrixXd& currentState ); /**< Applies the constraints to the state **/
    void (T::*state_constraint_jacobian_clptr)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian ); /** Returns the constraint error and constraint-error-jacobian but does not change state **/

    Eigen::MatrixXd (*ode_function_ptr)( double time, const Eigen::MatrixXd& currentState );
    Eigen::MatrixXd (*ode_jacobian_function_ptr)( double time, const Eigen::MatrixXd& currentState );
    void (*state_constraint_function_ptr)( double time, Eigen::MatrixXd& currentState ); /**< Applies the constraints to the state **/
    void (*state_constraint_jacobian_ptr)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian); /** Returns the constraint error and constraint-error-jacobian but does not change state **/


    Eigen::MatrixXd ode_function( double time, const Eigen::MatrixXd& currentState );
    Eigen::MatrixXd ode_jacobian_function( double time, const Eigen::MatrixXd& currentState );
    void state_constraint_function( double time, Eigen::MatrixXd& currentState ); /**< Applies the constraints to the state **/
    void state_constraint_jacobian( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian ); /** Returns the constraint error and constraint-error-jacobian but does not change state **/

    unsigned int numActiveNodes;
    void clearNodeList();

    // helper function
    static bool inRange(double bound1, double bound2, double point) ;

    bool solutionFound;

};

template< class T >
IVP_Solver<T>::IVP_Solver(const IntegrationMethod& integrationMethod, Eigen::MatrixXd (*ode_function_ptr)( double time, const Eigen::MatrixXd& currentState ), Eigen::MatrixXd (*ode_Jacobian_function_prt)( double time, const Eigen::MatrixXd& currentState ), void (*state_constraint_function_prt)( double time, Eigen::MatrixXd& currentState ), void (*state_constraint_jacobian_prt)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian ) ):
    methodPtr(&integrationMethod),
    NR_solver(this, &IVP_Solver<T>::NR_errorFunction, &IVP_Solver<T>::NR_errorJacobian, EVALUATE, -1, -1)
{
    solution_node_list.clear();
    assert( false ); // Using the wrong constructor for IVP_Solver. You need to pass the object that the functions belong to!
}

template< class T >
IVP_Solver<T>::IVP_Solver(const IntegrationMethod& integrationMethod, T& owningObject, Eigen::MatrixXd (T::*ode_function_ptr)( double time, const Eigen::MatrixXd& currentState ), Eigen::MatrixXd (T::*ode_Jacobian_function_prt)( double time, const Eigen::MatrixXd& currentState ), void (T::*state_constraint_function_prt)( double time, Eigen::MatrixXd& currentState ), void (T::*state_constraint_jacobian_prt)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian ) ):
    methodPtr(&integrationMethod),
    solution_node_list(500, (IVP_Solver<T>::NodeCoeff*) 0 ),
    NR_solver(this, &IVP_Solver<T>::NR_errorFunction, &IVP_Solver<T>::NR_errorJacobian, EVALUATE, -1, -1)
{
    //solution_node_list = std::vector< typename IVP_Solver<T>::NodeCoeff* >(500, (IVP_Solver<T>::NodeCoeff*)0);

    numActiveNodes = 0;

    solutionFound = false;

    this->owningObjectPtr = &owningObject;
    this->ode_function_clptr = ode_function_ptr;
    this->ode_jacobian_function_clptr = ode_Jacobian_function_prt;
    this->state_constraint_function_clptr = state_constraint_function_prt;
    this->state_constraint_jacobian_clptr = state_constraint_jacobian_prt;

    this->ode_function_ptr = 0;
    this->ode_jacobian_function_ptr = 0;
    this->state_constraint_function_ptr = 0;
    this->state_constraint_jacobian_ptr = 0;


    for( int i=0; i<500; i++ )
        solution_node_list[i] = new IVP_Solver<T>::NodeCoeff();


    assert( this->ode_function_clptr !=0 );


}

template<class T>
IVP_Solver<T>::IVP_Solver( const IVP_Solver<T>& other):
    methodPtr(other.methodPtr),
    solution_initial_state(other.solution_initial_state),
    solution_final_state(other.solution_final_state),
    NR_solver(other.NR_solver),
    owningObjectPtr(other.owningObjectPtr),
    ode_function_clptr(other.ode_function_clptr),
    ode_jacobian_function_clptr(other.ode_jacobian_function_clptr),
    state_constraint_function_clptr(other.state_constraint_function_clptr),
    state_constraint_jacobian_clptr(other.state_constraint_jacobian_clptr),
    ode_function_ptr(other.ode_function_ptr),
    ode_jacobian_function_ptr(other.ode_jacobian_function_ptr),
    state_constraint_function_ptr(other.state_constraint_function_ptr),
    state_constraint_jacobian_ptr(other.state_constraint_jacobian_ptr),
    numActiveNodes(other.numActiveNodes),
    solutionFound(other.solutionFound)
{
    for( int i=0; i<other.solution_node_list.size() && i < other.numActiveNodes; i++ )
    {
        IVP_Solver<T>::NodeCoeff*  newCoeff = 0;
        IVP_Solver<T>::NodeCoeff*  oldCoeff = other.solution_node_list[i];
        if( oldCoeff )
        {
            newCoeff = new IVP_Solver<T>::NodeCoeff();
            newCoeff->interpolationCoeff = oldCoeff->interpolationCoeff;
            newCoeff->t_f = oldCoeff->t_f;
            newCoeff->t_o = oldCoeff->t_o;
        }
        solution_node_list.push_back(newCoeff);
    }

}

template< class T >
IVP_Solver<T>::~IVP_Solver()
{
    clearNodeList();
}

template< class T >
IVP_Solver<T>& IVP_Solver<T>::operator=(const IVP_Solver<T>& other)
{
    methodPtr = other.methodPtr;
    solution_initial_state = (other.solution_initial_state);
    solution_final_state = other.solution_final_state;
    NR_solver=(other.NR_solver);
    owningObjectPtr=(other.owningObjectPtr);
    ode_function_clptr=(other.ode_function_clptr);
    ode_jacobian_function_clptr=(other.ode_jacobian_function_clptr);
    state_constraint_function_clptr=(other.state_constraint_function_clptr);
    state_constraint_jacobian_clptr=(other.state_constraint_jacobian_clptr);
    ode_function_ptr=(other.ode_function_ptr);
    ode_jacobian_function_ptr=(other.ode_jacobian_function_ptr);
    state_constraint_function_ptr=(other.state_constraint_function_ptr);
    state_constraint_jacobian_ptr=(other.state_constraint_jacobian_ptr);
    numActiveNodes=(other.numActiveNodes);
    solutionFound=(other.solutionFound);

    for( int i=0; i<other.solution_node_list.size() && i < other.numActiveNodes; i++ )
    {
        IVP_Solver<T>::NodeCoeff*  newCoeff = 0;
        IVP_Solver<T>::NodeCoeff*  oldCoeff = other.solution_node_list[i];
        if( oldCoeff )
        {
            newCoeff = new IVP_Solver<T>::NodeCoeff();
            newCoeff->interpolationCoeff = oldCoeff->interpolationCoeff;
            newCoeff->t_f = oldCoeff->t_f;
            newCoeff->t_o = oldCoeff->t_o;
        }
        solution_node_list.push_back(newCoeff);
    }


    return *this;
}


template< class T >
void IVP_Solver<T>::clearNodeList()
{
    solutionFound = false;
    while( solution_node_list.size() > 0 )
    {
        if( solution_node_list.back() )
            delete solution_node_list.back();

        solution_node_list.back() = 0;
        solution_node_list.pop_back();
    }
    numActiveNodes = 0;
}

template< class T >
Eigen::MatrixXd       IVP_Solver<T>::ode_function( double time, const Eigen::MatrixXd& currentState )
{

    if( owningObjectPtr )
    {
        return (owningObjectPtr->*ode_function_clptr)(time, currentState);
    }

    return (*ode_function_ptr)(time, currentState);
}

template< class T >
Eigen::MatrixXd       IVP_Solver<T>::ode_jacobian_function( double time, const Eigen::MatrixXd& currentState )
{
    if( ode_jacobian_function_clptr !=0 || ode_jacobian_function_ptr!=0 )
    {
        if( owningObjectPtr )
            return (owningObjectPtr->*ode_jacobian_function_clptr)(time, currentState);
        else
            return (*ode_jacobian_function_ptr)(time, currentState);
    }

    std::cout << "ERROR: IVP_Solver<T>::ode_jacobian_function numerical jacobian estimation not implemented yet" << std::endl;
    assert(false); // not implemented yet

    return Eigen::MatrixXd::Zero(0);
}

template< class T >
void      IVP_Solver<T>::state_constraint_function( double time, Eigen::MatrixXd& currentState )
{
    if( 0 != state_constraint_function_ptr || 0 !=  state_constraint_function_clptr )
    {
        if( owningObjectPtr )
            (owningObjectPtr->*state_constraint_function_clptr)(time, currentState);

        else
            (*state_constraint_function_ptr)(time, currentState);
    }
}

template< class T >
void       IVP_Solver<T>::state_constraint_jacobian( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian )
{
    if( 0 != state_constraint_jacobian_ptr || 0 !=  state_constraint_jacobian_clptr )
    {
        if( owningObjectPtr )
            (owningObjectPtr->*state_constraint_jacobian_clptr)(time, currentState, constraintError, constraintJacobian);
        else
            (*state_constraint_jacobian_ptr)(time, currentState, constraintError, constraintJacobian);
    }

}

template< class T >
void IVP_Solver<T>::NR_errorFunction( Eigen::VectorXd& /*collocatoinStateCurrent*/, Eigen::VectorXd& /*errorOut*/ )
{
    /**< Packs the implicit error vector **/
    std::cout << "ERROR: IVP_Solver<T>::NR_errorFunction not implemented yet" << std::endl;
    assert(false); // not implemented yet

}

template< class T >
void IVP_Solver<T>::NR_errorJacobian( Eigen::VectorXd& /*collocatoinStateCurrent*/, Eigen::MatrixXd& /*jacobOut*/ )
{
    /**< Packs the implicit error jacobian **/
    std::cout << "ERROR: IVP_Solver<T>::NR_errorJacobian not implemented yet" << std::endl;

    assert(false); // not implemented yet

}

template< class T >
double IVP_Solver<T>::packNode( double start, double finish, const Eigen::VectorXd& collocationState, IVP_Solver<T>::NodeCoeff* nodePtr ) const
{
    /**< returns estimated error **/
    std::cout << "ERROR: IVP_Solver<T>::packNode not implemented yet" << std::endl;
    assert(false); // not implemented yet

    return nan("");
}

template< class T >
double IVP_Solver<T>::integrateNode( double start, double finish, const Eigen::MatrixXd& IC, Eigen::MatrixXd& state_final_ ,IVP_Solver<T>::NodeCoeff* outputNode )
{
    /**< returns estimated error **/


    int num_row = IC.rows();
    int num_col = IC.cols();
    Eigen::MatrixXd state_final = Eigen::MatrixXd::Zero(num_row,num_col);

    // initialize internode slopes
    std::vector< Eigen::MatrixXd > K_list( methodPtr->numStages(), state_final );

    double h = finish - start;


    //cout << state_final <<endl <<endl;
    //cout << IC <<endl <<endl;

    for(unsigned int kInd = 0; kInd <  methodPtr->numStages(); kInd ++ )
    {
        state_final = IC;

        for( int bInd = 0; bInd < kInd; bInd ++ )
            if(  methodPtr->a(kInd,bInd) != 0 )
                state_final += h* methodPtr->a(kInd,bInd)*K_list[bInd];

        double station  = start +  methodPtr->c(kInd)*h;


        state_constraint_function( station, state_final );

//        cout << state_final << endl;
        K_list[kInd] = ode_function(station, state_final );

    }

//    cout << "\t\t----------------------\t\t\n";
//    cout << Eigen::Vector3d::Zero() << endl;

    if( ! methodPtr->first_same_as_last() )
    {
        state_final = IC;
        for( int bInd = 0; (int) bInd <  methodPtr->numStages(); bInd ++ )
            if(  methodPtr->b(bInd) != 0)
                state_final += (h* methodPtr->b(bInd))*K_list[bInd];

        state_constraint_function( finish, state_final );

    }

    double error_est = -1;
    if( methodPtr->isAdaptive() )
    {
        // calculate lower order solution
        Eigen::MatrixXd state_final_comp = IC;
        for( int bInd = 0; (int) bInd <  methodPtr->numStages(); bInd ++ )
            if(  methodPtr->b_errEst(bInd) != 0)
                state_final_comp += (h* methodPtr->b_errEst(bInd))*K_list[bInd];

        // constrain per problem constraints
        state_constraint_function( finish, state_final_comp );

        // calculate error
        state_final_comp -= state_final;

        // find the max( abs( element ))
        error_est = state_final_comp.lpNorm<Eigen::Infinity>();
    }

    // Resort K_List if there are duplicates for the c's, later calculations of the slope at a station are assumed to be  more accurate
    for( int i=0; i< methodPtr->numStages(); i++ )
    {
        if( i !=  methodPtr->c_duplication_index(i))
            K_list[i] = K_list[ methodPtr->c_duplication_index(i)];
    }

    // update node parameters
    outputNode->t_o = start;
    outputNode->t_f = finish;

    // population interpolation matrix

    assert( methodPtr->order() > 0 ); // check to see if it is constraint based

    Eigen::MatrixXd interpMatrix =  methodPtr->buildInterpolationMatrix(h);

    assert( interpMatrix.rows() == methodPtr->order()+1 );
    assert( interpMatrix.cols() == 2+methodPtr->numStages() );

    outputNode->interpolationCoeff = std::vector<Eigen::MatrixXd >(methodPtr->order()+1, Eigen::MatrixXd::Zero(num_row,num_col));

    for( int coeffNum = 0; (int) coeffNum < outputNode->interpolationCoeff.size(); coeffNum ++ )
    {

        outputNode->interpolationCoeff[coeffNum] = interpMatrix(coeffNum,0)*IC;
        outputNode->interpolationCoeff[coeffNum] += interpMatrix(coeffNum,1+ methodPtr->numStages())*state_final;

        for( unsigned int Kind =0; Kind< methodPtr->numStages(); Kind++)
        {
            outputNode->interpolationCoeff[coeffNum] += interpMatrix(coeffNum,1+Kind)*K_list[Kind];
        }
    }

    if( error_est == -1 )
    {
        double s_prime = h;
        Eigen::MatrixXd state = Eigen::MatrixXd::Zero(num_row,num_col); // initialize the state to zero
        Eigen::MatrixXd state0 = Eigen::MatrixXd::Zero(num_row,num_col); // initialize the state to zero

        state0 = outputNode->interpolationCoeff[0];

        // Calculate the state coeff_i*(s_prime)^i
        for( unsigned int i = 0; i<= methodPtr->order(); i++ )
            state += outputNode->interpolationCoeff[i]*pow(s_prime,i);

        error_est = (IC-state0).lpNorm<Eigen::Infinity>()+(state_final - state).lpNorm<Eigen::Infinity>();

        // i dont trust this yet
        if( error_est == 0 )
            error_est = -1;
    }

    state_final_ = state_final;

    return error_est;
}

template< class T >
bool IVP_Solver<T>::solve(const Eigen::MatrixXd& initialCondition, double to, double tf, double tolerance, unsigned int maxNumSteps )
{

    solution_initial_state = initialCondition;


    if( methodPtr->type() == methodPtr->EXPLICIT )
    {

        double minStepSize =  (tf-to)/maxNumSteps;
        double stepSize = (tf-to)/maxNumSteps;

        double sign = (tf < to)?-1.0:1.0;
        double t_current = to;

        Eigen::MatrixXd IC = initialCondition;
        Eigen::MatrixXd FC = IC;
        solution_final_state = IC;
        //NodeCoeff* tmpNode;

        numActiveNodes = 0;

        stepSize = sign*std::min( sign*stepSize,sign*(tf-t_current) );
        //cout << IC << endl << endl;

        while( sign*(tf-t_current) > 0 )
        {
            //stepSize = sign*std::min( sign*stepSize,sign*(tf-t_current) );

            if( numActiveNodes == solution_node_list.size() )
                solution_node_list.push_back(0);

            if( solution_node_list[numActiveNodes] == 0 )
                solution_node_list[numActiveNodes] = new IVP_Solver<T>::NodeCoeff();

            // integrate node
            //cout << IC << endl << endl;
            double error = integrateNode(t_current, t_current+stepSize, IC, FC, solution_node_list[numActiveNodes] );
            assert( !std::isnan(error) );

            //IVP_Solver<T>::NodeCoeff* activeNodePointer = solution_node_list[numActiveNodes];
            //assert( activeNodePointer != 0 );
            //double error = integrateNode(t_current, t_current+stepSize, IC, solution_final_state, activeNodePointer );

            //            if( error > tolerance && sign*minStepSize < sign*stepSize )
            //            {
            //                //stepSize = sign*std::max( sign*stepSize/2.0, sign*minStepSize);
            //            } else
            //            {
            //                numActiveNodes ++;
            //                t_current += stepSize;
            //                IC = solution_final_state;

            //                //minStepSize=  (tf-t_current)/(maxNumSteps - numActiveNodes);

            //                /*if( error > 0 && error < tolerance/10.0 )
            //                {
            //                    stepSize *= 2.0;
            //                }*/
            //            }

            if( error <= tolerance || std::abs(stepSize) <= std::abs(minStepSize) )
            {
                numActiveNodes ++;
                t_current += stepSize;
                IC = FC;
            }


            assert( maxNumSteps >= numActiveNodes );
            if( maxNumSteps <= numActiveNodes )
            {
                minStepSize = tf-t_current;
            } else
            {
                minStepSize=  (tf-t_current)/(maxNumSteps - numActiveNodes); // maxNumSteps must always be > numActiveNodes to get here so were ok with the unsigned int subtraction
            }
            double suggestedStep = 0.9*stepSize*pow(tolerance/error,1.0/(methodPtr->order()+1.0));

            stepSize = sign*std::min( sign*(tf-t_current), std::max( sign*suggestedStep, sign*minStepSize) );

            //cout << "Error: " << error << "\tMin Step Size: " << minStepSize << "\tSuggested Step Size: " << suggestedStep << "\tFinal Step Size: " << stepSize << endl;
        }



        unsigned int nodeIndex = numActiveNodes;
        while( nodeIndex < solution_node_list.size() && solution_node_list[nodeIndex] )
        {
            solution_node_list[nodeIndex]->t_o = tf;
            solution_node_list[nodeIndex]->t_f = tf;
            nodeIndex++;
        }

        // update final state
        solution_final_state = FC;
        solutionFound = true;

    } else
    {
        std::cout << "ERROR: IVP_Solver<T>::solve for implicit functions not implemented yet" << std::endl;
        assert(false); // not implemented yet
    }

    //    cout << " ------------------------------ \n\n" << endl;
    //    cout << solution_final_state << endl;

    return isSolved();
}

template< class T >
bool IVP_Solver<T>::isSolved() const
{
    return solutionFound;
}

template< class T>
std::string IVP_Solver<T>::getSolverName() const
{
    return methodPtr->name();
}

template< class T >
Eigen::MatrixXd IVP_Solver<T>::getSolution() const
{
    /**< Returns the solution at the final time **/

    return solution_final_state;
}


template< class T >
Eigen::MatrixXd IVP_Solver<T>::getSolution(double time) const
{
    /**< Returns the solution at an intermediate time based on the interpolation polynomials **/

    double to = solution_node_list.front()->t_o;
    double tf = solution_node_list[numActiveNodes-1]->t_f;

    assert( inRange(to,tf,time) );

    if( time == to )
        return solution_initial_state;
    else if( time == tf )
        return solution_final_state;
    else
    {

        unsigned int nodeIndex = 0;
        //cout << "In Range: " << inRange(solution_node_list[nodeIndex]->t_o,solution_node_list[nodeIndex]->t_f,time) << std::endl;
        for(; nodeIndex < solution_node_list.size() && nodeIndex < numActiveNodes && !inRange(solution_node_list[nodeIndex]->t_o,solution_node_list[nodeIndex]->t_f,time); nodeIndex++ );

        assert( nodeIndex < solution_node_list.size() );

        double s_prime = time - solution_node_list[nodeIndex]->t_o;

        Eigen::MatrixXd state = solution_final_state;
        state.setZero(state.rows(),state.cols());

        // Calculate the state coeff_i*(s_prime)^i
        for( unsigned int i = 0; i<= methodPtr->order(); i++ )
            state += solution_node_list[nodeIndex]->interpolationCoeff[i]*pow(s_prime,i);

        return state;
    }
}

template< class T >
Eigen::MatrixXd IVP_Solver<T>::getSolutionIntegral() const
{
    /**< Returns the integral of the solution from the initial time to the final time based on the interpolation polynomials **/

    Eigen::MatrixXd stateIntegral = solution_final_state;
    stateIntegral.setZero(stateIntegral.rows(),stateIntegral.cols());

    int nodeIndex = 0;
    for(; nodeIndex < numActiveNodes && nodeIndex < solution_node_list.size(); nodeIndex++ )
    {
        double h = solution_node_list[nodeIndex]->t_f - solution_node_list[nodeIndex]->t_o;

        // Calculate the state coeff_i/i*(s_prime)^(i+1)
        for( unsigned int i = 0; i<= methodPtr->order(); i++ )
            stateIntegral += solution_node_list[nodeIndex]->interpolationCoeff[i]/((double)(i+1))* (pow(h,i+1) );//- pow(0,i+1));
    }

    return stateIntegral;
}

template< class T >
Eigen::MatrixXd IVP_Solver<T>::getSolutionIntegral(double to, double tf)
{
    /**< Returns the integration of the solution from t0 to tf based on the interpolation polynomials **/
    double soln_to = solution_node_list.front()->t_o;

    assert( numActiveNodes > 0 && isSolved() );
    double soln_tf = solution_node_list[numActiveNodes-1]->t_f;

    double sign = (tf < to)?-1.0:1.0;
    Eigen::MatrixXd stateIntegral = solution_final_state;
    stateIntegral.setZero(stateIntegral.rows(),stateIntegral.cols());

    if( sign < 0 )
    {
        // flip around
        double tmp = to;
        to = tf;
        tf = tmp;
    }


    if( soln_tf < soln_to )
    {
        std::cout << "ERROR: IVP_Solver<T>::getSolutionIntegral not implemented for negative initial integral direction " << std::endl;
        assert(false && ("IVP_Solver<T>::getSolutionIntegral not implemented for negative initial integral direction.")); // not implemented yet
    } else
    {
        int initialNodeIndex = 0;
        for(; initialNodeIndex < numActiveNodes && initialNodeIndex < solution_node_list.size() && !inRange(solution_node_list[initialNodeIndex]->t_o,solution_node_list[initialNodeIndex]->t_f,to); initialNodeIndex++ );


        int nodeIndex = initialNodeIndex;
        for( ; nodeIndex < numActiveNodes && nodeIndex < solution_node_list.size() && !inRange(solution_node_list[nodeIndex]->t_o,solution_node_list[nodeIndex]->t_f,tf); nodeIndex++ )
        {
            double tStart = std::max( solution_node_list[nodeIndex]->t_o, to );
            double tFinal = std::min( solution_node_list[nodeIndex]->t_f, tf );

            // Calculate the state coeff_i/i*(s_prime)^(i+1)
            for( unsigned int i = 0; i<= methodPtr->order(); i++ )
                stateIntegral += solution_node_list[nodeIndex]->interpolationCoeff[i]/((double)(i+1))* (pow(tFinal,i+1) - pow(tStart,i+1));
        }
    }


    return sign*stateIntegral;
}

template< class T >
bool IVP_Solver<T>::inRange(double bound1, double bound2, double point)
{
    // bool a = ((bound1 > bound2 ) && ( point >= bound2 && point <= bound1 ));
    //  bool b =  ((bound2 > bound1 ) && ( point >= bound1 && point <= bound2 ));
    return ((bound1 > bound2 ) && ( point >= bound2 && point <= bound1 ))
            || ((bound2 > bound1 ) && ( point >= bound1 && point <= bound2 ));
}

template< class T >
int IVP_Solver<T>::getNumberOfSteps() const
{
    if( !isSolved() )
        return -1;
    else
        return numActiveNodes;
}




#endif // IVP_SOLVER_H
