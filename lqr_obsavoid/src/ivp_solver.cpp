#include "ivp_solver.h"

template< >
IVP_Solver<IVP_Solver_Default>::IVP_Solver(const IntegrationMethod& integrationMethod, Eigen::MatrixXd (*ode_function_ptr)( double time, const Eigen::MatrixXd& currentState ), Eigen::MatrixXd (*ode_Jacobian_function_prt)( double time, const Eigen::MatrixXd& currentState ), void (*state_constraint_function_prt)( double time, Eigen::MatrixXd& currentState ), void (*state_constraint_jacobian_prt)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian ) ):
    methodPtr(&integrationMethod),
    NR_solver(this, &IVP_Solver<IVP_Solver_Default>::NR_errorFunction, &IVP_Solver<IVP_Solver_Default>::NR_errorJacobian, EVALUATE, -1, -1),
    solution_node_list(500, (IVP_Solver<IVP_Solver_Default>::NodeCoeff*) 0 )
{
    //    solution_node_list = std::vector< typename IVP_Solver<IVP_Solver_Default>::NodeCoeff* >

    numActiveNodes = 0;
    solutionFound = false;

    this->owningObjectPtr = 0;
    this->ode_function_clptr = 0;
    this->ode_jacobian_function_clptr = 0;
    this->state_constraint_function_clptr = 0;
    this->state_constraint_jacobian_clptr = 0;


    this->ode_function_ptr = ode_function_ptr;
    this->ode_jacobian_function_ptr = ode_Jacobian_function_prt;
    this->state_constraint_function_ptr = state_constraint_function_prt;
    this->state_constraint_jacobian_ptr = state_constraint_jacobian_prt;

    assert( this->ode_function_ptr != 0 );
    solution_node_list.clear();
}

template< >
IVP_Solver<IVP_Solver_Default>::IVP_Solver(const IntegrationMethod& integrationMethod, IVP_Solver_Default& owningObject, Eigen::MatrixXd (IVP_Solver_Default::*ode_function_ptr)( double time, const Eigen::MatrixXd& currentState ), Eigen::MatrixXd (IVP_Solver_Default::*ode_Jacobian_function_prt)( double time, const Eigen::MatrixXd& currentState ), void (IVP_Solver_Default::*state_constraint_function_prt)( double time, Eigen::MatrixXd& currentState ), void (IVP_Solver_Default::*state_constraint_jacobian_prt)( double time, const Eigen::MatrixXd& currentState, Eigen::MatrixXd& constraintError, Eigen::MatrixXd& constraintJacobian ) ):
    methodPtr(&integrationMethod),
    NR_solver(this, &IVP_Solver<IVP_Solver_Default>::NR_errorFunction, &IVP_Solver<IVP_Solver_Default>::NR_errorJacobian, EVALUATE, -1, -1)
{
    solution_node_list.clear();
    assert( false ); // "Using the wrong constructor for IVP_Solver. You need to define the template type as your object type, eg. IVP_Solver<MyClassName>::mysolver(...)"
}
