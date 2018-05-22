#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <Eigen/Dense>
#include "dynamic_system.h"
#include "linearsystem.h"
#include "tspan.h"

template< int M, int N >
class Trajectory
{
public:
    Eigen::VectorXd target_output; ///<- Desired Target Location
    Eigen::VectorXd target_control; ///<- Nominal Control Input to get there 0 if we dont know yet
    Eigen::MatrixXd Q; ///<-  Target Error Weighting matrix.  1/max_target_error^2 along diagonal is a good starting point
    Eigen::MatrixXd R; ///<-  Control Penalty Matrix.  1/max_control_des^2 along diagonal is a good starting point
    bool includeTargetControlInCost; ///<- Allows us to decide if control cost is (target + U)'*R*(target + U) or just U'*R*U

    Eigen::VectorXd soln_state; ///<- Solved location of system at target time. Some error may/will exist!
    Eigen::VectorXd soln_output; ///<- Solved system output of system at target time. Some error relative to target will exist!
    Eigen::VectorXd soln_control; ///<- Delta control input control = nominal_control + soln_control

    T_Span t_span; ///<- Time span for this node

    Dynamic_System<M,N>* p_dynamic_system; ///<- Pointer to a dynamic system object



    Trajectory( const Dynamic_System<M,N>& system );
    Trajectory( Dynamic_System<M,N>* system );

    Trajectory( );
    ~Trajectory( );

//    Trajectory operator=(const Trajectory& other )
//    {
//        if( p_dynamic_system )
//            delete p_dynamic_system;


//        p_dynamic_system = other.p_dynamic_system->clone();
//        p_dynamic_system->setOutSelectMatrix(other.p_dynamic_system->getOutSelectMatrix());
//        target_output = other.target_output;
//        target_control = other.target_control;
//        includeTargetControlInCost  = other.includeTargetControlInCost;
//        Q = other.Q;
//        R = other.R;
//        t_span = other.t_span;

//        soln_state = other.soln_state;
//        soln_control = other.soln_control;
//        soln_output = other.soln_output;
//    }

    Trajectory(const Trajectory& other ); ///<- Creates a trajectory from another coying all parameters but making a new dynamic system pointer

    void clone(const Trajectory& other); ///<- Creates a trajectory from another coying all parameters but making a new dynamic system pointer
    Trajectory clone(); ///<- Creates a trajectory from another coying all parameters but making a new dynamic system pointer
};

template< int M, int N>
Trajectory<M,N> Trajectory<M,N>::clone()
{
    return Trajectory<M,N>(*this);
}

template< int M, int N>
void Trajectory<M,N>::clone(const Trajectory<M,N>& other)
{
    p_dynamic_system = other.p_dynamic_system->clone();
    p_dynamic_system->setOutSelectMatrix(other.p_dynamic_system->getOutSelectMatrix());
    target_output = other.target_output;
    target_control = other.target_control;
    includeTargetControlInCost  = other.includeTargetControlInCost;
    Q = other.Q;
    R = other.R;
    t_span = other.t_span;

    soln_state = other.soln_state;
    soln_control = other.soln_control;
    soln_output = other.soln_output;
}

template< int M, int N>
Trajectory<M,N>::Trajectory( const Dynamic_System<M,N>& system )
{
    p_dynamic_system = system.clone();
}

template< int M, int N>
Trajectory<M,N>::Trajectory( Dynamic_System<M,N>* system )
{
  p_dynamic_system = system;
}

template< int M, int N>
Trajectory<M,N>::Trajectory( )
{
    p_dynamic_system = 0;
}

template< int M, int N>
Trajectory<M,N>::~Trajectory<M,N>()
{
    if( p_dynamic_system )
        delete p_dynamic_system;

    p_dynamic_system = 0;
}

template< int M, int N >
Trajectory<M,N>::Trajectory(const Trajectory& other )
{
    p_dynamic_system = other.p_dynamic_system->clone();
    p_dynamic_system->setOutSelectMatrix(other.p_dynamic_system->getOutSelectMatrix());
    target_output = other.target_output;
    target_control = other.target_control;
    includeTargetControlInCost  = other.includeTargetControlInCost;
    Q = other.Q;
    R = other.R;
    t_span = other.t_span;

    soln_state = other.soln_state;
    soln_control = other.soln_control;
    soln_output = other.soln_output;
}




template< int M, int N >
Trajectory<M,N> makeInitialCondition(Eigen::Matrix<double,M,1> Xo)
{
    Trajectory<M,N>  retval( LinearSystem<M,N>(Eigen::Matrix<double,M,M>::Zeros(), Eigen::Matrix<double,M,N>::Zeros() ));
    retval.target_control = Eigen::Matrix<double,N,1>::Zeros();
    retval.target_output = Xo;
    retval.t_span = T_Span(0,0);
    retval.Q = Eigen::Matrix<double,M,M>::Zeros();
    retval.R = Eigen::Matrix<double,N,N>::Zeros();
    retval.soln_state = Xo;
    retval.soln_control = Eigen::Matrix<double,N,1>::Zeros();

    return retval;
}

#endif // TRAJECTORY_H

