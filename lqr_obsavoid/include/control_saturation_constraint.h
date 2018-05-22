#ifndef CONTROL_SATURATION_CONSTRAINT_H
#define CONTROL_SATURATION_CONSTRAINT_H

#include "dynamic_system_constraint.h"

template< int M, int N>
class Control_Saturation_Constraint: public Dynamic_System_Constraint<M,N>
{
    unsigned int controlChannel;
    double max_control;
    double min_control;
    double logistic_K;
    double max_threshold;
    double min_threshold;

public:
    Control_Saturation_Constraint(unsigned int channel, double max_control, double min_control);

    double constraint_cost(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const; ///<- The hard constraint. Should be zero when satisfied
    Eigen::Matrix<double,N,1> u_sens(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the derivative of the constraint with respect to control
    Eigen::Matrix<double,N,N> uu_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const;

    Eigen::Matrix<double,N,1> inforce_u_constraint(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u)const;

};

#include "control_saturation_constraint.hpp"

#endif // CONTROL_SATURATION_CONSTRAINT_H
