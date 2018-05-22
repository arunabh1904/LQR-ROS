#ifndef Control_Saturation_Constraint_HPP
#define Control_Saturation_Constraint_HPP

#include "control_saturation_constraint.h"

template<int M, int N>
Control_Saturation_Constraint<M,N>::Control_Saturation_Constraint(unsigned int channel, double max_control_, double min_control_):
    Dynamic_System_Constraint<M,N>(false,false,false,true,true)
{
    assert( channel < N && "Control_Saturation_Constraint must specify a control channel in range for the system");
    assert( max_control_ > min_control_ && "Control_Saturation_Constraint must specify a constraint range where max_control is greator than min_control");

    controlChannel = channel;
    max_control = max_control_;
    min_control = min_control_;

    logistic_K = 100.0;

    max_threshold = 0.5*this->logisitcDer2_kx_zero/logistic_K+max_control;
    min_threshold = 0.5*this->logisitcDer2_kx_zero/(-logistic_K)+min_control;
}


template<int M, int N>
double Control_Saturation_Constraint<M,N>::constraint_cost(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type) const
{
    //<- The hard constraint. Should be zero when satisfied
    double retVal = 0;

    if( u(controlChannel) >= max_threshold )
        retVal = this->logistic(u(controlChannel) - max_control,logistic_K);
    else if( u(controlChannel) <= min_threshold )
        retVal = this->logistic(u(controlChannel) - min_control,-logistic_K);

    return retVal;

}


template<int M, int N>
Eigen::Matrix<double,N,1> Control_Saturation_Constraint<M,N>::u_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type) const
{
    //<- the derivative of the constraint with respect to control
    Eigen::Matrix<double,N,1> retSens(Eigen::Matrix<double,N,1>::Zero());

    if( u(controlChannel) >= max_threshold )
        retSens(controlChannel) = this->logisticDer(u(controlChannel) - max_control,logistic_K);
    else if( u(controlChannel) <= min_threshold )
        retSens(controlChannel) = this->logisticDer(u(controlChannel) - min_control,-logistic_K);

    return retSens;
}

template<int M, int N>
Eigen::Matrix<double,N,N> Control_Saturation_Constraint<M,N>::uu_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type) const
{
    //<- the derivative of the constraint with respect to control
    Eigen::Matrix<double,N,N> retSens(Eigen::Matrix<double,N,N>::Zero());

    if( u(controlChannel) >= max_threshold )
        retSens(controlChannel,controlChannel) = this->logisitcDer2(u(controlChannel) - max_control,logistic_K);
    else if( u(controlChannel) <= min_threshold )
        retSens(controlChannel,controlChannel) = this->logisitcDer2(u(controlChannel) - min_control,-logistic_K);

    return retSens;
}

template<int M, int N>
Eigen::Matrix<double,N,1> Control_Saturation_Constraint<M,N>::inforce_u_constraint(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>& u)const
{
    Eigen::Matrix<double,N,1> u_new = u;

    if( u_new(controlChannel) > max_threshold )
        u_new(controlChannel) =  this->logisitcDer2_kx_zero/logistic_K  + max_control;
    else if( u_new(controlChannel) < min_threshold )
        u_new(controlChannel) = this->logisitcDer2_kx_zero/(-logistic_K) + min_control;

    return u_new;
}


#endif
