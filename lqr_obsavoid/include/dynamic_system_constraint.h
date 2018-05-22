#ifndef DYNAMIC_SYSTEM_CONSTRAINT
#define DYNAMIC_SYSTEM_CONSTRAINT

#include <Eigen/Core>

enum Dynamic_System_Constraint_Type
{
    SOFT, ///<- positive semidefinate cost function such as a potential field
    HARD  ///<- Legrange Multiplier enforced hard constraint
};

template<int M, int N>
class Dynamic_System_Constraint
{
public:
    Dynamic_System_Constraint(bool has_x_sens_, bool has_xx_sens_, bool has_ux_sens_, bool has_u_sens_, bool has_uu_sens_);
    bool has_x_sens() const; ///<- Returns if the constraint has a non-zero d/dx sensitivity
    bool has_xx_sens() const; ///<- Returns if the constraint has a non-zero d^2/dxdx sensitivity
    bool has_u_sens() const; ///<- Returns if the constraint has a non-zero d/du sensitivity
    bool has_uu_sens() const; ///<- Returns if the constraint has a non-zero d^2/dudu sensitivity
    bool has_ux_sens() const; ///<- Returns if the constraint has a non-zero d^2/dxdu sensitivity

    virtual double constraint_cost(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const = 0; ///<- The hard constraint. Should be zero when satisfied
    virtual Eigen::Matrix<double,M,1> x_sens(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the derivative of the constraint with respect to state
    virtual Eigen::Matrix<double,N,1> u_sens(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the derivative of the constraint with respect to control
    virtual Eigen::Matrix<double,M,M> xx_sens(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the second derivative of the constraint with respect to state^2
    virtual Eigen::Matrix<double,N,N> uu_sens(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the derivative of the constraint with respect to control^2
    virtual Eigen::Matrix<double,N,M> ux_sens(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the derivative of the constraint with respect to control then state

    virtual Eigen::Matrix<double,N,1> inforce_u_constraint(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u)const;

    bool verify_constraint_jacobians(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type);
protected:
    bool b_u_sens;
    bool b_x_sens;
    bool b_uu_sens;
    bool b_xx_sens;
    bool b_ux_sens;

    double logistic(double x, double k=100) const;
    double logisticDer(double x, double k=100) const;
    double logisitcDer2(double x, double k=100) const;

    static const double logisitcDer2_kx_zero = -2.3993572805154676678327396972823;

};


#include "dynamic_system_constraint.hpp"

#endif // DYNAMIC_SYSTEM_CONSTRAINT

