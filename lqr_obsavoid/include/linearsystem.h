#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "dynamic_system.h"
#include "tspan.h"

template<int M, int N>
class LinearSystem:public Dynamic_System<M,N>
{

public:
    LinearSystem(const Eigen::Matrix<double,M,M>& A_, const Eigen::Matrix<double,M,N>& B_, const Eigen::Matrix<double,M,M>& C_ = Eigen::Matrix<double,M,M>::Identity(), const Eigen::Matrix<double,M,N>& D_ = Eigen::Matrix<double,M,N>::Zero());
    LinearSystem(double dT, const Eigen::Matrix<double,M,M>& A_, const Eigen::Matrix<double,M,N>& B_, const Eigen::Matrix<double,M,M>& C_ = Eigen::Matrix<double,M,M>::Identity(), const Eigen::Matrix<double,M,N>& D_ = Eigen::Matrix<double,M,N>::Zero());

    virtual Dynamic_System<M,N>* clone() const;

    bool solve(T_Span t_span, bool, double, int); ///<- Solves the ODE

    Matrix<double,M,1> x( double t) const; ///<- returns the state at time t
    Matrix<double,M,1> x( ) const;

    Matrix<double,M,M> Jx_sens() const; ///<- returns final value of dXf/dXo
    Matrix<double,M,M> Jx_sens(double time) const; ///<- returns intermediate value of dXf/dXo

    Matrix<double,M,N> Ju_sens() const; ///<- returns final value of dXf/dU
    Matrix<double,M,N> Ju_sens(double time) const; ///<- returns intermediate value of dXf/dU



protected:
    virtual Matrix<double,M,1> g(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns the output y = g(t,x,u)  Multiply by C for output selection
    virtual Matrix<double,M,1> ode(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns the time derivative of x, f(t,x,u)

    virtual Matrix<double,M,M> g_x_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns dg/dx
    virtual Matrix<double,M,N> g_u_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns dg/du

    virtual Matrix<double,M,M> ode_x_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns df/dx
    virtual Matrix<double,M,N> ode_u_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns df/du

private:
    Eigen::Matrix<double,M,M> A;
    Eigen::Matrix<double,M,N> B;
    Eigen::Matrix<double,M,M> C;
    Eigen::Matrix<double,M,N> D;

    Eigen::Matrix<double,M,M> Ad;
    Eigen::Matrix<double,M,N> Bd;

    T_Span t_span;

    Matrix<double,M,M> matrixExp(const Matrix<double,M,M>& mat);

};

#include "linearsystem.hpp"

#endif // LINEARSYSTEM_H
