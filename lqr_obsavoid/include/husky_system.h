#ifndef HUSKYSYSTEM_H
#define HUSKYSYSTEM_H

#include "dynamic_system.h"

//X = [x y theta]
//U = [v w]
class HuskySystem:public Dynamic_System<3,2>//M = 3 states (x,y,theta),2 = 2 control inputs(v,w)
{

    double to;
    double tf;

public:

    /***** THESE FUNCTIONS SUPPORT CLOSED-LOOP ODE SOLUTION *******/
    virtual Matrix<double,3,1> x( double t) const; ///<- returns the state at time t
    virtual Matrix<double,3,1> x( ) const; ///<- returns the state at time t_f


    virtual bool solve(T_Span t_span, bool find_sensitivities = false, double tol = 1e-6, int maxSteps = 100); ///<- Solves the ODE

    virtual Matrix<double,3,3> Jx_sens() const; ///<- returns final value of dXf/dXo

    virtual Matrix<double,3,3> Jx_sens(double time) const; ///<- returns intermediate value of dXf/dXo


    virtual Matrix<double,3,2> Ju_sens() const; ///<- returns final value of dXf/dU

    virtual Matrix<double,3,2> Ju_sens(double time) const; ///<- returns intermediate value of dXf/dU


    virtual R3Tensor<3,3,3> Jxx_sens() const; ///<- returns final value of d^2Xf/dXo^2

    virtual R3Tensor<3,3,3> Jxx_sens(double time) const; ///<- returns intermediate value of d^2Xf/dXo^2

    virtual R3Tensor<3,2,2> Juu_sens() const; ///<- returns final value of d^2Xf/dU^2

    virtual R3Tensor<3,2,2> Juu_sens(double time) const; ///<- returns intermediate value of d^2Xf/dU^2

    virtual R3Tensor<3,3,2> Jxu_sens() const; ///<- returns final value of d^2Xf/dU^2

    virtual R3Tensor<3,3,2> Jxu_sens(double time) const; ///<- returns intermediate value of d^2Xf/dU^2


    /***** THESE FUNCTIONS SUPPORT NUMERIC ODE SOLUTION & General Dynamic System Stuff *******/
    virtual Matrix<double,3,1> g(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& ) const; ///<- returns the output y = X  Multiply by C for output selection
    virtual Matrix<double,3,1> ode(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& u) const; ///<- returns the time derivative of x, f(t,x,u)

    virtual Matrix<double,3,3> g_x_sens(double, const Matrix<double,3,1>& , const Matrix<double,2,1>& ) const; ///<- returns dg/dx
    //virtual Matrix<double,2,2> g_u_sens(double t, const Matrix<double,2,1>& x, const Matrix<double,2,1>& u) const; ///<- returns dg/du

    virtual Matrix<double,3,3> ode_x_sens(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& u) const; ///<- returns df/dx
    virtual Matrix<double,3,2> ode_u_sens(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& u) const; ///<- returns df/du

    //virtual R3Tensor<2,2,2> g_u2_sens(double t, const Matrix<double,2,1>& x, const Matrix<double,2,1>& u) const; ///<- returns d^2g/du^2
    //virtual R3Tensor<2,2,2> g_x2_sens(double t, const Matrix<double,2,1>& x, const Matrix<double,2,1>& u) const; ///<- returns d^2g/dx^2
    //virtual R3Tensor<2,2,2> g_ux_sens(double t, const Matrix<double,2,1>& x, const Matrix<double,2,1>& u) const; ///<- returns d^2g/dxu


    //virtual R3Tensor<3,2,2> ode_u2_sens(double, const Matrix<double,3,1>& , const Matrix<double,2,1>& u) const; ///<- returns d^2f/du^2
    virtual R3Tensor<3,3,3> ode_x2_sens(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& ) const; ///<- returns d^2f/dx^2

    virtual R3Tensor<3,3,2> ode_xu_sens(double, const Matrix<double,3,1>& , const Matrix<double,2,1>& ) const; ///<- returns d^2f/dxdu


    HuskySystem();
    virtual Dynamic_System<3,2>* clone() const;//cloning

};

#endif // HUSKYSYSTEM_H
