#ifndef DYNAMIC_SYSTEM_H
#define DYNAMIC_SYSTEM_H

#include <Eigen/Dense>
using Eigen::Matrix;

#include "ivp_solver.h"
#include "integrationmethod.h"
#include "r3tensor.h"
#include "tspan.h"

template<int M, int N >///<- M is number of states, N is number of controls
class Dynamic_System
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // need to do this because we dynamically allociate this object and it has a fixed sized Eigen Matrix

    Dynamic_System(bool non_zero_out_u_sens, bool non_zero_out_u2_sens, bool non_zero_out_x_sens, bool non_zero_out_x2_sens,  bool non_zero_out_ux_sens, bool non_zero_ode_x_sens, bool non_zero_ode_u_sens,bool non_zero_ode_xx_sense, bool non_zero_ode_xu_sens, bool non_zero_ode_uu_sens, Eigen::Matrix<double,Eigen::Dynamic,M> outputSelect = Eigen::Matrix<double,M,M>::Identity(), const IntegrationMethod& ode_method = Dormand_Prince_8th, const IntegrationMethod& jacob_method = Rk_10_8_Curtis );
    virtual ~Dynamic_System() {;}

    static int getNumStates() {return M;}
    static int getNumControls() {return N;}
    int getNumOutputs() const;

    Matrix<double,Eigen::Dynamic,M> getOutSelectMatrix() const; ///<- Returns the output selection matrix C the output y = C*g(x) This should be an identity matrix of sorts
    void setOutSelectMatrix(Matrix<double,Eigen::Dynamic,M> C); ///<- sets the output selection matrix C the output y = C*g(x) This should be an identity matrix of sorts

    virtual Matrix<double,M,1> y( double t) const; ///<- returns the output quantity g(t,x,u) at time t you need to right multiply by C for the outputs for this system.
    virtual Matrix<double,M,1> y( ) const;         ///<- returns the output quantity g(t,x,u) at time t_f

    virtual Matrix<double,M,1> x( double t) const; ///<- returns the state at time t
    virtual Matrix<double,M,1> x( ) const;         ///<- returns the state at time t_f



    void set_Xo(const Matrix<double,M,1>& Xo_); ///<- Sets the initial condition
    void set_U(const Matrix<double,N,1>& U_);   ///<- Sets the control input


    bool isSolved() const;
    virtual bool solve(T_Span t_span, bool find_sensitivities = false, double tol = 1e-6, int maxSteps = 100); ///<- Solves the ODE
    bool solve(const Matrix<double,M,1>& Xo_, const Matrix<double,N,1>& U_, T_Span t_span, bool find_sensitivities = false, double tol = 1e-6, int maxSteps = 100); ///<- sets new IC and U and solves the ODE

    bool has_out_x_sense() const;  ///<- Indicates if the output function y = g(t,x,u) has a nontrivial dy/dx
    bool has_out_u_sense() const;  ///<- Indicates if the output function y = g(t,x,u) has a nontrivial dy/du
    bool has_out_uu_sense() const; ///<- Indicates if the output function y = g(t,x,u) has a nontrivial d^2y/du^2
    bool has_out_xx_sense() const; ///<- Indicates if the output finction y = g(t,x,u) has a nontrivial d^2y/dx^2
    bool has_out_ux_sense() const; ///<- Indicates if the output function y = g(t,x,u) has a nontrivial d^2y/dudx

    bool has_ode_x_sense() const;   ///<- Indicates if Jx is non trivial
    bool has_ode_u_sense() const;   ///<- Indicates if Ju is non trivial
    bool has_ode_xx_sense() const; ///<- Indicates if Jxx is non trivial
    bool has_ode_xu_sense() const; ///<- Indicates if Jxu is non trivial
    bool has_ode_uu_sense() const; ///<- Indicates if Juu is non trivial

    virtual Matrix<double,M,M> Jx_sens() const; ///<- returns final value of dXf/dXo
    virtual Matrix<double,M,M> Jx_sens(double time) const; ///<- returns intermediate value of dXf/dXo


    virtual Matrix<double,M,N> Ju_sens() const; ///<- returns final value of dXf/dU
    virtual Matrix<double,M,N> Ju_sens(double time) const; ///<- returns intermediate value of dXf/dU


    virtual R3Tensor<M,M,M> Jxx_sens() const; ///<- returns final value of d^2Xf/dXo^2
    virtual R3Tensor<M,M,M> Jxx_sens(double time) const; ///<- returns intermediate value of d^2Xf/dXo^2

    virtual R3Tensor<M,N,N> Juu_sens() const; ///<- returns final value of d^2Xf/dU^2
    virtual R3Tensor<M,N,N> Juu_sens(double time) const; ///<- returns intermediate value of d^2Xf/dU^2

    virtual R3Tensor<M,M,N> Jxu_sens() const; ///<- returns final value of d^2Xf/dU^2
    virtual R3Tensor<M,M,N> Jxu_sens(double time) const; ///<- returns intermediate value of d^2Xf/dU^2


    virtual Matrix<double,M,M> Yx_sens() const; ///<- returns final value of dg/dx
    virtual Matrix<double,M,M> Yx_sens(double time) const ; ///<- returns intermediate value of dg/dx

    virtual Matrix<double,M,N> Yu_sens() const; ///<- returns final value of dg/du
    virtual Matrix<double,M,N> Yu_sens(double time) const;  ///<- returns intermediate value of dg/du


    virtual R3Tensor<M,M,M> Yxx_sens() const; ///<- returns final value of d^2g/dx^2
    virtual R3Tensor<M,M,M> Yxx_sens(double time) const;  ///<- returns intermediate value of d^2g/dx^2

    virtual R3Tensor<M,N,N> Yuu_sens() const; ///<- returns final value of d^2g/du^2
    virtual R3Tensor<M,N,N> Yuu_sens(double time) const;   ///<- returns intermediate value of d^2g/du^2

    virtual R3Tensor<M,N,M> Yux_sens() const; ///<- returns final value of d^2g/dxdu
    virtual R3Tensor<M,N,M> Yux_sens(double time) const;  ///<- returns intermediate value of d^2g/dxdu


/***********  VIRTUAL FUNCTIONS AND PARAMETERS TO BE IMPLEMENTED BY INHERTING CLASS  **************/
public:
    /**
     * @brief clone returns a pointer to a newly allocated dyamic system object to be implemented in inheriting class
     * @return Dynamic_System*
     *
     * If a class FOO inheritis Dyanmic system then this funciton should consist of:  return static_cast<Dynamic_System*>(new FOO());
     * This helps with the trajectory class and making new trajectory nodes.
     *
     *
     * \see trajectory.h
     *
     */
    virtual Dynamic_System<M,N>* clone() const = 0;

    bool verify_jacobians(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u);
    bool verify_jacobians();


    virtual Matrix<double,M,1> g(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const = 0; ///<- returns the output y = g(t,x,u)  Multiply by C for output selection
    virtual Matrix<double,M,M> g_x_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns dg/dx
    virtual Matrix<double,M,N> g_u_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns dg/du; default is zeros
    virtual R3Tensor<M,N,N>    g_u2_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns d^2g/du^2; default is zeros
    virtual R3Tensor<M,M,M>    g_x2_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns d^2g/dx^2; default is zeros
    virtual R3Tensor<M,N,M>    g_ux_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns d^2g/dx^2; default is zeros


protected:
    virtual Matrix<double,M,1> ode(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const = 0; ///<- returns the time derivative of x, f(t,x,u)

    virtual Matrix<double,M,M> ode_x_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const = 0; ///<- returns df/dx
    virtual Matrix<double,M,N> ode_u_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const = 0; ///<- returns df/du

    virtual R3Tensor<M,N,N> ode_u2_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns d^2f/du^2; default is zeros
    virtual R3Tensor<M,M,M> ode_x2_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns d^2f/dx^2; default is zeros

    virtual R3Tensor<M,M,N> ode_xu_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const; ///<- returns d^2f/dudx; default is zeros

/***** END VIRTUAL FUNCTIONS *****/

protected:  /* protected members owned by this class */
    bool out_x_sens;  ///<- Indicates if the output function y = g(t,x,u) has a nontrivial dy/dx
    bool out_u_sens;  ///<- Indicates if the output function y = g(t,x,u) has a nontrivial dy/du
    bool out_u2_sens; ///<- Indicates if the output function y = g(t,x,u) has a nontrivial d^2y/du^2
    bool out_x2_sens; ///<- Indicates if the output finction y = g(t,x,u) has a nontrivial d^2y/dx^2
    bool out_ux_sens; ///<- Indicates if the output function y = g(t,x,u) has a nontrivial d^2y/dxdu
    bool g_x_sens_bool;    ///<- Indicates if the ode has a nontrivial df/dx
    bool g_u_sens_bool;    ///<- Indicates if the ode has a nontrivial df/du
    bool g_xx_sens;   ///<- Indicates if the ode has a nontrivial d^2f/dx^2
    bool g_xu_sens;   ///<- Indicates if the ode has a nontrivial d^2f/dxdu
    bool g_uu_sens;   ///<- Indicates if the ode has a nontrivial d^2f/dudu

    bool solved;

    T_Span solTSpan;

    Eigen::Matrix<double,M,1> Xo; ///<- The intial condition
    Eigen::Matrix<double,N,1> U;  ///<- The Control Input
    Eigen::Matrix<double,Dynamic,M> C; ///<- The output selection matrix

    IVP_Solver< Dynamic_System<M,N> > x_int;
    IVP_Solver< Dynamic_System<M,N> > x_sens_int;
    IVP_Solver< Dynamic_System<M,N> > x2_sens_int;
    IVP_Solver< Dynamic_System<M,N> > u_sens_int;
    IVP_Solver< Dynamic_System<M,N> > u2_sens_int;
    IVP_Solver< Dynamic_System<M,N> > xu_sens_int;



 private:
    /* These are helper functions for the IVP solver */
    MatrixXd    ode_x_dot      (double t, const MatrixXd& x_cur );
    MatrixXd    ode_x_sens_dot (double t, const MatrixXd& Jx_cur);
    MatrixXd    ode_u_sens_dot (double t, const MatrixXd& Ju_cur);
    MatrixXd    ode_uu_sens_dot(double t, const MatrixXd& Juu_cur);
    MatrixXd    ode_xx_sens_dot(double t, const MatrixXd& Jxx_cur);
    MatrixXd    ode_xu_sens_dot(double t, const MatrixXd& Jxu_cur);


};


// include the template code
// if you have template specializations they go in a seporate cpp
#include "dynamic_system.hpp"


#endif // DYNAMIC_SYSTEM_H
