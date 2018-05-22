#ifndef DYNAMIC_SYSTEM_HPP
#define DYNAMIC_SYSTEM_HPP

#include "dynamic_system.h"

/* Template code associated with dynamic_system.h definitions */


template<int M, int N>
Dynamic_System<M,N>::Dynamic_System(bool non_zero_out_u_sens, bool non_zero_out_u2_sens, bool non_zero_out_x_sens, bool non_zero_out_x2_sens,  bool non_zero_out_ux_sens, bool non_zero_ode_x_sens, bool non_zero_ode_u_sens,bool non_zero_ode_xx_sense, bool non_zero_ode_xu_sens, bool non_zero_ode_uu_sens, Eigen::Matrix<double,Eigen::Dynamic,M> outputSelect, const IntegrationMethod& ode_method, const IntegrationMethod& jacob_method ):
    out_u_sens(non_zero_out_u_sens),
    out_u2_sens(non_zero_out_u2_sens),
    out_x_sens(non_zero_out_x_sens),
    out_x2_sens(non_zero_out_x2_sens),
    out_ux_sens(non_zero_out_ux_sens),
    g_x_sens_bool(non_zero_ode_x_sens),
    g_u_sens_bool(non_zero_ode_u_sens),
    g_xx_sens(non_zero_ode_xx_sense),
    g_xu_sens(non_zero_ode_xu_sens),
    g_uu_sens(non_zero_ode_uu_sens),
    solved(false),
    C(outputSelect),
    x_int       (ode_method,   *this, &Dynamic_System<M,N>::ode_x_dot),
    x_sens_int  (jacob_method, *this, &Dynamic_System<M,N>::ode_x_sens_dot),
    x2_sens_int (jacob_method, *this, &Dynamic_System<M,N>::ode_xx_sens_dot),
    u_sens_int  (jacob_method, *this, &Dynamic_System<M,N>::ode_u_sens_dot),
    u2_sens_int (jacob_method, *this, &Dynamic_System<M,N>::ode_uu_sens_dot),
    xu_sens_int (jacob_method, *this, &Dynamic_System<M,N>::ode_xu_sens_dot)
{
    if( !g_x_sens_bool )
        assert( !g_xx_sens && !g_xu_sens && "Cannot have xx sensitivity or xu sensitivity without x sensitivity");

    if( !g_u_sens_bool )
        assert( !g_uu_sens && !g_xu_sens && "Cannot have uu sensitivity or xu sensitivity without u sensitivity");

    if( !out_u_sens )
        assert( !out_u2_sens && !out_ux_sens && "Cannot have uu sensitivity or xu sensitivity without u sensitivity");

    if( !out_x_sens )
        assert( !out_x2_sens && !out_ux_sens && "Cannot have xx sensitivity or xu sensitivity without x sensitivity");

    return;
}

template<int M, int N>
void Dynamic_System<M,N>::set_Xo(const Matrix<double,M,1>& Xo_)
{
    ///<- Sets the initial condition
    Xo = Xo_;
    solved = false;
}


template<int M, int N>
void Dynamic_System<M,N>::set_U(const Matrix<double,N,1>& U_)
{
    ///<- Sets the control input
    U = U_;
    solved = false;
}

template<int M, int N>
int Dynamic_System<M,N>::getNumOutputs() const
{
    return C.rows();
}

template<int M, int N>
Matrix<double,Eigen::Dynamic,M> Dynamic_System<M,N>::getOutSelectMatrix() const
{
    ///<- Returns the output selection matrix C the output y = C*g(x) This should be an identity matrix of sorts
    return C;
}

template<int M, int N>
void Dynamic_System<M,N>::setOutSelectMatrix(Matrix<double,Eigen::Dynamic,M> newOutputSelectionMat)
{
    ///<- sets the output selection matrix C the output y = C*g(x) This should be an identity matrix of sorts
    assert(newOutputSelectionMat.rows() <= M && "Output Selection matrix may not have more rows than the number of states");
    C = newOutputSelectionMat;
}

template<int M, int N>
bool  Dynamic_System<M,N>::isSolved() const
{
    return solved;
}

template<int M, int N>
bool  Dynamic_System<M,N>::solve(T_Span t_span, bool find_sensitivities, double tol, int maxSteps )
{

    solved = true;
    solTSpan = t_span;

    solved &= x_int.solve(Xo,t_span.t0,t_span.tf,tol,maxSteps);

    if( find_sensitivities )
    {
        if( has_ode_x_sense() )
        {
            Eigen::Matrix<double,M,M> Jx_0(Eigen::Matrix<double,M,M>::Identity());
            solved &= x_sens_int.solve(Jx_0,t_span.t0,t_span.tf,tol,maxSteps);
        }

        if( has_ode_u_sense() )
        {
            Eigen::Matrix<double,M,N> Ju_0(Eigen::Matrix<double,M,N>::Zero());
            solved &= u_sens_int.solve(Ju_0,t_span.t0,t_span.tf,tol,maxSteps);
        }

        if( has_ode_uu_sense() )
        {
            solved &= u2_sens_int.solve(R3Tensor<M,N,N>::Zero().toMatrixXd(),t_span.t0,t_span.tf,tol,maxSteps);
        }

        if( has_ode_xx_sense() )
        {
            solved &= x2_sens_int.solve(R3Tensor<M,M,M>::Zero().toMatrixXd(),t_span.t0,t_span.tf,tol,maxSteps);
        }

        if( has_ode_xx_sense() )
        {
            solved &= xu_sens_int.solve(R3Tensor<M,M,N>::Zero().toMatrixXd(),t_span.t0,t_span.tf,tol,maxSteps);
        }

    }

    return solved;
}

template<int M, int N>
bool  Dynamic_System<M,N>::solve(const Matrix<double,M,1>& Xo_, const Matrix<double,N,1>& U_, T_Span t_span, bool find_sensitivities, double tol, int maxSteps )
{ ///<- sets new IC and U and solves the ODE

    set_Xo(Xo_);
    set_U(U_);

    return solve(t_span,find_sensitivities, tol, maxSteps);

}




template<int M, int N>
MatrixXd    Dynamic_System<M,N>::ode_x_dot(double t, const MatrixXd& x_cur )
{
    assert( x_cur.rows() == M && "Dynamic_System::ode_x_dot failed because x_cur is the wrong size");
    return ode(t,x_cur,U);
}

template<int M, int N>
MatrixXd    Dynamic_System<M,N>::ode_x_sens_dot (double t, const MatrixXd& Jx_cur)
{
    assert( Jx_cur.rows() == M && Jx_cur.cols()==M && "Dynamic_System::ode_x_sens_dot failed because Jx_cur is the wrong size");
    assert( x_int.isSolved() && "Dynamic_System::ode_x_sens_dot failed because x has not be solved for");

    Eigen::Matrix<double,M,1> x_cur = x_int.getSolution(t);

    // Jx_dot = df/dx*Jx(t)
    return ode_x_sens(t,x_cur,U)*Jx_cur;
}

template<int M, int N>
MatrixXd    Dynamic_System<M,N>::ode_u_sens_dot (double t, const MatrixXd& Ju_cur)
{
    assert( Ju_cur.rows() == M && Ju_cur.cols() == N && "Dynamic_System::ode_u_sens_dot failed because Ju_cur is the wrong size");
    assert( x_int.isSolved() && "Dynamic_System::ode_u_sens_dot failed because x has not be solved for");

    Eigen::Matrix<double,M,1> x_cur = x_int.getSolution(t);

    // Ju_dot = df/dx*Ju + df/du
    return ode_x_sens(t,x_cur,U)*Ju_cur + ode_u_sens(t,x_cur,U);
}

template<int M, int N>
MatrixXd    Dynamic_System<M,N>::ode_uu_sens_dot(double t, const MatrixXd& Juu_cur)
{
    // nonzero if ode_xu or ode_uu or ode_xx are non zero

    assert( Juu_cur.rows() == M*N && Juu_cur.cols() == N && "Dynamic_System::ode_uu_sens_dot failed because Juu_cur is the wrong size");
    assert( x_int.isSolved() && "Dynamic_System::ode_uu_sens_dot failed because x has not be solved for");
    assert( u_sens_int.isSolved() && "Dynamic_System::ode_uu_sens_dot failed because Ju has not be solved for");

    R3Tensor<M,N,N> Juu_tens(Juu_cur);

    Eigen::Matrix<double,M,1> x_cur = x_int.getSolution(t);
    Matrix<double,M,N> Ju(u_sens_int.getSolution(t));


    // Juu_dot = df/dxx*Ju Ju + df/dxu Ju + df/duu + df/dxJuu
    R3Tensor<M,N,N> Juu_dot(ode_x_sens(t,x_cur,U)*Juu_tens);

    if( g_uu_sens )
        Juu_dot += ode_u2_sens(t,x_cur,U);

    if( g_xx_sens )
        Juu_dot += ode_x2_sens(t,x_cur,U).depthMultiply(Ju)*Ju;

    if( g_xu_sens )
    {
        R3Tensor<M,N,N> odeXU_Ju = ode_xu_sens(t,x_cur,U)*Ju;
        R3Tensor<M,N,N> odeUX_Ju;
        for( int i=0; i< N; i++ )
        {
            odeUX_Ju.rd_slice(i) = odeXU_Ju.rc_slice(i);
        }
        //Juu_dot = odeUXsens.depthMultiply(Ju) + ode_xu_sens(t,x_cur,U)*Ju;
        Juu_dot += odeXU_Ju + odeUX_Ju;
    }


    // cout << Juu_dot + odeUXsens.depthMultiply(Ju) << endl << endl;
    return Juu_dot.toMatrixXd();
}

template<int M, int N>
MatrixXd    Dynamic_System<M,N>::ode_xx_sens_dot(double t, const MatrixXd& Jxx_cur)
{
    // nonzero if xx sensitivity

    assert( Jxx_cur.rows() == M*M && Jxx_cur.cols() == M && "Dynamic_System::ode_xx_sens_dot failed because Jxx_cur is the wrong size");
    assert( x_int.isSolved() && "Dynamic_System::ode_xx_sens_dot failed because x has not be solved for");
    assert( x_sens_int.isSolved() && "Dynamic_System::ode_xx_sens_dot failed because Jx has not be solved for");


    Eigen::Matrix<double,M,1> x_cur = x_int.getSolution(t);
    Eigen::Matrix<double,M,M> Jx_cur = x_sens_int.getSolution(t);


    // Jxx_dot = d^fdx^2 Jx Jx + df_dx*Jxx_cur
    R3Tensor<M,M,M> Jxx_cur_tens(Jxx_cur);
    R3Tensor<M,M,M> df_dxx(ode_x2_sens(t,x_cur,U));

    return (df_dxx.depthMultiplyInPlace(Jx_cur)*Jx_cur + ode_x_sens(t,x_cur,U)*Jxx_cur_tens).toMatrixXd();
}

template<int M, int N>
MatrixXd    Dynamic_System<M,N>::ode_xu_sens_dot(double t, const MatrixXd& Jxu_cur)
{
    // nonzero if xu or xx sensitivity

    assert( Jxu_cur.rows() == M*N && Jxu_cur.cols() == M && "Dynamic_System::ode_xu_sens_dot failed because Jxu_cur is the wrong size");
    assert( x_int.isSolved() && "Dynamic_System::ode_xu_sens_dot failed because x has not be solved for");
    assert( x_sens_int.isSolved() && "Dynamic_System::ode_xu_sens_dot failed because Jx has not be solved for");
    assert( u_sens_int.isSolved() && "Dynamic_System::ode_xu_sens_dot failed because Ju has not be solved for");

    Eigen::Matrix<double,M,1> x_cur = x_int.getSolution(t);
    Eigen::Matrix<double,M,M> Jx_cur = x_sens_int.getSolution(t);
    Eigen::Matrix<double,M,N> Ju_cur = u_sens_int.getSolution(t);

    R3Tensor<M,M,N> Jxu_tens(Jxu_cur);
    R3Tensor<M,M,N> df_dxu;
    df_dxu.setZero();

    if( g_xu_sens )
        df_dxu += ode_xu_sens(t,x_cur,U);
    if( g_xx_sens )
        df_dxu += ode_x2_sens(t,x_cur,U).depthMultiply(Ju_cur);


    //Jxu_dot = d2f_dxx*Jx*Ju_in_depth + dfdx*Jxu + d2f_dxdu*Jx
    return (df_dxu*Jx_cur + ode_x_sens(t,x_cur,U)*Jxu_tens ).toMatrixXd();
}

template<int M, int N>
Matrix<double,M,1>      Dynamic_System<M,N>::y(double t) const
{
    return g(t,x(t),U);
}

template<int M, int N>
Matrix<double,M,1>      Dynamic_System<M,N>::x(double t) const
{
    return x_int.getSolution(t);
}

template<int M, int N>
Matrix<double,M,1>    Dynamic_System<M,N>::y() const
{
    return g(solTSpan.tf,x(),U);
}

template<int M, int N>
Matrix<double,M,1>      Dynamic_System<M,N>::x() const
{
    return x_int.getSolution();
}


template<int M, int N>
bool        Dynamic_System<M,N>::has_out_x_sense() const
{
    //<- Indicates if the output function y = g(t,x,u) has a nontrivial dy/dx
    return out_x_sens;
}


template<int M, int N>
bool        Dynamic_System<M,N>::has_out_u_sense() const
{
    return out_u_sens;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_out_uu_sense() const
{
    return out_u2_sens;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_out_xx_sense() const
{
    return out_x2_sens;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_out_ux_sense() const
{
    return out_ux_sens;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_ode_x_sense() const
{
    return g_x_sens_bool;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_ode_u_sense() const
{
    return g_u_sens_bool;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_ode_xx_sense() const
{
    return g_xx_sens;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_ode_xu_sense() const
{
    return g_xu_sens || g_xx_sens;
}

template<int M, int N>
bool        Dynamic_System<M,N>::has_ode_uu_sense() const
{
    return g_uu_sens || g_xx_sens || g_xu_sens;
}

template<int M, int N>
Matrix<double,M,M>      Dynamic_System<M,N>::Jx_sens() const
{
    if( has_ode_x_sense() )
        return x_sens_int.getSolution();
    else
        return Matrix<double,M,M>::Identity();
}

template<int M, int N>
Matrix<double,M,M>      Dynamic_System<M,N>::Jx_sens(double time) const
{
    if( has_ode_x_sense() )
        return x_sens_int.getSolution(time);
    else
        return Matrix<double,M,M>::Identity();
}

template<int M, int N>
Matrix<double,M,N>      Dynamic_System<M,N>::Ju_sens() const
{
    if( has_ode_u_sense() )
        return u_sens_int.getSolution();
    else
        return Matrix<double,M,N>::Zero();
}

template<int M, int N>
Matrix<double,M,N>      Dynamic_System<M,N>::Ju_sens(double time) const
{
    if( has_ode_u_sense() )
        return u_sens_int.getSolution(time);
    else
        return Matrix<double,M,N>::Zero();
}

template<int M, int N>
R3Tensor<M,M,M>    Dynamic_System<M,N>::Jxx_sens() const
{
    if( has_ode_xx_sense() )
        return R3Tensor<M,M,M>(x2_sens_int.getSolution());
    else
        return R3Tensor<M,M,M>::Zero();
}

template<int M, int N>
R3Tensor<M,M,M>    Dynamic_System<M,N>::Jxx_sens(double time) const
{
    if( has_ode_xx_sense() )
        return R3Tensor<M,M,M>(x2_sens_int.getSolution(time));
    else
        return R3Tensor<M,M,M>::Zero();
}

template<int M, int N>
R3Tensor<M,M,N> Dynamic_System<M,N>::Jxu_sens() const
{
    ///<- returns final value of d^2Xf/dU^2
    if( has_ode_xu_sense() )
        return R3Tensor<M,M,N>(xu_sens_int.getSolution());
    else
        return R3Tensor<M,M,N>::Zero();
}


template<int M, int N>
R3Tensor<M,M,N> Dynamic_System<M,N>::Jxu_sens(double time) const
{
    ///<- returns intermediate value of d^2Xf/dU^2

    if( has_ode_xu_sense() )
        return R3Tensor<M,M,N>(xu_sens_int.getSolution(time));
    else
        return R3Tensor<M,M,N>::Zero();
}


template<int M, int N>
R3Tensor<M,N,N>    Dynamic_System<M,N>::Juu_sens() const
{
    if( has_ode_uu_sense() )
        return R3Tensor<M,N,N>(u2_sens_int.getSolution());
    else
        return R3Tensor<M,N,N>::Zero();
}

template<int M, int N>
R3Tensor<M,N,N>    Dynamic_System<M,N>::Juu_sens(double time) const
{
    if( has_ode_uu_sense() )
        return R3Tensor<M,N,N>(u2_sens_int.getSolution(time));
    else
        return R3Tensor<M,N,N>::Zero();
}

template<int M, int N>
Matrix<double,M,M>      Dynamic_System<M,N>::Yx_sens() const
{
    return g_x_sens(solTSpan.tf,x(),U);
}

template<int M, int N>
Matrix<double,M,M>      Dynamic_System<M,N>::Yx_sens(double time) const
{
    return g_x_sens(time,x(time),U);
}

template<int M, int N>
Matrix<double,M,N>      Dynamic_System<M,N>::Yu_sens() const
{
    if( has_out_u_sense() )
        return g_u_sens(solTSpan.tf,x(),U);
    else
        return Matrix<double,M,N>::Zero();
}

template<int M, int N>
Matrix<double,M,N>      Dynamic_System<M,N>::Yu_sens(double time) const
{
    if( has_out_u_sense() )
        return g_u_sens(time,x(time),U);
    else
        return Matrix<double,M,N>::Zero();
}

template<int M, int N>
R3Tensor<M,M,M>    Dynamic_System<M,N>::Yxx_sens() const
{
    if( has_out_xx_sense() )
        return g_x2_sens(solTSpan.tf,x(),U);
    else
        return R3Tensor<M,M,M>::Zero();
}

template<int M, int N>
R3Tensor<M,M,M>    Dynamic_System<M,N>::Yxx_sens(double time) const
{
    if( has_out_xx_sense() )
        return  g_x2_sens(time,x(time),U);
    else
        return R3Tensor<M,M,M>::Zero();
}

template<int M, int N>
R3Tensor<M,N,N>    Dynamic_System<M,N>::Yuu_sens() const
{
    if( has_out_uu_sense() )
        return  g_u2_sens(solTSpan.tf,x(),U);
    else
        return R3Tensor<M,N,N>::Zero();
}

template<int M, int N>
R3Tensor<M,N,N>    Dynamic_System<M,N>::Yuu_sens(double time) const
{
    if( has_out_uu_sense() )
        return g_u2_sens(time, x(time), U);
    else
        return R3Tensor<M,N,N>::Zero();
}

template<int M, int N>
R3Tensor<M,N,M> Dynamic_System<M,N>::Yux_sens() const
{
    //<- returns final value of d^2g/dxdu
    if( has_out_ux_sense() )
        return g_ux_sens(solTSpan.tf, x(), U);
    else
        return R3Tensor<M,N,M>::Zero();
}

template<int M, int N>
R3Tensor<M,N,M> Dynamic_System<M,N>::Yux_sens(double time) const
{
    //<- returns intermediate value of d^2g/dxdu
    if( has_out_ux_sense() )
        return g_ux_sens(time, x(time), U);
    else
        return R3Tensor<M,N,M>::Zero();
}


template<int M, int N>
Matrix<double,M,N> Dynamic_System<M,N>::g_u_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>&) const
{
    //<- returns dg/du; default is zeros
    assert(!out_u_sens && "g_u_sens Undefined for inheriting class that says it should exist.");
    return Matrix<double,M,N>::Zero();
}

template<int M, int N>
Matrix<double,M,M> Dynamic_System<M,N>::g_x_sens(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const
{
    assert(!out_x_sens && "g_x_sens not defined for inheriting class that says it should exist.");
    return Matrix<double,M,M>::Zero();
}

template<int M, int N>
R3Tensor<M,N,N> Dynamic_System<M,N>::g_u2_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>& ) const
{
    //<- returns d^2g/du^2; default is zeros
    assert(!out_u2_sens && "g_u2_sens Undefined for inheriting class that says it shoudl exist.");
    return   R3Tensor<M,N,N>::Zero();
}

template<int M, int N>
R3Tensor<M,M,M> Dynamic_System<M,N>::g_x2_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>& ) const
{
    //<- returns d^2g/dx^2; default is zeros
    assert(!out_x2_sens && "g_x2_sens Undefined for inheriting class that says it shoudl exist.");
    return   R3Tensor<M,M,M>::Zero();
}

template<int M, int N>
R3Tensor<M,N,M> Dynamic_System<M,N>::g_ux_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>& ) const
{
    //<- returns d^2g/dx^2; default is zeros
    assert(!out_ux_sens && "g_xu_sens Undefined for inheriting class that says it shoudl exist.");
    return   R3Tensor<M,N,M>::Zero();
}

template<int M, int N>
R3Tensor<M,N,N> Dynamic_System<M,N>::ode_u2_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>& ) const
{

    //<- returns d^2f/du^2; default is zeros
    assert(!g_uu_sens && "ode_u2_sens Undefined for inheriting class that says it shoudl exist.");
    return   R3Tensor<M,N,N>::Zero();
}

template<int M, int N>
R3Tensor<M,M,M> Dynamic_System<M,N>::ode_x2_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>& ) const
{

    //<- returns d^2f/dx^2; default is zeros
    assert(!g_xx_sens && "ode_x2_sens Undefined for inheriting class that says it shoudl exist.");
    return   R3Tensor<M,M,M>::Zero();
}

template<int M, int N>
R3Tensor<M,M,N> Dynamic_System<M,N>::ode_xu_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>& ) const
{
    //<- returns d^2f/dudx; default is zeros
    assert(!g_xu_sens && "ode_xu_sens Undefined for inheriting class that says it shoudl exist.");
    return   R3Tensor<M,M,N>::Zero();
}

template<int M, int N>
bool Dynamic_System<M,N>::verify_jacobians()
{
    Eigen::Matrix<double,M,1> X_save = Xo;
    T_Span tSpanSave = solTSpan;
    Eigen::Matrix<double,N,1> U_save = U;

    bool retval = verify_jacobians(solTSpan.tf-solTSpan.t0, X_save, U_save);

    Xo = X_save;
    U = U_save;
    solTSpan = tSpanSave;

    solve(solTSpan,true);

    return retval;
}


template<int M, int N>
bool Dynamic_System<M,N>::verify_jacobians(double t, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u)
{
    double intagration_tol = 1e-11;
    Matrix<double,M,1> g_nom;
    Matrix<double,M,1> ode_nom;

    Matrix<double,M,M> g_x_sens_nom, g_x_sens_num;
    Matrix<double,M,N> g_u_sens_nom, g_u_sens_num;

    Matrix<double,M,M> ode_x_sens_nom, ode_x_sens_num;
    Matrix<double,M,N> ode_u_sens_nom, ode_u_sens_num;

    R3Tensor<M,N,N> g_u2_sens_nom, g_u2_sens_num;
    R3Tensor<M,M,M> g_x2_sens_nom, g_x2_sens_num;
    R3Tensor<M,N,M> g_ux_sens_nom, g_ux_sens_num;

    R3Tensor<M,N,N> ode_u2_sens_nom, ode_u2_sens_num;
    R3Tensor<M,M,M> ode_x2_sens_nom, ode_x2_sens_num;

    R3Tensor<M,M,N> ode_xu_sens_nom, ode_xu_sens_num;

    g_nom =  g(t, x, u); ///<- returns the output y = g(t,x,u)  Multiply by C for output selection
    ode_nom = ode(t, x, u); ///<- returns the time derivative of x, f(t,x,u)

    solve(x,u,T_Span(0,t),true,intagration_tol,10000);
    Matrix<double,M,1> x_nom = this->x();
    Matrix<double,M,1> y_nom = y();
    Matrix<double,M,M> Jx_nom,Jx_num,Yx_nom,Yx_num;
    Matrix<double,M,N> Ju_nom,Ju_num,Yu_nom,Yu_num;
    R3Tensor<M,M,M> Jxx_nom,Jxx_num,Yxx_nom,Yxx_num;
    R3Tensor<M,N,N> Juu_nom,Juu_num,Yuu_nom,Yuu_num;
    R3Tensor<M,M,N> Jxu_nom,Jxu_num;
    R3Tensor<M,N,M> Yux_nom,Yux_num;

    Jx_nom = Jx_sens();
    Yx_nom = Yx_sens();
    Ju_nom = Ju_sens();

    if( has_out_u_sense() )
        Yu_nom = Yu_sens();
    else
    {
        Yu_nom.setZero();
        Yux_num.setZero();
    }

    if( has_ode_xx_sense() )
        Jxx_nom = Jxx_sens();
    else
        Jxx_nom.setZero();

    if( has_ode_uu_sense() )
        Juu_nom = Juu_sens();
    else
        Juu_nom.setZero();

    if( has_ode_xu_sense() )
        Jxu_nom = Jxu_sens();
    else Jxu_nom.setZero();

    if( has_out_uu_sense() )
        Yuu_nom = Yuu_sens();
    else
        Yuu_nom.setZero();

    if( has_out_ux_sense() )
        Yux_nom = Yux_sens();
    else
        Yux_nom.setZero();

    if( has_out_xx_sense())
        Yxx_nom = Yxx_sens();
    else
        Yxx_nom.setZero();



    /*
    bool out_u_sens;  ///<- Indicates if the output function y = g(t,x,u) has a nontrivial dy/du
    bool out_u2_sens; ///<- Indicates if the output function y = g(t,x,u) has a nontrivial d^2y/du^2
    bool out_x2_sens; ///<- Indicates if the output finction y = g(t,x,u) has a nontrivial d^2y/dx^2
    bool out_ux_sens; ///<- Indicates if the output function y = g(t,x,u) has a nontrivial d^2y/dxdu
    bool g_xx_sens; ///<- Indicates if the ode has a nontrivial d^2f/dx^2
    bool g_xu_sens; ///<- Indicates if the ode has a nontrivial d^2f/dxdu
    bool g_uu_sens;
    */

    g_x_sens_nom = g_x_sens(t,x,u); ///<- returns dg/dx

    if( out_u_sens )
        g_u_sens_nom = g_u_sens(t, x, u); ///<- returns dg/du; default is zeros
    else
    {
        g_u_sens_nom.setZero();
        g_u2_sens_num.setZero();
    }

    ode_x_sens_nom = ode_x_sens(t, x, u); ///<- returns df/dx
    ode_u_sens_nom = ode_u_sens(t, x, u); ///<- returns df/du

    if( out_u2_sens )
        g_u2_sens_nom =  g_u2_sens(t, x, u); ///<- returns d^2g/du^2; default is zeros
    else
        g_u2_sens_nom.setZero();

    if( out_x2_sens )
        g_x2_sens_nom =  g_x2_sens(t, x, u); ///<- returns d^2g/dx^2; default is zeros
    else
        g_x2_sens_nom.setZero();

    if( out_ux_sens )
        g_ux_sens_nom =  g_ux_sens(t, x, u); ///<- returns d^2g/dx^2; default is zeros
    else
        g_ux_sens_nom.setZero();

    if( g_uu_sens )
        ode_u2_sens_nom = ode_u2_sens( t, x, u); ///<- returns d^2f/du^2; default is zeros
    else
        ode_u2_sens_nom.setZero();

    if( g_xx_sens )
        ode_x2_sens_nom = ode_x2_sens( t, x, u); ///<- returns d^2f/dx^2; default is zeros
    else
        ode_x2_sens_nom.setZero();

    if( g_xu_sens )
        ode_xu_sens_nom = ode_xu_sens( t, x, u); ///<- returns d^2f/dudx; default is zeros
    else
        ode_xu_sens_nom.setZero();

    double delta = 1e-6;
    Matrix<double,M,1> dx = Matrix<double,M,1>::Zero();
    Matrix<double,N,1> du = Matrix<double,N,1>::Zero();

    for( int i=0; i<M; i++ )
    {
        dx(i) = std::max(std::abs(x(i))*delta,delta);
        g_x_sens_num.block(0,i,M,1) = (g(t,x+dx,u)-g(t,x-dx,u))/(2.0*dx(i));
        g_x2_sens_num.rc_slice(i) = (g_x_sens(t,x+dx,u)-g_x_sens(t,x-dx,u))/(2.0*dx(i));

        ode_x_sens_num.block(0,i,M,1) = (ode(t, x+dx, u)-ode(t, x-dx, u))/(2.0*dx(i)); ///<- returns df/dx
        ode_x2_sens_num.rc_slice(i) = ( ode_x_sens(t,x+dx,u) - ode_x_sens(t,x-dx,u))/(2.0*dx(i));
        ode_xu_sens_num.rd_slice(i) = ( ode_u_sens(t,x+dx,u) - ode_u_sens(t,x-dx,u))/(2.0*dx(i));

        Eigen::Matrix<double,M,1> x_shift = x+dx;
        solve(x_shift,u,T_Span(0,t),true,intagration_tol,10000);
        Jx_num.block(0,i,M,1) = this->x();
        Yx_num.block(0,i,M,1) = this->y();
        Jxx_num.rc_slice(i) = this->Jx_sens();
        Yxx_num.rc_slice(i) = this->Yx_sens();
        if( has_out_u_sense() )
            Yux_num.rc_slice(i) = Yu_sens();

        x_shift = x-dx;
        solve(x_shift,u,T_Span(0,t),true,intagration_tol,10000);
        Jx_num.block(0,i,M,1) -= this->x();
        Yx_num.block(0,i,M,1) -= this->y();
        Jxx_num.rc_slice(i) -= this->Jx_sens();
        Yxx_num.rc_slice(i) -= this->Yx_sens();
        if( has_out_u_sense() )
            Yux_num.rc_slice(i) -= Yu_sens();

        Jx_num.block(0,i,M,1) /= 2.0*dx(i);
        Yx_num.block(0,i,M,1) /= 2.0*dx(i);
        Jxx_num.rc_slice(i) /= 2.0*dx(i);
        Yxx_num.rc_slice(i) /= 2.0*dx(i);
        if( has_out_u_sense() )
            Yux_num.rc_slice(i) /= 2.0*dx(i);

        dx(i) = 0;
    }

    for( int i=0; i<N; i++ )
    {
        du(i) = std::max(std::abs(u(i))*delta,delta);

        g_u_sens_num.block(0,i,M,1) = (g(t,x,u+du)-g(t,x,u-du))/(2.0*du(i));

        if(out_u_sens)
            g_u2_sens_num.rc_slice(i) = (g_u_sens(t,x,u+du)-g_u_sens(t,x,u-du))/(2.0*du(i));

        g_ux_sens_num.rd_slice(i) = (g_x_sens(t,x,u+du)-g_x_sens(t,x,u-du))/(2.0*du(i));

        ode_u_sens_num.block(0,i,M,1) = (ode(t, x, u+du)-ode(t, x, u-du))/(2.0*du(i)); ///<- returns df/dx
        ode_u2_sens_num.rc_slice(i) = (ode_u_sens(t,x,u+du) - ode_u_sens(t,x,u-du))/(2.0*du(i));


        solve(x,u+du,T_Span(0,t),true,intagration_tol,10000);
        Ju_num.block(0,i,M,1) = this->x();
        Yu_num.block(0,i,M,1) = y();
        Juu_num.rc_slice(i) = Ju_sens();
        Yuu_num.rc_slice(i) = Yu_sens();
        Jxu_num.rc_slice(i) = Jx_sens();

        solve(x,u-du,T_Span(0,t),true,intagration_tol,10000);
        Ju_num.block(0,i,M,1) -= this->x();
        Yu_num.block(0,i,M,1) -= y();
        Juu_num.rc_slice(i) -= Ju_sens();
        Yuu_num.rc_slice(i) -= Yu_sens();
        Jxu_num.rc_slice(i) -= Jx_sens();

        Ju_num.block(0,i,M,1) /= 2.0*du(i);
        Yu_num.block(0,i,M,1) /= 2.0*du(i);
        Juu_num.rc_slice(i) /= 2.0*du(i);
        Yuu_num.rc_slice(i) /= 2.0*du(i);
        Jxu_num.rc_slice(i) /= 2.0*du(i);

        du(i) = 0;
    }

    double g_x_sens_err = (g_x_sens_nom - g_x_sens_num).norm()/g_x_sens_nom.norm();
    double g_u_sens_err = (g_u_sens_nom - g_u_sens_num).norm()/g_u_sens_nom.norm();
    double ode_x_sens_err = (ode_x_sens_nom - ode_x_sens_num).norm()/ode_x_sens_nom.norm();
    double ode_u_sens_err = (ode_u_sens_nom - ode_u_sens_num).norm()/ode_u_sens_nom.norm();
    double g_u2_sens_err = (g_u2_sens_nom - g_u2_sens_num).norm()/g_u2_sens_nom.norm();
    double g_x2_sens_err = (g_x2_sens_nom - g_x2_sens_num).norm()/g_x2_sens_nom.norm();
    double g_ux_sens_err = (g_ux_sens_nom - g_ux_sens_num).norm()/g_ux_sens_nom.norm();
    double ode_u2_sens_err = (ode_u2_sens_nom - ode_u2_sens_num).norm()/ode_u2_sens_nom.norm();
    double ode_x2_sens_err = (ode_x2_sens_nom - ode_x2_sens_num).norm()/ode_x2_sens_nom.norm();
    double ode_xu_sens_err = (ode_xu_sens_nom - ode_xu_sens_num).norm()/ode_xu_sens_nom.norm();

    double Jx_err =  (Jx_nom-Jx_num).norm()/Jx_nom.norm();
    double Jxx_err = (Jxx_nom-Jxx_num).norm()/Jxx_nom.norm();
    double Jxu_err = (Jxu_nom-Jxu_num).norm()/Jxu_nom.norm();
    double Ju_err =  (Ju_nom-Ju_num).norm()/Ju_nom.norm();
    double Juu_err = (Juu_nom-Juu_num).norm()/Juu_nom.norm();

    /*
    double Yx_err =  (Yx_nom-Yx_num).norm()/Yx_nom.norm();
    double Yxx_err = (Yxx_nom-Yxx_num).norm()/Yxx_nom.norm();
    double Yux_err = (Yux_nom-Yux_num).norm()/Yux_nom.norm();
    double Yu_err =  (Yu_nom-Yu_num).norm()/Yu_nom.norm();
    double Yuu_err = (Yuu_nom-Yuu_num).norm()/Yuu_nom.norm();
    */

    bool retVal = true;

    double tol = 1e-4;

    if( g_x_sens_err > tol )
    {
        retVal = false;
        cout << "g_x_sens err: " << g_x_sens_err * 100 << "%" << endl;
        cout << "g_x_sens nom: " << endl << g_x_sens_nom << endl;
        cout << "g_x_sens num: " << endl << g_x_sens_num << endl;
        cout << "g_x_sens nom-num: " << endl << g_x_sens_nom-g_x_sens_num << endl << endl;
    }

    if( g_u_sens_err > tol )
    {
        retVal = false;
        cout << "g_u_sens err: " << g_u_sens_err * 100 << "%" << endl;
        cout << "g_u_sens nom: " << endl << g_u_sens_nom << endl;
        cout << "g_u_sens_num: " << endl << g_u_sens_num << endl;
        cout << "g_u_sens nom-num: " << endl << g_u_sens_nom-g_u_sens_num << endl << endl;

    }

    if( g_u2_sens_err > tol )
    {
        retVal = false;
        cout << "g_u2_sens err: " << g_u2_sens_err * 100 << "%" << endl;
        cout << "g_u2_sens nom: " << endl << g_u2_sens_nom << endl;
        cout << "g_u2_sens_num: " << endl << g_u2_sens_num << endl;
        cout << "g_u2_sens nom-num: " << endl << g_u2_sens_nom-g_u2_sens_num << endl << endl;

    }

    if( g_ux_sens_err > tol )
    {
        retVal = false;
        cout << "g_ux_sens err: " << g_ux_sens_err * 100 << "%" << endl;
        cout << "g_ux_sens nom: " << endl << g_ux_sens_nom << endl;
        cout << "g_ux_sens_num: " << endl << g_ux_sens_num << endl;
        cout << "g_ux_sens nom-num: " << endl << g_ux_sens_nom-g_ux_sens_num << endl << endl;

    }

    if( g_x2_sens_err > tol )
    {
        retVal = false;
        cout << "g_x2_sens err: " << g_x2_sens_err * 100 << "%" << endl;
        cout << "g_x2_sens nom: " << endl << g_x2_sens_nom << endl;
        cout << "g_x2_sens_num: " << endl << g_x2_sens_num << endl;
        cout << "g_x2_sens nom-num: " << endl << g_x2_sens_nom-g_x2_sens_num << endl << endl;

    }



    if( ode_x_sens_err > tol )
    {
        retVal = false;
        cout << "ode_x_sens err: " << ode_x_sens_err * 100 << "%" << endl;
        cout << "ode_x_sens nom: " << endl << ode_x_sens_nom << endl;
        cout << "ode_x_sens_num: " << endl << ode_x_sens_num << endl;
        cout << "ode_x_sens nom-num: " << endl << ode_x_sens_nom-ode_x_sens_num << endl << endl;

    }

    if( ode_u_sens_err > tol )
    {
        retVal = false;
        cout << "ode_u_sens err: " << ode_u_sens_err * 100 << "%" << endl;
        cout << "ode_u_sens nom: " << endl << ode_u_sens_nom << endl;
        cout << "ode_u_sens_num: " << endl << ode_u_sens_num << endl;
        cout << "ode_u_sens nom-num: " << endl << ode_u_sens_nom-ode_u_sens_num << endl << endl;

    }

    if( ode_u2_sens_err > tol )
    {
        retVal = false;
        cout << "ode_u2_sens err: " << ode_u2_sens_err * 100 << "%" << endl;
        cout << "ode_u2_sens nom: " << endl << ode_u2_sens_nom << endl;
        cout << "ode_u2_sens_num: " << endl << ode_u2_sens_num << endl;
        cout << "ode_u2_sens nom-num: " << endl << ode_u2_sens_nom-ode_u2_sens_num << endl << endl;

    }

    if( ode_xu_sens_err > tol )
    {
        retVal = false;
        cout << "ode_xu_sens err: " << ode_xu_sens_err * 100 << "%" << endl;
        cout << "ode_xu_sens nom: " << endl << ode_xu_sens_nom << endl;
        cout << "ode_xu_sens_num: " << endl << ode_xu_sens_num << endl;
        cout << "ode_xu_sens nom-num: " << endl << ode_xu_sens_nom-ode_xu_sens_num << endl << endl;

    }

    if( ode_x2_sens_err > tol )
    {
        retVal = false;
        cout << "ode_x2_sens err: " << ode_x2_sens_err * 100 << "%" << endl;
        cout << "ode_x2_sens nom: " << endl << ode_x2_sens_nom << endl;
        cout << "ode_x2_sens_num: " << endl << ode_x2_sens_num << endl;
        cout << "ode_x2_sens nom-num: " << endl << ode_x2_sens_nom-ode_x2_sens_num << endl << endl;

    }


    if( Jx_err > tol )
    {
        retVal = false;
        cout << "Jx err: " << Jx_err * 100 << "%" << endl;
        cout << "Jx nom: " << endl << Jx_nom << endl;
        cout << "Jx num: " << endl << Jx_num << endl;
        cout << "Jx nom-num: " << endl << Jx_nom-Jx_num << endl << endl;

    }
    if( Jxx_err > tol )
    {
        retVal = false;
        cout << "Jxx err: " << Jxx_err * 100 << "%" << endl;
        cout << "Jxx nom: " << endl << Jxx_nom << endl;
        cout << "Jxx num: " << endl << Jxx_num << endl;
        cout << "Jxx nom-num: " << endl << Jxx_nom-Jxx_num << endl << endl;

    }
    if( Ju_err > tol )
    {
        retVal = false;
        cout << "Ju err: " << Ju_err * 100 << "%" << endl;
        cout << "Ju nom: " << endl << Ju_nom << endl;
        cout << "Ju num: " << endl << Ju_num << endl;
        cout << "Ju nom-num: " << endl << Ju_nom-Ju_num << endl << endl;

    }
    if( Juu_err > tol )
    {
        retVal = false;
        cout << "Juu err: " << Juu_err * 100 << "%" << endl;
        cout << "Juu nom: " << endl << Juu_nom << endl;
        cout << "Juu num: " << endl << Juu_num << endl;
        cout << "Juu nom-num: " << endl << Juu_nom-Juu_num << endl << endl;

    }
    if( Jxu_err > tol )
    {
        retVal = false;
        cout << "Jxu err: " << Jxu_err * 100 << "%" << endl;
        cout << "Jxu nom: " << endl << Jxu_nom << endl;
        cout << "Jxu num: " << endl << Jxu_num << endl;
        cout << "Jxu nom-num: " << endl << Jxu_nom-Jxu_num << endl << endl;

    }


    /*
    if( Yx_err > tol )
    {
        retVal = false;
        cout << "Yx err: " << Yx_err * 100 << "%" << endl;
        cout << "Yx nom: " << endl << Yx_nom << endl;
        cout << "Yx num: " << endl << Yx_num << endl << endl;
    }
    if( Yxx_err > tol )
    {
        retVal = false;
        cout << "Yxx err: " << Yxx_err * 100 << "%" << endl;
        cout << "Yxx nom: " << endl << Yxx_nom << endl;
        cout << "Yxx num: " << endl << Yxx_num << endl << endl;
    }
    if( Yu_err > tol )
    {
        retVal = false;
        cout << "Yu err: " << Yu_err * 100 << "%" << endl;
        cout << "Yu nom: " << endl << Yu_nom << endl;
        cout << "Yu num: " << endl << Yu_num << endl << endl;
    }
    if( Yuu_err > tol )
    {
        retVal = false;
        cout << "Yuu err: " << Yuu_err * 100 << "%" << endl;
        cout << "Yuu nom: " << endl << Yuu_nom << endl;
        cout << "Yuu num: " << endl << Yuu_num << endl << endl;
    }
    if( Yux_err > tol )
    {
        retVal = false;
        cout << "Yux err: " << Yux_err * 100 << "%" << endl;
        cout << "Yux nom: " << endl << Yux_nom << endl;
        cout << "Yux num: " << endl << Yux_num << endl << endl;
    }
    */
    return retVal;

}

#endif
