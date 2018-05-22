#ifndef LINEARSYSTEM_HPP
#define LINEARSYSTEM_HPP

#include "linearsystem.h"


template<int M, int N>
LinearSystem<M,N>::LinearSystem(const Eigen::Matrix<double,M,M>& A_, const Eigen::Matrix<double,M,N>& B_, const Eigen::Matrix<double,M,M>& C_, const Eigen::Matrix<double,M,N>& D_):
    Dynamic_System<M,N>(true,false,true,false,false,true,true,false,false,false,Eigen::Matrix<double,M,M>::Identity()),
    A(A_),
    B(B_),
    C(C_),
    D(D_),
    Ad(Eigen::Matrix<double,M,M>::Identity()),
    Bd(Eigen::Matrix<double,M,N>::Zero()),
    t_span(0,0)
{
    ;
}

template<int M, int N>
LinearSystem<M,N>::LinearSystem(double dt, const Eigen::Matrix<double,M,M>& A_, const Eigen::Matrix<double,M,N>& B_, const Eigen::Matrix<double,M,M>& C_, const Eigen::Matrix<double,M,N>& D_):
    Dynamic_System<M,N>(true,false,true,false,false,true,true,false,false,false,Eigen::Matrix<double,M,M>::Identity()),
    A(A_),
    B(B_),
    C(C_),
    D(D_),
    Ad(Eigen::Matrix<double,M,M>::Identity()),
    Bd(Eigen::Matrix<double,M,N>::Zero()),
    t_span(0,0)
{
    solve(T_Span(0,dt));
}


template<int M, int N>
 Dynamic_System<M,N>* LinearSystem<M,N>::clone() const
{
     LinearSystem<M,N>* newSys = new LinearSystem<M,N>(A,B,C,D);
     newSys->Ad = Ad;
     newSys->Bd = Bd;
     newSys->t_span = t_span;

    return static_cast<Dynamic_System<M,N>*>(newSys);
}

template<int M, int N>
 Matrix<double,M,1> LinearSystem<M,N>::g(double, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const
{
    ///<- returns the output y = g(t,x,u)  Multiply by C for output selection
    return C*x+D*u;
}

template<int M, int N>
Matrix<double,M,1> LinearSystem<M,N>::ode(double, const Matrix<double,M,1>& x, const Matrix<double,N,1>& u) const
{
    ///<- returns the time derivative of x, f(t,x,u)
    return A*x+B*u;
}

template<int M, int N>
Matrix<double,M,M> LinearSystem<M,N>::g_x_sens(double, const Matrix<double,M,1>& , const Matrix<double,N,1>& ) const
{
    ///<- returns dg/dx
    return C;
}


template<int M, int N>
Matrix<double,M,N> LinearSystem<M,N>::g_u_sens(double , const Matrix<double,M,1>& , const Matrix<double,N,1>& ) const
{
    ///<- returns dg/du
    return D;
}

template<int M, int N>
Matrix<double,M,M> LinearSystem<M,N>::ode_x_sens(double, const Matrix<double,M,1>&, const Matrix<double,N,1>& ) const
{
    ///<- returns df/dx
    return A;
}

template<int M, int N>
Matrix<double,M,N> LinearSystem<M,N>::ode_u_sens(double , const Matrix<double,M,1>& , const Matrix<double,N,1>& ) const
{   ///<- returns df/du
    return B;
}

template<int M, int N>
bool LinearSystem<M,N>::solve(T_Span t_span, bool, double, int)
{
    //<- Solves the ODE
    double dt = t_span.tf-t_span.t0;
    assert( dt != 0 && "Error in LinearSystem::solve: t_span must have a non-zero time span");

    if( std::abs(dt - this->t_span.tf+this->t_span.t0)/std::abs(dt) > 1e-8 )
    {
        Ad = (dt*A).exp();
        Bd = A.ColPivHouseholderQR.solve((Eigen::Matrix<double,M,M>::Identity-Ad)*B);
    }
    this->t_span = t_span;

    return true;
}

template<int M, int N>
Matrix<double,M,M> LinearSystem<M,N>::Jx_sens() const
{
    //<- returns final value of dXf/dXo
    return Ad;
}

template<int M, int N>
Matrix<double,M,M> LinearSystem<M,N>::Jx_sens(double time) const
{
    //<- returns intermediate value of dXf/dXo
    return (A*(time-t_span.t0)).exp();
}


template<int M, int N>
Matrix<double,M,N> LinearSystem<M,N>::Ju_sens() const
{
    //<- returns final value of dXf/dU
    return Bd;
}

template<int M, int N>
Matrix<double,M,N> LinearSystem<M,N>::Ju_sens(double time) const
{
    //<- returns intermediate value of dXf/dU
    double dt = time-t_span.t0;
    Eigen::Matrix<double,M,M> Ad_tmp = matrixExp(dt*A);
    Eigen::Matrix<double,M,N> Bd_tmp = A.ColPivHouseholderQR.solve((Eigen::Matrix<double,M,M>::Identity-Ad_tmp)*B);

    return Bd_tmp;
}

template<int M, int N>
Matrix<double,M,1> LinearSystem<M,N>::x( double t) const
{
    //<- returns the state at time t
    if( t == t_span.tf)
        return Ad*this->Xo+Bd*this->U;
    else
    {
        double dt = t-t_span.t0;
        Eigen::Matrix<double,M,M> Ad_tmp = matrixExp(dt*A);
        Eigen::Matrix<double,M,N> Bd_tmp = A.ColPivHouseholderQR.solve((Eigen::Matrix<double,M,M>::Identity-Ad_tmp)*B);

        return Ad_tmp*this->Xo+Bd_tmp*this->U;
    }

}

template<int M, int N>
Matrix<double,M,1> LinearSystem<M,N>::x( ) const
{
    return Ad*this->Xo+Bd*this->U;
}

template<int M, int N>
Matrix<double,M,M> LinearSystem<M,N>::matrixExp(const Matrix<double,M,M>& mat)
{
    Matrix<double,M,M> exp(Matrix<double,M,M>::Identity());
    Matrix<double,M,M> mat_pow_k( Matrix<double,M,M>::Identity());
    double k = 1;

    do
    {
        mat_pow_k *= mat/k;
        exp += mat_pow_k;
        k++;
    }while( mat_pow_k.norm()/exp.norm() > 1e-14 );

    return exp;
}



#endif

