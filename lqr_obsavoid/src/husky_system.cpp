#include "husky_system.h"
#include "math.h"
#include "r3tensor.h"
using namespace Eigen;

HuskySystem::HuskySystem():Dynamic_System<3,2>(false,false,true,false,false,true,true,true,true,false,Matrix<double,3,3>::Identity(), Dormand_Prince_8th, Dormand_Prince_8th)
{
    return;
}

//done
Matrix<double,3,1> HuskySystem::g(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& ) const
{
    ///<- returns the output y = g(t,x,u)  Multiply by C for output selection

    return x;
}

//done
Matrix<double,3,1> HuskySystem::ode(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& u) const
{
    //<- returns the time derivative of x, f(t,x,u)

    Matrix<double,3,1> x_dot;

    x_dot << u(0)*std::cos(x(2)),u(0)*std::sin(x(2)),u(1);

    return x_dot;
}

//done
Matrix<double,3,3> HuskySystem::g_x_sens(double, const Matrix<double,3,1>& , const Matrix<double,2,1>& ) const//didnt initialize x and u because of warnings
{
    //<- returns dg/dx

    return Matrix<double,3,3>::Identity();
}


Matrix<double,3,3> HuskySystem::ode_x_sens(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& u) const
{
    //<- returns df/dx

    Matrix<double,3,3> df_dx;
    df_dx.setZero(3,3);
    df_dx(0,2) = -u(0)*std::sin(x(2));
    df_dx(1,2) =  u(0)*std::cos(x(2));



    return df_dx;
}
Matrix<double,3,2> HuskySystem::ode_u_sens(double, const Matrix<double,3,1>& x, const Matrix<double,2,1>& u) const
{
    //<- returns df/du
    Matrix<double,3,2> df_du = Matrix<double,3,2>::Zero();


    df_du(0,0) = std::cos(x(2));
    df_du(1,0) = std::sin(x(2));
    df_du(2,1) = 1;

    return df_du;
}


R3Tensor<3,3,3> HuskySystem::ode_x2_sens(double , const Matrix<double,3,1>& x, const Matrix<double,2,1>& u) const
{
    //<- returns d^2f/dx^2
    R3Tensor<3,3,3> df_dx2;
    df_dx2.setZero();
    df_dx2(0,2,2) = -u(0)*std::cos(x(2));
    df_dx2(1,2,2) = -u(0)*std::sin(x(2));

    return df_dx2;

}

R3Tensor<3,3,2> HuskySystem::ode_xu_sens(double , const Matrix<double,3,1>& x, const Matrix<double,2,1>& ) const
{
    //<- returns d^2f/dxdu
    R3Tensor<3,3,2> df_dxu;
    df_dxu.setZero();

    df_dxu(0,2,0) = -std::sin(x(2));
    df_dxu(1,2,0) =  std::cos(x(2));

    return df_dxu;
}


Dynamic_System<3,2>* HuskySystem::clone() const
{
    return new HuskySystem();
}



Matrix<double,3,1> HuskySystem::x( double t) const
{
    //<- returns the state at time t
    Matrix<double,3,1> state;

    double dT = t-to;
    double dT2 = dT*dT;//std::pow(t-to,2);
    double dT3 = dT2*dT;//std::pow(t-to,3);
    double dT4 = dT3*dT;//std::pow(dT,4);
    double dT5 = dT4*dT;

    double cTh_f = std::cos(Xo(2)+U(1)*dT);
    double sTh_f = std::sin(Xo(2)+U(1)*dT);
    double sTh_o = std::sin(Xo(2));
    double cTh_o = std::cos(Xo(2));
    double u1_sq = std::pow(U(1),2);
    double u1_cub = u1_sq*U(1);

    if( u1_sq > 8e-2 )
    {
        state(0) = Xo(0) + U(0)*(sTh_f-sTh_o)/U(1);
        state(1) = Xo(1) - U(0)*(cTh_f-cTh_o)/U(1) ;
        state(2) = Xo(2) + U(1)*dT;


    }
    else
    {
        state(0) = Xo(0) + U(0)*cTh_o*dT;
        state(1) = Xo(1) + U(0)*sTh_o*dT;
        state(2) = Xo(2) + dT*U(1);


       // cout << "cost taylor Mag: " << cTh_o*dT5*u1_cub/120.0  <<" " << sTh_o*dT5*u1_cub/120.0  << endl;

        state(0) += U(0)*U(1)*(-sTh_o*dT2/2.0  - U(1)*cTh_o*dT3/6.0 + u1_sq*sTh_o*dT4/24.0 + cTh_o*dT5*u1_cub/120.0 - sTh_o*dT*dT5*U(1)*u1_cub/720.0);
        state(1) += U(0)*U(1)*( cTh_o*dT2/2.0  - U(1)*sTh_o*dT3/6.0 - u1_sq*cTh_o*dT4/24.0 + sTh_o*dT5*u1_cub/120.0 + cTh_o*dT*dT5*U(1)*u1_cub/720.0);

    }

    return state;
}
Matrix<double,3,1> HuskySystem::x( ) const
{
    //<- returns the state at time t_f
    return x(tf);
}


bool HuskySystem::solve(T_Span t_span, bool, double, int)
{
    //<- Solves the ODE
    to = t_span.t0;
    tf = t_span.tf;
    solved = true;
    return true;
}

Matrix<double,3,3> HuskySystem::Jx_sens() const
{
    //<- returns final value of dXf/dXo
    return(Jx_sens(tf));
}

Matrix<double,3,3> HuskySystem::Jx_sens(double t) const
{
    //<- returns intermediate value of dXf/dXo
    Matrix<double,3,3> Jx = Matrix<double,3,3>::Identity();
    double dT = t-to;
    double dT2 = std::pow(t-to,2);
    double dT3 = dT2*dT;//std::pow(t-to,3);
    double dT4 = dT3*dT;
    double dT5 = dT4*dT;

    double cTh_f = std::cos(Xo(2)+U(1)*dT);
    double sTh_f = std::sin(Xo(2)+U(1)*dT);
    double sTh_o = std::sin(Xo(2));
    double cTh_o = std::cos(Xo(2));
    double u1_sq = std::pow(U(1),2);
    double u1_cub = u1_sq*U(1);

    if( u1_sq > 8e-2 )
    {
        Jx(0,2) = U(0)*(cTh_f-cTh_o)/U(1);
        Jx(1,2) = U(0)*(sTh_f-sTh_o)/U(1);
    }
    else//linearized by taylor theorem about 0
    {
        Jx(0,2) = U(0)*(-sTh_o*dT - U(1)*cTh_o*dT2/2.0 + u1_sq*sTh_o*dT3/6.0 + cTh_o*dT4*u1_cub/24.0 - sTh_o*dT5*U(1)*u1_cub/120.0);
        Jx(1,2) = U(0)*( cTh_o*dT - U(1)*sTh_o*dT2/2.0 - u1_sq*cTh_o*dT3/6.0 + sTh_o*dT4*u1_cub/24.0 + cTh_o*dT5*U(1)*u1_cub/120.0);
    }
    return Jx;
}

Matrix<double,3,2> HuskySystem::Ju_sens() const
{
    //<- returns final value of dXf/dU
    return Ju_sens(tf);
}

Matrix<double,3,2> HuskySystem::Ju_sens(double t) const ///<- returns intermediate value of dXf/dU
{
    //<- returns intermediate value of dXf/dXo
    Matrix<double,3,2> Ju = Matrix<double,3,2>::Zero();
    double dT = t-to;
    double dT2 = dT*dT;//std::pow(dT,2);
    double dT3 = dT2*dT;//std::pow(dT,3);
    double dT4 = dT3*dT;//std::pow(dT,4);
    double dT5 = dT4*dT;
    double dT6 = dT5*dT;

    double cTh_f = std::cos(Xo(2)+U(1)*dT);
    double sTh_f = std::sin(Xo(2)+U(1)*dT);
    double sTh_o = std::sin(Xo(2));
    double cTh_o = std::cos(Xo(2));
    double u1_sq = std::pow(U(1),2);
    double u1_cub = u1_sq*U(1);


   // state(0) = Xo(0) + U(0)*(sTh_f-sTh_o)/U(1);
   // state(1) = Xo(1) - U(0)*(cTh_f-cTh_o)/U(1) ;

    if( u1_sq > 8e-2 )
    {
        Ju(0,0) = (sTh_f-sTh_o)/U(1);
        Ju(1,0) = -(cTh_f-cTh_o)/U(1);

        Ju(0,1) = U(0)*(cTh_f*dT/U(1) - (sTh_f-sTh_o)/u1_sq);
        Ju(1,1) = U(0)*(sTh_f*dT/U(1) + (cTh_f-cTh_o)/u1_sq);

        Ju(2,1) = dT;

    }
    else//linearized by taylor theorem about 0
    {
        Ju(0,0) = cTh_o*dT - U(1)*sTh_o*dT2/2.0 - (u1_sq*cTh_o*dT3/6.0) + sTh_o*dT4*u1_cub/24.0 + cTh_o*dT5*u1_cub*U(1)/120.0 - sTh_o*dT6*u1_cub*u1_sq/720.0;
        Ju(1,0) = sTh_o*dT + U(1)*cTh_o*dT2/2.0 - (u1_sq*sTh_o*dT3/6.0) - cTh_o*dT4*u1_cub/24.0 + sTh_o*dT5*u1_cub*U(1)/120.0 + cTh_o*dT6*u1_cub*u1_sq/720.0;

        //cout << "taylor Mag: " << u1_cub*cTh_o*dT5/30.0  <<" " << u1_cub*sTh_o*dT5/30.0  << endl;
        Ju(0,1) = U(0)*(-sTh_o*dT2/2.0 - U(1)*cTh_o*dT3/3.0 + u1_sq*sTh_o*dT4/8.0 + u1_cub*cTh_o*dT5/30.0 - u1_cub*U(1)*dT6*sTh_o/144.0);
        Ju(1,1) = U(0)*( cTh_o*dT2/2.0 - U(1)*sTh_o*dT3/3.0 - u1_sq*cTh_o*dT4/8.0 + u1_cub*sTh_o*dT5/30.0 + u1_cub*U(1)*dT6*cTh_o/144.0);

        Ju(2,1) = dT;
    }
    return Ju;
}

R3Tensor<3,3,3> HuskySystem::Jxx_sens() const
{
    //<- returns final value of d^2Xf/dXo^2
    return Jxx_sens(tf);
}

R3Tensor<3,3,3> HuskySystem::Jxx_sens(double t) const
{
    //<- returns intermediate value of d^2Xf/dXo^2
    //R3Tensor<3,3,3> Jxx = R3Tensor<3,3,3>::Zero();
    R3Tensor<3,3,3> Jxx;
    Jxx.setZero();

    double dT = t-to;
    double dT2 = dT*dT;//std::pow(dT,2);
    double dT3 = dT2*dT;//std::pow(dT,3);
    double dT4 = dT3*dT;//std::pow(dT,4);
    double dT5 = dT4*dT;
    double dT6 = dT5*dT;

    double cTh_f = std::cos(Xo(2)+U(1)*dT);
    double sTh_f = std::sin(Xo(2)+U(1)*dT);
    double sTh_o = std::sin(Xo(2));
    double cTh_o = std::cos(Xo(2));
    double u1_sq = std::pow(U(1),2);
    double u1_cub = u1_sq*U(1);


    if( u1_sq > 8e-2 )
    {
        Jxx(0,2,2) = -U(0)*(sTh_f-sTh_o)/U(1);
        Jxx(1,2,2) =  U(0)*(cTh_f-cTh_o)/U(1);
    }
    else//jxxu*u1
    {
        Jxx(0,2,2) = U(0)*(-cTh_o*dT + U(1)*sTh_o*dT2/2.0 + cTh_o*dT3*u1_sq/6.0 - sTh_o*dT4*u1_cub/24.0 - cTh_o*dT5*u1_sq*u1_sq/120.0 + sTh_o*dT6*u1_sq*u1_cub/720.0);
        Jxx(1,2,2) = U(0)*(-sTh_o*dT - U(1)*cTh_o*dT2/2.0 + sTh_o*dT3*u1_sq/6.0 + cTh_o*dT4*u1_cub/24.0 - sTh_o*dT5*u1_sq*u1_sq/120.0 - cTh_o*dT6*u1_sq*u1_cub/720.0);
    }
    return Jxx;
}

R3Tensor<3,2,2> HuskySystem::Juu_sens() const
{
    //<- returns final value of d^2Xf/dU^2
    return Juu_sens(tf);
}

R3Tensor<3,2,2> HuskySystem::Juu_sens(double t) const
{
    //<- returns intermediate value of d^2Xf/dU^2
    //R3Tensor<3,2,2> Juu = R3Tensor<3,2,2>::Zero();
    R3Tensor<3,2,2> Juu;
    Juu.setZero();

    double dT = t-to;
    double dT2 = dT*dT;
    double dT3 = dT2*dT;
    double dT4 = dT3*dT;
    double dT5 = dT4*dT;
    double dT6 = dT5*dT;
    double cTh_f = std::cos(Xo(2)+U(1)*dT);
    double sTh_f = std::sin(Xo(2)+U(1)*dT);
    double sTh_o = std::sin(Xo(2));
    double cTh_o = std::cos(Xo(2));
    double u1_sq = std::pow(U(1),2);
    double u1_cu = std::pow(U(1),3);


    if( u1_sq > 8e-2 )
    {
        Juu(0,1,0) = cTh_f*dT/U(1)-(sTh_f-sTh_o)/u1_sq;
        Juu(0,0,1) = Juu(0,1,0);

        Juu(1,0,1) = sTh_f*dT/U(1) + (cTh_f-cTh_o)/u1_sq;
        Juu(1,1,0) = Juu(1,0,1);

        Juu(1,1,1) = U(0)*( cTh_f*dT2/U(1) - 2.0*sTh_f*dT/u1_sq - 2.0*(cTh_f-cTh_o)/u1_cu);
        Juu(0,1,1) = U(0)*(-sTh_f*dT2/U(1) - 2.0*cTh_f*dT/u1_sq + 2.0*(sTh_f-sTh_o)/u1_cu);

    }
    else//juuu*u1
    {
        Juu(0,1,0) = -sTh_o*dT2/2.0 - U(1)*cTh_o*dT3/3.0  + u1_sq*sTh_o*dT4/8.0 + u1_cu*cTh_o*dT5/30.0 - u1_cu*U(1)*sTh_o*dT6/144.0;
        Juu(0,0,1) = Juu(0,1,0);

        Juu(1,1,0) = cTh_o*dT2/2.0 - U(1)*sTh_o*dT3/3.0 - u1_sq*cTh_o*dT4/8.0 + u1_cu*sTh_o*dT5/30.0 + u1_cu*U(1)*cTh_o*dT6/144.0;
        Juu(1,0,1) = Juu(1,1,0);

        Juu(0,1,1) = U(0)*(-cTh_o*dT3/3.0 + U(1)*sTh_o*dT4/4.0 + u1_sq*cTh_o*dT5/10.0 - u1_cu*sTh_o*dT6/36.0);
        Juu(1,1,1) = U(0)*(-sTh_o*dT3/3.0 - U(1)*cTh_o*dT4/4.0 + u1_sq*sTh_o*dT5/10.0 + u1_cu*cTh_o*dT6/36.0);

        //cout << " below tol " << endl;
    }
    return Juu;
}

R3Tensor<3,3,2> HuskySystem::Jxu_sens() const
{
    //<- returns final value of d^2Xf/dU^2
    return Jxu_sens(tf);
}

R3Tensor<3,3,2> HuskySystem::Jxu_sens(double t) const
{
    //<- returns intermediate value of d^2Xf/dU^2
    //R3Tensor<3,3,2> Jxu = R3Tensor<3,3,2>::Zero();
    R3Tensor<3,3,2> Jxu;
    Jxu.setZero();

    double dT = t-to;
    double dT2 = dT*dT;
    double dT3 = dT2*dT;
    double dT4 = dT3*dT;
    double dT5 = dT4*dT;
    double dT6 = dT5*dT;
    double cTh_f = std::cos(Xo(2)+U(1)*dT);
    double sTh_f = std::sin(Xo(2)+U(1)*dT);
    double sTh_o = std::sin(Xo(2));
    double cTh_o = std::cos(Xo(2));
    double u1_sq = std::pow(U(1),2);
    double u1_cu = std::pow(U(1),3);


    if( u1_sq > 8e-2 )
    {
        Jxu(0,2,0) = (cTh_f-cTh_o)/U(1);
        Jxu(1,2,0) = (sTh_f-sTh_o)/U(1);

        Jxu(0,2,1) = -U(0)*(sTh_f*dT/U(1) -(cTh_f-cTh_o)/u1_sq);
        Jxu(1,2,1) =  U(0)*(cTh_f*dT/U(1) -(sTh_f-sTh_o)/u1_sq);
    }
    else//jxuu*u1
    {
        Jxu(0,2,0) = -sTh_o*dT - U(1)*cTh_o*dT2/2.0 + u1_sq*sTh_o*dT3/6.0 + u1_cu*cTh_o*dT4/24.0 - u1_cu*U(1)*sTh_o*dT5/120.0 - u1_cu*u1_sq*cTh_o*dT6/720.0;
        Jxu(1,2,0) =  cTh_o*dT - U(1)*sTh_o*dT2/2.0 - u1_sq*cTh_o*dT3/6.0 + u1_cu*sTh_o*dT4/24.0 + u1_cu*U(1)*cTh_o*dT5/120.0 - u1_cu*u1_sq*sTh_o*dT6/720.0;

        Jxu(0,2,1) = U(0)*(-cTh_o*dT2/2.0 + U(1)*sTh_o*dT3/3.0 + u1_sq*cTh_o*dT4/8.0 - u1_cu*sTh_o*dT5/30.0 - u1_cu*U(1)*cTh_o*dT6/144.0);
        Jxu(1,2,1) = U(0)*(-sTh_o*dT2/2.0 - U(1)*cTh_o*dT3/3.0 + u1_sq*sTh_o*dT4/8.0 + u1_cu*cTh_o*dT5/30.0 - u1_cu*U(1)*sTh_o*dT6/144.0);
    }

    return Jxu;
}
