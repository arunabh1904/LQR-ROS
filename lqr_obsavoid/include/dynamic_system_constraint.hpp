#ifndef DYNAMIC_SYSTEM_CONSTRAINT_HPP
#define DYNAMIC_SYSTEM_CONSTRAINT_HPP

#include "dynamic_system_constraint.h"
#include <iostream>
using std::cout;
using std::endl;

template<int M, int N>
Dynamic_System_Constraint<M,N>::Dynamic_System_Constraint(bool has_x_sens_, bool has_xx_sens_, bool has_ux_sens_, bool has_u_sens_, bool has_uu_sens_):
    b_u_sens(has_u_sens_),
    b_x_sens(has_x_sens_),
    b_uu_sens(has_uu_sens_),
    b_xx_sens(has_xx_sens_),
    b_ux_sens(has_ux_sens_)
{
    ;
}

//template<int M, int N>
//double Dynamic_System_Constraint<M,N>::constraint_cost(double t, const Eigen::Matrix<double,M,1>& x, Eigen::Matrix<double,N,1> u) const
//{
//    //<- The hard constraint. Should be zero when satisfied
//    assert(false && "Dynamic_System_Constraint<M,N>::constraint_cost must be implemented by an inheriting function!");
//    return 0;
//}

template<int M, int N>
bool Dynamic_System_Constraint<M,N>::has_x_sens() const
{
    //<- Returns if the constraint has a non-zero d/dx sensitivity
    return b_x_sens;
}

template<int M, int N>
bool Dynamic_System_Constraint<M,N>::has_xx_sens() const
{
    //<- Returns if the constraint has a non-zero d^2/dxdx sensitivity
    return b_xx_sens;
}

template<int M, int N>
bool Dynamic_System_Constraint<M,N>::has_u_sens() const
{
    //<- Returns if the constraint has a non-zero d/du sensitivity
    return b_u_sens;
}

template<int M, int N>
bool Dynamic_System_Constraint<M,N>::has_uu_sens() const
{
    //<- Returns if the constraint has a non-zero d^2/dudu sensitivity
    return b_uu_sens;
}

template<int M, int N>
bool Dynamic_System_Constraint<M,N>::has_ux_sens() const
{
    //<- Returns if the constraint has a non-zero d^2/dxdu sensitivity
    return b_ux_sens;
}

template<int M, int N>
Eigen::Matrix<double,M,1> Dynamic_System_Constraint<M,N>::x_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>&, Dynamic_System_Constraint_Type) const
{
    //<- the derivative of the constraint with respect to state
    assert( !b_x_sens && "Dynamic_System_Constraint calling ZERO x_sens when inheriting class indicates it should be nonzero" );
    return Eigen::Matrix<double,M,1>::Zero();
}

template<int M, int N>
Eigen::Matrix<double,N,1> Dynamic_System_Constraint<M,N>::u_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>&, Dynamic_System_Constraint_Type) const
{
    //<- the derivative of the constraint with respect to control
    assert( !b_u_sens && "Dynamic_System_Constraint calling ZERO u_sens when inheriting class indicates it should be nonzero" );
    return Eigen::Matrix<double,N,1>::Zero();
}

template<int M, int N>
Eigen::Matrix<double,M,M> Dynamic_System_Constraint<M,N>::xx_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>&, Dynamic_System_Constraint_Type) const
{
    //<- the second derivative of the constraint with respect to state^2
    assert( !b_xx_sens && "Dynamic_System_Constraint calling ZERO xx_sens when inheriting class indicates it should be nonzero" );
    return Eigen::Matrix<double,M,M>::Zero();
}

template<int M, int N>
Eigen::Matrix<double,N,N> Dynamic_System_Constraint<M,N>::uu_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>&, Dynamic_System_Constraint_Type) const
{
    //<- the derivative of the constraint with respect to control^2
    assert( !b_uu_sens && "Dynamic_System_Constraint calling ZERO uu_sens when inheriting class indicates it should be nonzero" );
    return Eigen::Matrix<double,N,N>::Zero();
}

template<int M, int N>
Eigen::Matrix<double,N,M> Dynamic_System_Constraint<M,N>::ux_sens(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>&, Dynamic_System_Constraint_Type) const
{
    //<- the derivative of the constraint with respect to control then state
    assert( !b_ux_sens && "Dynamic_System_Constraint calling ZERO ux_sens when inheriting class indicates it should be nonzero" );
    return Eigen::Matrix<double,N,M>::Zero();
}

template<int M, int N>
bool Dynamic_System_Constraint<M,N>::verify_constraint_jacobians(double t, const Eigen::Matrix<double,M,1>& x, const Eigen::Matrix<double,N,1>& u, Dynamic_System_Constraint_Type type)
{
    Eigen::Matrix<double,M,1> dx; dx.setZero(M,1);
    Eigen::Matrix<double,N,1> du; du.setZero(N,1);

    double costNom;
    Eigen::Matrix<double,M,1> x_sens_nom, x_sens_num;
    Eigen::Matrix<double,N,1> u_sens_nom, u_sens_num;
    Eigen::Matrix<double,M,M> xx_sens_nom, xx_sens_num;
    Eigen::Matrix<double,N,N> uu_sens_nom, uu_sens_num;
    Eigen::Matrix<double,N,M> ux_sens_nom, ux_sens_num;

    costNom = constraint_cost(t,x,u,type);
    x_sens_nom = x_sens(t,x,u,type);
    u_sens_nom = u_sens(t,x,u,type);
    xx_sens_nom = xx_sens(t,x,u,type);
    uu_sens_nom = uu_sens(t,x,u,type);
    ux_sens_nom = ux_sens(t,x,u,type);

    double delta = 1e-8;
    for(int i=0; i<M; i++)
    {
        dx(i) = delta;
        x_sens_num(i) = (constraint_cost(t,x+dx,u,type)-costNom)/delta;
        xx_sens_num.block(0,i,M,1) = (x_sens(t,x+dx,u,type)-x_sens_nom)/delta;
        ux_sens_num.block(0,i,N,1) = (u_sens(t,x+dx,u,type)-u_sens_nom)/delta;

        dx(i) = 0;
    }

    for(int i=0; i<N; i++)
    {
        du(i) = delta;
        u_sens_num(i) = (constraint_cost(t,x,u+du,type)-costNom)/delta;
        uu_sens_num.block(0,i,N,1) = (u_sens(t,x,u+du,type)-u_sens_nom)/delta;

        du(i) = 0;
    }

    double x_sens_err = (x_sens_nom-x_sens_num).norm();
    if( has_x_sens() && x_sens_err > 1e-5)
        x_sens_err /= x_sens_num.norm();

    double u_sens_err = (u_sens_nom-u_sens_num).norm();
    if( has_u_sens() && u_sens_err > 1e-5)
        u_sens_err /= u_sens_num.norm();


    double xx_sens_err = (xx_sens_nom-xx_sens_num).norm();
    if( has_xx_sens()  && xx_sens_err > 1e-5)
        xx_sens_err /= xx_sens_num.norm();

    double uu_sens_err = (uu_sens_nom-uu_sens_num).norm();
    if( has_uu_sens() && uu_sens_err > 1e-5)
        uu_sens_err /= uu_sens_num.norm();

    double ux_sens_err = (ux_sens_nom-ux_sens_num).norm();
    if( has_ux_sens() && ux_sens_err > 1e-5)
        ux_sens_err /= ux_sens_num.norm();


    double tol = 1e-7;
    bool retVal = true;
    if( x_sens_err > tol )
    {
        retVal = false;
        cout << "x_sens_err: " << x_sens_err << endl;
        cout << "x_sens analitical: " << x_sens_nom.transpose() << endl;
        cout << "x_sens numeric:    " << x_sens_num.transpose() << endl;
        cout << "difference:        " << (x_sens_nom - x_sens_num).transpose() << endl;
    }

    if( u_sens_err > tol )
    {
        retVal = false;
        cout << "u_sens_err: " << u_sens_err << endl;
        cout << "u_sens analitical: " << u_sens_nom.transpose() << endl;
        cout << "u_sens numeric:    " << u_sens_num.transpose() << endl;
        cout << "difference:        " << (u_sens_nom - u_sens_num).transpose() << endl;
    }

    if( xx_sens_err > tol )
    {
        retVal = false;
        cout << "xx_sens_err: " << xx_sens_err << endl;
        cout << "xx_sens analitical:\n" << xx_sens_nom << endl<<endl;
        cout << "xx_sens numeric:\n" << xx_sens_num << endl <<endl;
        cout << "difference:\n" << (xx_sens_nom - xx_sens_num) << endl<<endl;
    }

    if( uu_sens_err > tol )
    {
        retVal = false;
        cout << "uu_sens_err: " << uu_sens_err << endl;
        cout << "uu_sens analitical:\n" << uu_sens_nom << endl<<endl;
        cout << "uu_sens numeric:\n" << uu_sens_num << endl <<endl;
        cout << "difference:\n" << (uu_sens_nom - uu_sens_num) << endl<<endl;
    }

    if( ux_sens_err > tol )
    {
        retVal = false;
        cout << "ux_sens_err: " << ux_sens_err << endl;
        cout << "ux_sens analitical:\n" << ux_sens_nom << endl<<endl;
        cout << "ux_sens numeric:\n" << ux_sens_num << endl <<endl;
        cout << "difference:\n" << (ux_sens_nom - ux_sens_num) << endl<<endl;
    }

    if(retVal)
        cout<< "All derivatives check out: % errors: "
            << x_sens_err*100 << ", "
            << u_sens_err*100 << ", "
            << xx_sens_err*100 << ", "
            << uu_sens_err*100 << ", "
            << ux_sens_err*100 << "%" << endl;


    return retVal;

}

template<int M, int N>
double Dynamic_System_Constraint<M,N>::logistic(double x, double k) const
{
    double kx = k*x;

    if( kx > 100 )
        return x;
    if( kx < -100 )
        return 0;

    double mekx = std::exp(-kx);
    return x/(1.0+mekx);
}

template<int M, int N>
double Dynamic_System_Constraint<M,N>::logisticDer(double x, double k) const
{
    double kx = k*x;

    if( kx > 100 )
        return 1;
    if( kx < -100 )
        return 0;

    double ekx = std::exp(kx);


    return ekx*(k*x+ekx+1.0)/std::pow(ekx+1.0,2);
}

template<int M, int N>
double Dynamic_System_Constraint<M,N>::logisitcDer2(double x, double k) const
{
    double kx = k*x;

    if( kx > 100 )
        return 0;
    if( kx < -100 )
        return 0;

    double ekx = std::exp(kx);

    return -k*ekx*(-kx+ekx*(kx-2.0)-2.0)/std::pow(ekx+1.0,3);
}

template<int M, int N>
Eigen::Matrix<double,N,1> Dynamic_System_Constraint<M,N>::inforce_u_constraint(double, const Eigen::Matrix<double,M,1>&, const Eigen::Matrix<double,N,1>&u)const
{
    assert( !has_u_sens() && "U constriant not being inforced" );
    return u;
}


#endif // DYNAMIC_SYSTEM_CONSTRAINT_HPP

