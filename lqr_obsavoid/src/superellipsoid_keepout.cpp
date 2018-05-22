#ifndef SUPERELLIPSOID_KEEPOUT_HPP
#define SUPERELLIPSOID_KEEPOUT_HPP

#include "superellipsoid_keepout.h"

SuperEllipsoidKeepout::SuperEllipsoidKeepout( double xLength, double yLength, double zLength, double r_, double t_, TYPE constraint_type_)
{

    p_field_mag = 0;
    n_pow = 0;

    xLen = xLength;
    yLen = yLength;
    zLen = zLength;
    r = r_;
    t = t_;

    constraint_type = constraint_type_;

    logistic_K = 100.0;
    max_threshold = 0.5*this->logisitcDer2_kx_zero/logistic_K;
    min_threshold = 0.5*this->logisitcDer2_kx_zero/(-logistic_K);
}

void SuperEllipsoidKeepout::setMagnitude( double newMag )
{
    p_field_mag = newMag;

}

void SuperEllipsoidKeepout::setExponent( double newExponent )
{
    n_pow = newExponent;
}

void SuperEllipsoidKeepout::setSuperEllipseT(double newT )
{

    t = newT;
}

void SuperEllipsoidKeepout::setSuperEllipseR(double newR )
{
    r = newR;
}

SuperEllipsoidKeepout::SuperEllipsoidKeepout( double p_field_mag_, double p_field_pow_, double xLength, double yLength, double zLength, double r_, double t_, TYPE constraint_type_ )
{
    p_field_mag = p_field_mag_;
    n_pow = p_field_pow_;

    xLen = xLength;
    yLen = yLength;
    zLen = zLength;
    r = r_;
    t = t_;

    constraint_type = constraint_type_;

    logistic_K = 100.0;
    max_threshold = 0.5*this->logisitcDer2_kx_zero/logistic_K;
    min_threshold = 0.5*this->logisitcDer2_kx_zero/(-logistic_K);
}


void SuperEllipsoidKeepout::setPointList(const std::vector<Eigen::Vector3d>& pointList)
{
    point_list = pointList;
}


double SuperEllipsoidKeepout::cost(const Eigen::Vector3d& robot_pos_inPointFrame, const Eigen::Matrix3d& rotation_robot_intoPointFrame) const
{
    //<- The hard constraint. Should be zero when satisfied

    Eigen::Vector3d disp_in_robotFrame;
    Eigen::Matrix3d Rinv = rotation_robot_intoPointFrame.transpose();


    double ret_cost = 0;

    for(unsigned int i=0; i<point_list.size(); i++ )
    {
        disp_in_robotFrame = Rinv*(robot_pos_inPointFrame - point_list[i]);
        double xDif = std::abs(disp_in_robotFrame.x())*2.0/xLen;
        double yDif = std::abs(disp_in_robotFrame.y())*2.0/yLen;
        double zDif = std::abs(disp_in_robotFrame.z())*2.0/zLen;

        double cost_tmp = calc_cost(xDif, yDif, zDif); //std::pow( std::pow(2.0*(x-x_c)/xLen,r) + std::pow(2.0*(y-y_c)/yLen,r), t/r) + std::pow(2.0*(z-z_c)/zLen,t) - 1.0;

        switch(constraint_type)
        {
        case KEEP_OUT:
            if( cost_tmp > min_threshold )
                cost_tmp = 0;
            else
                cost_tmp = logistic(cost_tmp,-logistic_K);
            break;
        case KEEP_IN:
            if( cost_tmp < max_threshold )
                cost_tmp = 0;
            else
                cost_tmp = logistic(cost_tmp,logistic_K);
            break;
        case KEEP_OUT_POTENTIAL_FIELD:
            if( cost_tmp < 0 )
                cost_tmp = 0;
            else
                cost_tmp = p_field_mag/std::pow(cost_tmp,n_pow);
            break;

        default:
            break;
        }

        ret_cost += cost_tmp;
    }

    return ret_cost;
}

Eigen::Matrix<double,6,1> SuperEllipsoidKeepout::cost_jacobian( const Eigen::Vector3d& robot_pos_inPointFrame, const Eigen::Matrix3d& R) const
{
    //<- the derivative of the constraint with respect to state
    Eigen::Vector3d disp_in_robotFrame, disp_in_pointFrame;
    Eigen::Matrix3d Rinv = R.transpose();

    Eigen::Matrix<double,6,1> dcdx(Eigen::Matrix<double,6,1>::Zero());

    for(unsigned int i=0; i<point_list.size(); i++ )
    {
        disp_in_pointFrame = (robot_pos_inPointFrame - point_list[i]);
        disp_in_robotFrame = Rinv*disp_in_pointFrame;


        double xDif = disp_in_robotFrame.x()*2.0/xLen;
        double yDif = disp_in_robotFrame.y()*2.0/yLen;
        double zDif = disp_in_robotFrame.z()*2.0/zLen;
        int xDif_s = sgn(xDif);
        int yDif_s = sgn(yDif);
        int zDif_s = sgn(zDif);
        xDif = std::abs(xDif);
        yDif = std::abs(yDif);
        zDif = std::abs(zDif);

        double p_cost = calc_cost(xDif, yDif, zDif);
        double logisticDerVal = 1;
        switch(constraint_type)
        {
        case KEEP_OUT:
            if( p_cost > min_threshold )
                continue;
            else
                logisticDerVal = logisticDer(p_cost,-logistic_K);
            break;
        case KEEP_IN:
            if( p_cost < max_threshold )
                continue;
            else
                logisticDerVal = logisticDer(p_cost,logistic_K);
            break;
        case KEEP_OUT_POTENTIAL_FIELD:
            if( p_cost < 0 )
                continue;
            else
                logisticDerVal = -n_pow*p_field_mag/std::pow(p_cost, n_pow+1.0);
            break;

        default:
            break;
        }

        Eigen::Vector3d dfdp = Eigen::Vector3d::Zero();
        using std::pow;
        double xr = pow(xDif,r);
        double yr = pow(yDif,r);
        double xr1 = pow(xDif,r-1);
        double yr1 = pow(yDif,r-1);
        dfdp(0) = t*xr1*pow(xr+yr,t/r-1);
        dfdp(1) = t*yr1*pow(xr+yr,t/r-1);
        dfdp(2) = t*pow(zDif,t-1);

        Eigen::Matrix3d sgnDif;
        sgnDif<< 2.0/xLen*xDif_s,0,0,0,2.0/yLen*yDif_s,0,0,0,2.0/zLen*zDif_s; //diagonal of sgn_dif


        dcdx.block(0,0,3,1) += logisticDerVal*R*sgnDif*dfdp;
        dcdx.block(3,0,3,1) += logisticDerVal*-cross(disp_in_pointFrame)*R*sgnDif*dfdp;


    }

    return dcdx;
}

Eigen::Matrix<double,6,6> SuperEllipsoidKeepout::cost_hessian( const Eigen::Vector3d& robot_pos_inPointFrame, const Eigen::Matrix3d& R) const
{
    //<- the second derivative of the constraint with respect to state^2
    Eigen::Vector3d disp_in_robotFrame,disp_in_pointFrame;
    Eigen::Matrix3d Rinv = R.transpose();

    Eigen::Matrix<double,6,6> dcdxx(Eigen::Matrix<double,6,6>::Zero());
    Eigen::Matrix3d dcdxw_tmp(Eigen::Matrix3d::Zero());

    for(unsigned int i=0; i<point_list.size(); i++ )
    {
        disp_in_pointFrame = (robot_pos_inPointFrame - point_list[i]);
        disp_in_robotFrame = Rinv*disp_in_pointFrame;

        double xDif = disp_in_robotFrame.x()*2.0/xLen;
        double yDif = disp_in_robotFrame.y()*2.0/yLen;
        double zDif = disp_in_robotFrame.z()*2.0/zLen;
        int xDif_s = sgn(xDif);
        int yDif_s = sgn(yDif);
        int zDif_s = sgn(zDif);
        xDif = std::abs(xDif);
        yDif = std::abs(yDif);
        zDif = std::abs(zDif);


        double p_cost = calc_cost(xDif, yDif, zDif);
        double logisticDerVal = 1;
        double logisticDer2Val = 0;

        switch(constraint_type)
        {
        case KEEP_OUT:
            if( p_cost > min_threshold )
                continue;
            else
            {
                logisticDerVal = logisticDer(p_cost,-logistic_K);
                logisticDer2Val = logisitcDer2(p_cost,-logistic_K);
            }
            break;
        case KEEP_IN:
            if( p_cost < max_threshold )
                continue;
            else
            {
                logisticDerVal = logisticDer(p_cost,logistic_K);
                logisticDer2Val = logisitcDer2(p_cost,logistic_K);
            }
            break;

        case KEEP_OUT_POTENTIAL_FIELD:
            if( p_cost < 0 )
                continue;
            else
            {
                logisticDerVal = -n_pow*p_field_mag/std::pow(p_cost, n_pow+1.0);
                logisticDer2Val = p_field_mag*n_pow*(n_pow+1.0)/std::pow(p_cost, n_pow+2.0);
            }
            break;

        default:
            break;
        }


        Eigen::Matrix3d dfdp2 = Eigen::Matrix3d::Zero();
        Eigen::Vector3d dfdp = Eigen::Vector3d::Zero();
        /*
            dxx = t*x^(r - 2)*(x^r + y^r)^(t/r - 2)*(r*y^r + t*x^r - x^r - y^r)
            dyy = t*y^(r - 2)*(x^r + y^r)^(t/r - 2)*(r*x^r + t*y^r - x^r - y^r)
            dxy = r*t*x^(r - 1)*y^(r - 1)*(t/r - 1)*(x^r + y^r)^(t/r - 2)
            dzz = t*z^(t - 2)*(t - 1)
        */
        using std::pow;
        double xr = pow(xDif,r);
        double yr = pow(yDif,r);
        double xr1 = pow(xDif,r-1);
        double yr1 = pow(yDif,r-1);
        double xryr_tr2 = pow(xr+yr,t/r-2);
        dfdp2(0,0) = t*pow(xDif,r-2)*xryr_tr2*(r*yr+t*xr-xr-yr);
        dfdp2(1,1) = t*pow(yDif,r-2)*xryr_tr2*(r*xr+t*yr-xr-yr);
        dfdp2(0,1) = r*t*xr1*yr1*(t/r-1)*xryr_tr2;
        dfdp2(1,0) = dfdp2(0,1);
        dfdp2(2,2) = t*(t-1)*pow(zDif,t-2);

        /*
            dx = t*x^(r - 1)*(x^r + y^r)^(t/r - 1)
            dy = t*y^(r - 1)*(x^r + y^r)^(t/r - 1)
            dz = t*z^(t - 1)

            f(abs(R'*(x-p)))
            df/dxx*sgn(R'*(x-p))*cross(x-p)
            dfdw = -dfdx*sgn()*lambda*cross(x-p)*R'
                  R*cross(x-p)*lambda*sgn()*dfdx
        */
        dfdp(0) = t*xr1*pow(xr+yr,t/r-1);
        dfdp(1) = t*yr1*pow(xr+yr,t/r-1);
        dfdp(2) = t*pow(zDif,t-1);

        Eigen::Matrix3d sgnDif;
        sgnDif<< 2.0/xLen*xDif_s,0,0,0,2.0/yLen*yDif_s,0,0,0,2.0/zLen*zDif_s; //diagonal of sgn_dif

        dcdxw_tmp = cross(R*sgnDif*dfdp) -cross(disp_in_pointFrame)*R*sgnDif*dfdp2*sgnDif*Rinv;

        Eigen::Matrix<double,6,1> dcdx;
        dcdx.segment(0,3) = R*sgnDif*dfdp;
        dcdx.segment(3,3) = -cross(disp_in_pointFrame)*R*sgnDif*dfdp;

        dcdxx.block(0,0,3,3) += logisticDerVal*R*sgnDif*dfdp2*sgnDif*Rinv;
        dcdxx.block(3,3,3,3) += logisticDerVal*(cross(disp_in_pointFrame)*R*cross(sgnDif*dfdp)*Rinv - cross(disp_in_pointFrame)*R*sgnDif*dfdp2*sgnDif*Rinv*cross(disp_in_pointFrame));
        dcdxx.block(3,0,3,3) += logisticDerVal*dcdxw_tmp;
        dcdxx.block(0,3,3,3) += logisticDerVal*dcdxw_tmp.transpose();

        dcdxx += logisticDer2Val*(dcdx*dcdx.transpose());

    }

    return dcdxx;

}



double SuperEllipsoidKeepout::calc_cost(double xdif, double ydif, double zdif) const
{
    return std::pow( std::pow(xdif,r) + std::pow(ydif,r), t/r) + std::pow(zdif,t) - 1.0;
}

int SuperEllipsoidKeepout::sgn(double val)
{
    if( val > 0 )
        return 1;
    else if( val < 0 )
        return -1;

    return 0;
}

Eigen::Matrix3d SuperEllipsoidKeepout::cross(const Eigen::Vector3d& other)
{
    Eigen::Matrix3d retMat(Eigen::Matrix3d::Zero());
    retMat(0,1) = -other.z();
    retMat(0,2) = other.y();
    retMat(1,0) = other.z();
    retMat(1,2) = -other.x();
    retMat(2,0) = -other.y();
    retMat(2,1) = other.x();

    return retMat;
}

#include "Eigen/Geometry"
#include "iostream"
bool SuperEllipsoidKeepout::verify_numerically(bool print_to_screen)
{
    double delta = 1e-8;
    std::vector<Eigen::Vector3d> pointListSave = point_list;
    point_list.clear();

    Eigen::Vector3d pos(.1,-.2,.3);

    point_list.push_back(Eigen::Vector3d(xLen/6.0+pos.x(),yLen/6.0+pos.y(), pos.z()-zLen/7.0));



    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();

    double cost_ = cost(pos,R);
    Eigen::Matrix<double,6,1> analiticJacob, numJacob;
    Eigen::Matrix<double,6,6> analiticHessian, numHessian;
    analiticJacob = cost_jacobian(pos,R);
    analiticHessian = cost_hessian(pos,R);

    numJacob(0) = (cost(pos+Eigen::Vector3d(delta,0,0),R)-cost_)/delta;
    numJacob(1) = (cost(pos+Eigen::Vector3d(0,delta,0),R)-cost_)/delta;
    numJacob(2) = (cost(pos+Eigen::Vector3d(0,0,delta),R)-cost_)/delta;
    numJacob(3) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-cost_)/delta;
    numJacob(4) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-cost_)/delta;
    numJacob(5) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-cost_)/delta;

    numHessian.block(0,0,6,1) = (cost_jacobian(pos+Eigen::Vector3d(delta,0,0),R)-analiticJacob)/delta;
    numHessian.block(0,1,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,delta,0),R)-analiticJacob)/delta;
    numHessian.block(0,2,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,delta),R)-analiticJacob)/delta;
    numHessian.block(0,3,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-analiticJacob)/delta;
    numHessian.block(0,4,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-analiticJacob)/delta;
    numHessian.block(0,5,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-analiticJacob)/delta;

    using namespace std;
    double jError = (analiticJacob - numJacob).norm()/analiticJacob.norm();
    double hError = (analiticHessian-numHessian).norm()/analiticHessian.norm();
    if( cost_ == 0)
    {
        jError = (analiticJacob - numJacob).norm();
        hError = (analiticHessian-numHessian).norm();
    }
    if( print_to_screen)
    {
        cout << "************ superellipsoid verification **************"<< endl;
        cout << "Point Inside Rot 1" << endl;
        cout << "Analitic Jacobian:\t" << analiticJacob.transpose() << endl;
        cout << "Numeric Jacobian:\t" << numJacob.transpose() << endl;
        cout << "Difference:\t" << (analiticJacob - numJacob).transpose() <<endl;
        cout << "Analitic Hessian:\t" << endl << analiticHessian << endl << endl;
        cout << "numeric Hessian:\t" << endl << numHessian << endl <<endl;
        cout << "Difference:\t" << endl << (analiticHessian-numHessian)<<endl<<endl;
        cout << "Analitic Jacobian PE: " << jError*100 << "%\tAnalitic  Hessian PE: " << hError*100 << "%"<<endl<<endl;
    }
    bool retVal = jError < delta*5*6;
    retVal &= hError < delta*5*36;

    R = Eigen::AngleAxisd(.6,Eigen::Vector3d::UnitY())*Eigen::AngleAxisd(-.4,Eigen::Vector3d::UnitX())*Eigen::AngleAxisd(.2,Eigen::Vector3d::UnitZ());

    cost_ = cost(pos,R);
    analiticJacob = cost_jacobian(pos,R);
    analiticHessian = cost_hessian(pos,R);

    numJacob(0) = (cost(pos+Eigen::Vector3d(delta,0,0),R)-cost_)/delta;
    numJacob(1) = (cost(pos+Eigen::Vector3d(0,delta,0),R)-cost_)/delta;
    numJacob(2) = (cost(pos+Eigen::Vector3d(0,0,delta),R)-cost_)/delta;
    numJacob(3) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-cost_)/delta;
    numJacob(4) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-cost_)/delta;
    numJacob(5) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-cost_)/delta;

    numHessian.block(0,0,6,1) = (cost_jacobian(pos+Eigen::Vector3d(delta,0,0),R)-analiticJacob)/delta;
    numHessian.block(0,1,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,delta,0),R)-analiticJacob)/delta;
    numHessian.block(0,2,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,delta),R)-analiticJacob)/delta;
    numHessian.block(0,3,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-analiticJacob)/delta;
    numHessian.block(0,4,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-analiticJacob)/delta;
    numHessian.block(0,5,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-analiticJacob)/delta;

    jError = (analiticJacob - numJacob).norm()/analiticJacob.norm();
    hError = (analiticHessian-numHessian).norm()/analiticHessian.norm();
    if( cost_ == 0)
    {
        jError = (analiticJacob - numJacob).norm();
        hError = (analiticHessian-numHessian).norm();
    }
    if( print_to_screen)
    {
        cout << "Point Inside Rot 2" << endl;
        cout << "Analitic Jacobian:\t" << analiticJacob.transpose() << endl;
        cout << "Numeric Jacobian:\t" << numJacob.transpose() << endl;
        cout << "Difference:\t" << (analiticJacob - numJacob).transpose() <<endl;
        cout << "Analitic Hessian:\t" << endl << analiticHessian << endl << endl;
        cout << "numeric Hessian:\t" << endl << numHessian << endl <<endl;
        cout << "Difference:\t" << endl << (analiticHessian-numHessian)<<endl<<endl;
        cout << "Analitic Jacobian PE: " << jError*100 << "%\tAnalitic  Hessian PE: " << hError*100 << "%"<<endl<<endl;
    }
    retVal &= jError < delta*5*6;
    retVal &= hError < delta*5*36;

    point_list.clear();
    point_list.push_back(Eigen::Vector3d(4*xLen/3.0+pos.x(),3*yLen/2.0+pos.y(), pos.z()-2*zLen/3.0));

    R = Eigen::Matrix3d::Identity();

    cost_ = cost(pos,R);
    analiticJacob = cost_jacobian(pos,R);
    analiticHessian = cost_hessian(pos,R);

    numJacob(0) = (cost(pos+Eigen::Vector3d(delta,0,0),R)-cost_)/delta;
    numJacob(1) = (cost(pos+Eigen::Vector3d(0,delta,0),R)-cost_)/delta;
    numJacob(2) = (cost(pos+Eigen::Vector3d(0,0,delta),R)-cost_)/delta;
    numJacob(3) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-cost_)/delta;
    numJacob(4) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-cost_)/delta;
    numJacob(5) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-cost_)/delta;

    numHessian.block(0,0,6,1) = (cost_jacobian(pos+Eigen::Vector3d(delta,0,0),R)-analiticJacob)/delta;
    numHessian.block(0,1,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,delta,0),R)-analiticJacob)/delta;
    numHessian.block(0,2,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,delta),R)-analiticJacob)/delta;
    numHessian.block(0,3,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-analiticJacob)/delta;
    numHessian.block(0,4,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-analiticJacob)/delta;
    numHessian.block(0,5,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-analiticJacob)/delta;

    jError = (analiticJacob - numJacob).norm()/analiticJacob.norm();
    hError = (analiticHessian-numHessian).norm()/analiticHessian.norm();
    if( cost_ == 0)
    {
        jError = (analiticJacob - numJacob).norm();
        hError = (analiticHessian-numHessian).norm();
    }
    if( print_to_screen)
    {
        cout << endl << endl;
        cout << "2 Points Outside Rot 1" << endl;
        cout << "Analitic Jacobian:\t" << analiticJacob.transpose() << endl;
        cout << "Numeric Jacobian:\t" << numJacob.transpose() << endl;
        cout << "Difference:\t" << (analiticJacob - numJacob).transpose() <<endl;
        cout << "Analitic Hessian:\t" << endl << analiticHessian << endl << endl;
        cout << "numeric Hessian:\t" << endl << numHessian << endl <<endl;
        cout << "Difference:\t" << endl << (analiticHessian-numHessian)<<endl<<endl;
        cout << "Analitic Jacobian PE: " << jError*100 << "%\tAnalitic  Hessian PE: " << hError*100 << "%"<<endl<<endl;
    }
    retVal &= jError < delta*5*6;
    retVal &= hError < delta*5*36;

    R = Eigen::AngleAxisd(.6,Eigen::Vector3d::UnitY())*Eigen::AngleAxisd(-.4,Eigen::Vector3d::UnitX())*Eigen::AngleAxisd(.2,Eigen::Vector3d::UnitZ());

    cost_ = cost(pos,R);
    analiticJacob = cost_jacobian(pos,R);
    analiticHessian = cost_hessian(pos,R);

    numJacob(0) = (cost(pos+Eigen::Vector3d(delta,0,0),R)-cost_)/delta;
    numJacob(1) = (cost(pos+Eigen::Vector3d(0,delta,0),R)-cost_)/delta;
    numJacob(2) = (cost(pos+Eigen::Vector3d(0,0,delta),R)-cost_)/delta;
    numJacob(3) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-cost_)/delta;
    numJacob(4) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-cost_)/delta;
    numJacob(5) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-cost_)/delta;

    numHessian.block(0,0,6,1) = (cost_jacobian(pos+Eigen::Vector3d(delta,0,0),R)-analiticJacob)/delta;
    numHessian.block(0,1,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,delta,0),R)-analiticJacob)/delta;
    numHessian.block(0,2,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,delta),R)-analiticJacob)/delta;
    numHessian.block(0,3,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-analiticJacob)/delta;
    numHessian.block(0,4,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-analiticJacob)/delta;
    numHessian.block(0,5,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-analiticJacob)/delta;

    jError = (analiticJacob - numJacob).norm()/analiticJacob.norm();
    hError = (analiticHessian-numHessian).norm()/analiticHessian.norm();
    if( cost_ == 0)
    {
        jError = (analiticJacob - numJacob).norm();
        hError = (analiticHessian-numHessian).norm();
    }
    if( print_to_screen)
    {
        cout << endl << endl;
        cout << "2 Points Outside Rot 2" << endl;
        cout << "Analitic Jacobian:\t" << analiticJacob.transpose() << endl;
        cout << "Numeric Jacobian:\t" << numJacob.transpose() << endl;
        cout << "Difference:\t" << (analiticJacob - numJacob).transpose() <<endl;
        cout << "Analitic Hessian:\t" << endl << analiticHessian << endl << endl;
        cout << "numeric Hessian:\t" << endl << numHessian << endl <<endl;
        cout << "Difference:\t" << endl << (analiticHessian-numHessian)<<endl<<endl;
        cout << "Analitic Jacobian PE: " << jError*100 << "%\tAnalitic  Hessian PE: " << hError*100 << "%"<<endl<<endl;
    }
    retVal &= jError < delta*5*6;
    retVal &= hError < delta*5*36;

    point_list.clear();
    point_list.push_back(Eigen::Vector3d(4.0*xLen/3.0+pos.x(),5.0*yLen/4.0+pos.y(), pos.z()+zLen/7.0));
    point_list.push_back(Eigen::Vector3d(-4.0*xLen/3.0+pos.x(),-5.0*yLen/4.0+pos.y(), pos.z()-zLen/7.0));

    R = Eigen::Matrix3d::Identity();

    cost_ = cost(pos,R);
    analiticJacob = cost_jacobian(pos,R);
    analiticHessian = cost_hessian(pos,R);

    numJacob(0) = (cost(pos+Eigen::Vector3d(delta,0,0),R)-cost_)/delta;
    numJacob(1) = (cost(pos+Eigen::Vector3d(0,delta,0),R)-cost_)/delta;
    numJacob(2) = (cost(pos+Eigen::Vector3d(0,0,delta),R)-cost_)/delta;
    numJacob(3) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-cost_)/delta;
    numJacob(4) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-cost_)/delta;
    numJacob(5) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-cost_)/delta;

    numHessian.block(0,0,6,1) = (cost_jacobian(pos+Eigen::Vector3d(delta,0,0),R)-analiticJacob)/delta;
    numHessian.block(0,1,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,delta,0),R)-analiticJacob)/delta;
    numHessian.block(0,2,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,delta),R)-analiticJacob)/delta;
    numHessian.block(0,3,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-analiticJacob)/delta;
    numHessian.block(0,4,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-analiticJacob)/delta;
    numHessian.block(0,5,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-analiticJacob)/delta;

    jError = (analiticJacob - numJacob).norm()/analiticJacob.norm();
    hError = (analiticHessian-numHessian).norm()/analiticHessian.norm();
    if( cost_ == 0)
    {
        jError = (analiticJacob - numJacob).norm();
        hError = (analiticHessian-numHessian).norm();
    }
    if( print_to_screen)
    {
        cout << endl << endl;
        cout << "2 Points Outside Rot 1" << endl;
        cout << "Analitic Jacobian:\t" << analiticJacob.transpose() << endl;
        cout << "Numeric Jacobian:\t" << numJacob.transpose() << endl;
        cout << "Difference:\t" << (analiticJacob - numJacob).transpose() <<endl;
        cout << "Analitic Hessian:\t" << endl << analiticHessian << endl << endl;
        cout << "numeric Hessian:\t" << endl << numHessian << endl <<endl;
        cout << "Difference:\t" << endl << (analiticHessian-numHessian)<<endl<<endl;
        cout << "Analitic Jacobian PE: " << jError*100 << "%\tAnalitic  Hessian PE: " << hError*100 << "%"<<endl<<endl;
    }
    retVal &= jError < delta*5*6;
    retVal &= hError < delta*5*36;

    R = Eigen::AngleAxisd(.6,Eigen::Vector3d::UnitY())*Eigen::AngleAxisd(-.4,Eigen::Vector3d::UnitX())*Eigen::AngleAxisd(.2,Eigen::Vector3d::UnitZ());

    cost_ = cost(pos,R);
    analiticJacob = cost_jacobian(pos,R);
    analiticHessian = cost_hessian(pos,R);

    numJacob(0) = (cost(pos+Eigen::Vector3d(delta,0,0),R)-cost_)/delta;
    numJacob(1) = (cost(pos+Eigen::Vector3d(0,delta,0),R)-cost_)/delta;
    numJacob(2) = (cost(pos+Eigen::Vector3d(0,0,delta),R)-cost_)/delta;
    numJacob(3) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-cost_)/delta;
    numJacob(4) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-cost_)/delta;
    numJacob(5) = (cost(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-cost_)/delta;

    numHessian.block(0,0,6,1) = (cost_jacobian(pos+Eigen::Vector3d(delta,0,0),R)-analiticJacob)/delta;
    numHessian.block(0,1,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,delta,0),R)-analiticJacob)/delta;
    numHessian.block(0,2,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,delta),R)-analiticJacob)/delta;
    numHessian.block(0,3,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitX())*R)-analiticJacob)/delta;
    numHessian.block(0,4,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitY())*R)-analiticJacob)/delta;
    numHessian.block(0,5,6,1) = (cost_jacobian(pos+Eigen::Vector3d(0,0,0),Eigen::AngleAxisd(delta,Eigen::Vector3d::UnitZ())*R)-analiticJacob)/delta;

    jError = (analiticJacob - numJacob).norm()/analiticJacob.norm();
    hError = (analiticHessian-numHessian).norm()/analiticHessian.norm();
    if( cost_ == 0)
    {
        jError = (analiticJacob - numJacob).norm();
        hError = (analiticHessian-numHessian).norm();
    }
    if( print_to_screen)
    {
        cout << endl << endl;
        cout << "2 Points Outside Rot 2" << endl;
        cout << "Analitic Jacobian:\t" << analiticJacob.transpose() << endl;
        cout << "Numeric Jacobian:\t" << numJacob.transpose() << endl;
        cout << "Difference:\t" << (analiticJacob - numJacob).transpose() <<endl;
        cout << "Analitic Hessian:\t" << endl << analiticHessian << endl << endl;
        cout << "numeric Hessian:\t" << endl << numHessian << endl <<endl;
        cout << "Difference:\t" << endl << (analiticHessian-numHessian)<<endl<<endl;
        cout << "Analitic Jacobian PE: " << jError*100 << "%\tAnalitic  Hessian PE: " << hError*100 << "%"<<endl<<endl;
    }
    retVal &= jError < delta*5*6;
    retVal &= hError < delta*5*36;

    pos.setZero();
    R.setIdentity();

    point_list.clear();
    if( constraint_type == KEEP_IN )
    {
        point_list.push_back(Eigen::Vector3d(xLen/2.0-.01,0,0));
        retVal &= cost(pos,R) == 0;
        point_list.push_back(Eigen::Vector3d(0,yLen/2.0-.01,0));
        retVal &= cost(pos,R) == 0;
        point_list.push_back(Eigen::Vector3d(0,0,zLen/2.0-.01));
        retVal &= cost(pos,R) == 0;

        point_list.push_back(Eigen::Vector3d(xLen/2.0+.01,0,0));
        retVal &= cost(pos,R) != 0;
        point_list.clear();
        point_list.push_back(Eigen::Vector3d(0,yLen/2.0+.01,0));
        retVal &= cost(pos,R) != 0;
        point_list.clear();
        point_list.push_back(Eigen::Vector3d(0,0,zLen/2.0+.01));
        retVal &= cost(pos,R) != 0;
        point_list.clear();
    }else if( constraint_type == KEEP_OUT )
    {
        point_list.push_back(Eigen::Vector3d(xLen/2.0-.01,0,0));
        retVal &= cost(pos,R) != 0;
        point_list.clear();
        point_list.push_back(Eigen::Vector3d(0,yLen/2.0-.01,0));
        retVal &= cost(pos,R) != 0;
        point_list.clear();
        point_list.push_back(Eigen::Vector3d(0,0,zLen/2.0-.01));
        retVal &= cost(pos,R) != 0;
        point_list.clear();

        point_list.push_back(Eigen::Vector3d(xLen/2.0+.01,0,0));
        retVal &= cost(pos,R) == 0;
        point_list.clear();
        point_list.push_back(Eigen::Vector3d(0,yLen/2.0+.01,0));
        retVal &= cost(pos,R) == 0;
        point_list.clear();
        point_list.push_back(Eigen::Vector3d(0,0,zLen/2.0+.01));
        retVal &= cost(pos,R) == 0;
    }


    point_list = pointListSave;
    return retVal;

}


double SuperEllipsoidKeepout::logistic(double x, double k) const
{
    double kx = k*x;

    if( kx > 100 )
        return x;
    if( kx < -100 )
        return 0;

    double mekx = std::exp(-kx);
    return x/(1.0+mekx);
}

double SuperEllipsoidKeepout::logisticDer(double x, double k) const
{
    double kx = k*x;

    if( kx > 100 )
        return 1;
    if( kx < -100 )
        return 0;

    double ekx = std::exp(kx);


    return ekx*(k*x+ekx+1.0)/std::pow(ekx+1.0,2);
}

double SuperEllipsoidKeepout::logisitcDer2(double x, double k) const
{
    double kx = k*x;

    if( kx > 100 )
        return 0;
    if( kx < -100 )
        return 0;

    double ekx = std::exp(kx);

    return -k*ekx*(-kx+ekx*(kx-2.0)-2.0)/std::pow(ekx+1.0,3);
}

#endif // SUPERELLIPSOID_KEEPOUT_HPP

