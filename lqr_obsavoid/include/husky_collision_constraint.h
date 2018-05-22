#ifndef HUSKY_COLLISION_CONSTRAINT_H
#define HUSKY_COLLISION_CONSTRAINT_H

#include "dynamic_system_constraint.h"
#include "superellipsoid_keepout.h"

#include <vector>

// All units in meters!!!

class Husky_Collision_Constraint: public Dynamic_System_Constraint<3,2>
{
public:
    Husky_Collision_Constraint(double potential_field_mag, double potential_field_exp, double husky_length = 1.0, double husky_width = 0.7, double husky_height = 0.85, const Eigen::Matrix<double,3,1>& sensor_location = Eigen::Vector3d(0.125,0,0)); ///<- defaults based on a slightly bigger measured envelope

    void updateStateAndObstacle(const Eigen::Matrix<double,3,1>& newState, const std::vector< Eigen::Matrix<double,3,1> >& newObstacle ); ///<- Updates our understanding of the robot at measurement time and the obstacle locations
    void updateStateAndObstacle(const Eigen::Matrix<double,3,1>& newState, const std::vector< Eigen::Matrix<double,2,1> >& newObstacle ); ///<- Updates our understanding of the robot at measurement time and the obstacle locations

    double constraint_cost(double t, const Eigen::Matrix<double,3,1>& x, const Eigen::Matrix<double,2,1>& u, Dynamic_System_Constraint_Type type) const; ///<- The hard constraint. Should be zero when satisfied
    Eigen::Matrix<double,3,1> x_sens(double t, const Eigen::Matrix<double,3,1>& x, const Eigen::Matrix<double,2,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the derivative of the constraint with respect to state
    Eigen::Matrix<double,3,3> xx_sens(double t, const Eigen::Matrix<double,3,1>& x, const Eigen::Matrix<double,2,1>& u, Dynamic_System_Constraint_Type type) const; ///<- the second derivative of the constraint with respect to state^2

    void setPotentialFieldMag( double newMag );
    void setPotentialFieldExponent( double newExponent );

private:
    SuperEllipsoidKeepout husky_envelope;
    SuperEllipsoidKeepout husky_envelope_soft;

    double husky_half_height;
    double eps;
    double n_pow;
    double a;

    Eigen::Matrix<double,3,1> measurement_state;
    Eigen::Matrix<double,3,1> sensor_location;
};

#endif // HUSKY_COLLISION_CONSTRAINT_H
