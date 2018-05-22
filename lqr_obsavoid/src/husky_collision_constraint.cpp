#include "husky_collision_constraint.h"

Husky_Collision_Constraint::Husky_Collision_Constraint(double potential_field_mag, double potential_field_exp, double husky_length, double husky_width, double husky_height, const Eigen::Matrix<double,3,1>& sensor_location_input):
    Dynamic_System_Constraint<3,2>(true,true,false,false,false),
    husky_envelope(husky_length,husky_width,husky_height,2,2,SuperEllipsoidKeepout::KEEP_OUT),
    husky_envelope_soft(potential_field_mag,potential_field_exp,husky_length*(1.0-1e-2),husky_width*(1.0-1e-2),husky_height*(1.0-1e-2),2,2,SuperEllipsoidKeepout::KEEP_OUT_POTENTIAL_FIELD)
{
    sensor_location = sensor_location_input;
    husky_half_height = husky_height/2.0;

    return;
}


void Husky_Collision_Constraint::setPotentialFieldMag( double newMag )
{
    husky_envelope.setMagnitude(newMag);
    husky_envelope_soft.setMagnitude(newMag);
}

void Husky_Collision_Constraint::setPotentialFieldExponent( double newExponent )
{
    husky_envelope.setExponent( newExponent);
    husky_envelope_soft.setExponent(newExponent);
}

void Husky_Collision_Constraint::updateStateAndObstacle(const Eigen::Matrix<double,3,1>& newState, const std::vector< Eigen::Matrix<double,3,1> >& newObstacle )
{
    // Assign measurement_state with offset to account for the difference in lidar placement to center of robot

    std::vector< Eigen::Matrix<double,3,1> > shifted_obstacle_data = newObstacle;

    measurement_state = newState;

    for (int i = 0; i < shifted_obstacle_data.size();i++)
    {
        shifted_obstacle_data[i] += sensor_location;
    }

    husky_envelope.setPointList(shifted_obstacle_data);
    husky_envelope_soft.setPointList(shifted_obstacle_data);

    return;
}

void Husky_Collision_Constraint::updateStateAndObstacle(const Eigen::Matrix<double,3,1>& newState, const std::vector< Eigen::Matrix<double,2,1> >& newObstacle )
{
    // Assign measurement_state with offset to account for the difference in lidar placement to center of robot

    std::vector< Eigen::Matrix<double,3,1> > shifted_obstacle_data;

    for (int i = 0; i < newObstacle.size();i++)
    {
        Eigen::Vector3d newObs;
        newObs.x() = newObstacle[i].x()+sensor_location.x();
        newObs.y() = newObstacle[i].y()+sensor_location.y();
        newObs.z() = husky_half_height+sensor_location.z();
        shifted_obstacle_data.push_back(newObs);
    }

    husky_envelope.setPointList(shifted_obstacle_data);
    husky_envelope_soft.setPointList(shifted_obstacle_data);

    measurement_state = newState;

    return;
}


double Husky_Collision_Constraint::constraint_cost(double, const Eigen::Matrix<double,3,1>& x, const Eigen::Matrix<double,2,1>&, Dynamic_System_Constraint_Type type) const
{
    //<- The hard constraint. Should be zero when satisfied

    // Calculate Robot rotation relative to measurement frame 3x3 matrix

    Eigen::Matrix<double,3,3> rot = Eigen::Matrix<double,3,3>::Identity();

    rot(0,0) = std::cos(x(2)-measurement_state(2));
    rot(1,0) = std::sin(x(2)-measurement_state(2));
    rot(0,1) = -rot(1,0);
    rot(1,1) = rot(0,0);

    // Calculate Robot position relative to measurement frame

    Eigen::Matrix<double,3,1> pos = Eigen::Matrix<double,3,1>::Zero();

    pos(0) = x(0)-measurement_state(0);
    pos(1) = x(1)-measurement_state(1);
    pos(2) = husky_half_height;

    // cost call...
    double cost = 0;
    if(type==HARD)
    {
        //cout<<"hard cost : "<< husky_envelope.cost(pos, rot);
        cost = husky_envelope.cost(pos, rot);
    }
    else
    {
        //cout<<"soft cost : "<< husky_envelope_soft.cost(pos, rot);
        cost = husky_envelope_soft.cost(pos, rot);
    }

    return cost;
}

Eigen::Matrix<double,3,1> Husky_Collision_Constraint::x_sens(double, const Eigen::Matrix<double,3,1>& x, const Eigen::Matrix<double,2,1>&, Dynamic_System_Constraint_Type type) const
{
    //<- the derivative of the constraint with respect to state

    // Calculate Robot rotation relative to measurement frame 3x3 matrix

    Eigen::Matrix<double,3,3> rot = Eigen::Matrix<double,3,3>::Identity();

    rot(0,0) = std::cos(x(2)-measurement_state(2));
    rot(1,0) = std::sin(x(2)-measurement_state(2));
    rot(0,1) = -rot(1,0);
    rot(1,1) = rot(0,0);

    // Calculate Robot position relative to measurement frame

    Eigen::Matrix<double,3,1> pos = Eigen::Matrix<double,3,1>::Zero();

    pos(0) = x(0)-measurement_state(0);
    pos(1) = x(1)-measurement_state(1);
    pos(2) = husky_half_height;

    Eigen::Matrix<double,6,1> x_sens_6D;

    if(type==HARD)
    {
        x_sens_6D = husky_envelope.cost_jacobian(pos, rot);
    }
    else
    {
        x_sens_6D = husky_envelope_soft.cost_jacobian(pos, rot);
    }

    Eigen::Matrix<double,3,1> x_sens_return;
    x_sens_return(0) = x_sens_6D(0); // motion in x
    x_sens_return(1) = x_sens_6D(1); // motion in y
    x_sens_return(2) = x_sens_6D(5); // motion in theta

    return x_sens_return;

}


Eigen::Matrix<double,3,3> Husky_Collision_Constraint::xx_sens(double, const Eigen::Matrix<double,3,1>& x, const Eigen::Matrix<double,2,1>&, Dynamic_System_Constraint_Type type) const
{
    //<- the second derivative of the constraint with respect to state^2

    // Calculate Robot rotation relative to measurement frame 3x3 matrix

    Eigen::Matrix<double,3,3> rot = Eigen::Matrix<double,3,3>::Identity();

    rot(0,0) = std::cos(x(2)-measurement_state(2));
    rot(1,0) = std::sin(x(2)-measurement_state(2));
    rot(0,1) = -rot(1,0);
    rot(1,1) = rot(0,0);

    // Calculate Robot position relative to measurement frame

    Eigen::Matrix<double,3,1> pos = Eigen::Matrix<double,3,1>::Zero();

    pos(0) = x(0)-measurement_state(0);
    pos(1) = x(1)-measurement_state(1);
    pos(2) = husky_half_height;

    Eigen::Matrix<double,6,6> xx_sens_6D;

    if(type==HARD)
    {
        xx_sens_6D = husky_envelope.cost_hessian(pos, rot);
    }
    else
    {
        xx_sens_6D = husky_envelope_soft.cost_hessian(pos, rot);
    }


    Eigen::Matrix<double,3,3> xx_sens_return;

    xx_sens_return.topLeftCorner(2,2) = xx_sens_6D.topLeftCorner(2,2);
    xx_sens_return(2,2) = xx_sens_6D(5,5);//theta,theta

    xx_sens_return.topRightCorner(2,1) = xx_sens_6D.topRightCorner(2,1);
    xx_sens_return.bottomLeftCorner(1,2) = xx_sens_6D.bottomLeftCorner(1,2);

    return  xx_sens_return;
}
