#ifndef LQR_OBSAVOID_H
#define LQR_OBSAVOID_H

#include <ros/ros.h>
#include <eigen3/Eigen/Core>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <assert.h>
#include <stdio.h>
#include "math.h"

//************************************* Message_types includes***************

#include <geometry_msgs/TransformStamped.h> //for odom to base_link transform
#include <geometry_msgs/Twist.h> //twist type for velocity
#include <nav_msgs/Odometry.h> //husky pose twist data from odom filtered
#include <velodyne_msgs/VelodyneScan.h> //lidar scan as lidar packet
#include <sensor_msgs/PointCloud2.h> //pcl2 format to extract xyz points
#include <std_msgs/Float64.h> //velocity command type
//#include <sensor_msgs/point_cloud2_iterator.h> //iterator for reading data
#include <std_msgs/Float32.h>
//#include <velodyne_pointcloud/point_types.h>//POINTXYZIR
#include <sensor_msgs/PointField.h>

//************************************* tf2 includes***************

#include <tf2/transform_datatypes.h>//quaternion functions
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include "tf2_geometry_msgs/tf2_geometry_msgs.h"

//************************************* tf includes***************

#include <tf/transform_listener.h>
#include <tf/tf.h>
#include <tf/LinearMath/Quaternion.h>
#include <tf/LinearMath/Matrix3x3.h>

//************************************* PCL includes***************
#include "pcl_conversions/pcl_conversions.h"
#include "pcl_ros/point_cloud.h"
#include "pcl/point_types.h"

//************************************* Solver includes***************

#include "control_saturation_constraint.h"
#include "husky_collision_constraint.h"
//#include "husky_potential_field.h"
#include "husky_system.h"
//#include "mpc_solver.h"
#include "trajectory.h"
#include "are_solver.h"

class LqrObsAvoid
{
private:
  //struct for pose data
  struct Pose{
    double x;
    double y;
    double heading;
  };

  //struct for storing pointcloud2 data
  struct Pcl2Data{
    float x;
    float y;
    float z;
    float buffer;
    float intensity;
    short ring;
    char buffer2[10];
  };

  //initialize the listener
  tf::TransformListener listener;//initialization should never be in a callback

public:
  LqrObsAvoid();

  //flag to check if global plan received
  bool flag;

  Eigen::Matrix<double, 3, 4> robot2World_Transform;//identity at first, keep updating for every frame as we go.

  tf::StampedTransform odom_base_transform;

  Husky_Collision_Constraint obstacles;//collision constraint object

  //robot2World_Transform.setIdentity();
  std::vector<Eigen::Vector2d> obstacle_list;//list for storing obstacle points

  std::vector<Eigen::Vector3d> world_obstacle_list;//obstacle points

  Pose pose;

  double target_velocity;

  //to get current pose data
  bool acquire_robotpose();

  //Callback for lidar points
  void points_received_cb(const sensor_msgs::PointCloud2::ConstPtr& pcloud);

  //Global planner callback
  void global_plan_cb(const std_msgs::Float64 vel);
};

#endif //LQR_OBSAVOID_H
