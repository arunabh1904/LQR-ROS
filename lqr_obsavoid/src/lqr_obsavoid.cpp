#include "lqr_obsavoid.h"

#define LIDAR_DISPLACEMENT_DISTANCE -0.125 //0.125m
#define OBSTACLE_FROM_GROUND 0.02 // 2 cm
#define WHEEL_RADIUS 0.17775

#define PRINT_CONVERGED
//#define PRINT_DEBUG
//#define DONT_MOVE
//#define SOLVED_VELOCITY

LqrObsAvoid::LqrObsAvoid():obstacles(1/800.0, 2, 3, 2, 1.0625, Eigen::Vector3d(0,0,0))
{
  return;
}

void LqrObsAvoid::global_plan_cb(const std_msgs::Float64 vel)
{
  flag = true;
  target_velocity = vel.data;
}

void LqrObsAvoid::points_received_cb(const sensor_msgs::PointCloud2::ConstPtr& pCloud)
{
  float distThreshold = 10; // in meters...

  //clear obstacle list
  obstacle_list.clear();

  //cout<<"pcloud height: "<<pCloud->height<<endl;//1
  //cout<<"pcloud pointstep: "<<pCloud->point_step<<endl;//32
  //cout<<"pcloud frame id: "<<pCloud->header.frame_id;//frame - velodyne
  //cout<<"pcloud rowstep: "<<pCloud->row_step<<endl;//32*numPoints
  //cout<<"pcloud is_bigendian: "<<pCloud->is_bigendian;
  //cout<<"pcloud width: "<<pCloud->width;//numPoints
  world_obstacle_list.clear();

  int ignored  = 0;

  for (int numPoints = 0; numPoints<pCloud->row_step; numPoints+=32)
  {
    //declare pcl2 pointer for data copy
    const Pcl2Data *pcl2Data = (Pcl2Data*)(&pCloud->data[numPoints]);
    //calculate r
    float radial_distance = std::sqrt(pcl2Data->x * pcl2Data->x + pcl2Data->y * pcl2Data->y + pcl2Data->z * pcl2Data->z);

    //cout << "Ring no: " << pcl2Data->ring << "\tz: " << pcl2Data->z << endl;

    float angle = atan2(pcl2Data->y, pcl2Data->x);

    //Azimuthal filtering + removing dummy points on lidar + restricting field of view to forward 180 degrees
    if(pcl2Data->z < (-0.45) || radial_distance > distThreshold|| radial_distance < 0.001 || angle>5*M_PI/8.0 || angle<-5*M_PI/8.0)
    {
      ignored++;
      continue;
    }


    Eigen::Vector3d worldPoint = robot2World_Transform * Eigen::Vector4d(pcl2Data->x,pcl2Data->y,pcl2Data->z,1);
#ifdef PRINT_DEBUG
    //cout << "pcl point " << pclPoint.transpose() << "\tworld Point: " << worldPoint.transpose() << "\twp norm: " <<worldPoint.norm() << endl;
#endif
    world_obstacle_list.push_back(worldPoint);

#ifdef PRINT_DEBUG
    cout << "Number of obstacles: " << world_obstacle_list.size() << "\tNumber Ignored: " << ignored << endl;
#endif
  }
}

bool LqrObsAvoid::acquire_robotpose()
{
  //geometry_msgs::TransformStamped odom_base_transform;//odom2basetransform
  //tf::StampedTransform odom_base_transform;
  tf::Matrix3x3 mat;
  double roll;
  double pitch;
  double yaw;

  try
  {
    // this is needed for first frame incase the EKF has not started yet
    listener.waitForTransform("odom", "base_link", ros::Time(0), ros::Duration(1.0));
    listener.lookupTransform("odom", "base_link", ros::Time(0), odom_base_transform);
  }
  catch(tf::TransformException ex)
  {
    ROS_ERROR("%s",ex.what());
    ROS_INFO("odom transform not received");
    return false;
  }

  mat.setRotation(odom_base_transform.getRotation());

  for( int i=0;i<3;i++ )
    for( int j=0;j<3;j++ )
      robot2World_Transform(i,j) = mat[i][j]; // copy rotation matrix over to eigen version

  // copy over vector... Maybe we should have just used tf types... :)
  robot2World_Transform.block(0,3,3,1) = Eigen::Vector3d(odom_base_transform.getOrigin().x(),odom_base_transform.getOrigin().y(),odom_base_transform.getOrigin().z()+WHEEL_RADIUS);


  //mat.getRPY(roll, pitch, yaw);

  pose.x = odom_base_transform.getOrigin().getX();
  pose.y = odom_base_transform.getOrigin().getY();
  pose.heading = atan2(robot2World_Transform(1,0),robot2World_Transform(0,0));

  Eigen::Vector3d curr_state( pose.x, pose.y, pose.heading );

  obstacles.updateStateAndObstacle(curr_state, world_obstacle_list);

  return true;
}

int main(int argc, char **argv)
{

  /*****   SETUP ROS STUFF *******/
  //initialize ros node
  ros::init(argc, argv, "lqr_obsavoid");

  //create node handle
  ros::NodeHandle node_handle;

  LqrObsAvoid lqr_obsAvoid;
  //ROS_INFO("MAIN");

  // loop update rate
  double node_update_rate_Hz = 100;
  ros::Rate rate( node_update_rate_Hz );//setup ros rate (command line check: rostopic hz)

  //create velocity publisher object
  int velocity_command_que_size = 1;
  //ros::Publisher velocity_cmdPublisher = node_handle.advertise<geometry_msgs::Twist>("cmd_vel", velocity_command_que_size);
  ros::Publisher velocity_cmdPublisher = node_handle.advertise<geometry_msgs::Twist>("husky_velocity_controller/cmd_vel", velocity_command_que_size);


  //create subscriber object for pose
  //ros::Subscriber pose_Subscriber = node_handle.subscribe("odometry/filtered", pose_sub_rate_hz,&  LqrObsAvoid::pose_message_received_cb, &call_back);//odometry/filtered husky topic which publishes odometry data as nav_msgs/Odometry

  // create subscriber object for LIDAR
  int lidar_buffer_que_size = 1;
  //initialize lidar scan subscriber
  ros::Subscriber points_Subscriber = node_handle.subscribe("/velodyne_points", lidar_buffer_que_size,&  LqrObsAvoid::points_received_cb, &lqr_obsAvoid);//velodyne_points topic which publishes lidar data as sensor_msgs/PointCloud2

  int global_plan_buffer_que_size = 20;
  //initialize global plan subscriber
  ros::Subscriber globalPlan_Subscriber = node_handle.subscribe("/global_plan", global_plan_buffer_que_size,&  LqrObsAvoid::global_plan_cb, &lqr_obsAvoid);
  //*****************************************************************




  /*****************************   Initialization of MPC Objects   *****************************/

  Trajectory<3,2> pathTrajectory(  new HuskySystem() );//trajectory object

  // create constraints
  //Control_Saturation_Constraint<3,2> linVel_Constraint(0,1.0,-1.0);//lin vel saturation
  //Control_Saturation_Constraint<3,2> angVel_Constraint(1, M_PI/6.0, -M_PI/6.0);//ang vel saturation
  //call_back.obstacles(1/64.0, 2, 2.5, 1, 1.0625, Eigen::Vector3d(0,0,0));

  // create constraint list
  std::vector< Dynamic_System_Constraint<3,2>* > constraintList;//hard constraint list
  std::vector< Dynamic_System_Constraint<3,2>* > soft_ConstraintList;//soft constraint list
  // append constraints to list

  //constraintList.push_back(& linVel_Constraint);
  //constraintList.push_back(& angVel_Constraint);
  constraintList.push_back(& lqr_obsAvoid.obstacles);
  soft_ConstraintList.push_back(& lqr_obsAvoid.obstacles);


  // initialize MPC Weighting
  Eigen::MatrixXd R(2,2);//control penalty matrix
  R.setIdentity();
  //R *= 100;
  R(0,0) = 100000.0;//v
  R(1,1) = 90000.0;//w

  Eigen::MatrixXd Q(3,3);//error penalty matrix
  Q.setIdentity();
  Q(0,0) = Q(1,1) = 1.0;//x,y
  Q(2,2) = 5.0;//theta

  Eigen::MatrixXd C(3,3);//selection matrix
  C.setIdentity(); // pick off theta only

  double constant_target_velocity = 0;//desired target velocity
#ifdef DONT_MOVE
  constant_target_velocity = 0;
#endif

  double t = 1.0/node_update_rate_Hz;//time span
  Eigen::Vector3d des_target(0,0,0);//desired target pos and angle

  pathTrajectory.target_output = des_target;
  pathTrajectory.target_control = Eigen::Vector2d(constant_target_velocity,0);
  pathTrajectory.Q = Q;
  pathTrajectory.R = R;
  pathTrajectory.includeTargetControlInCost = false;
  pathTrajectory.t_span = T_Span(0, t);
  pathTrajectory.p_dynamic_system->setOutSelectMatrix(C);

  pathTrajectory.p_dynamic_system->solve(pathTrajectory.target_output, pathTrajectory.target_control, pathTrajectory.t_span, true);

  //cout << "Target Control: " << pathTrajectory.front().target_control.transpose() << endl;

  //    //get lidar points and set new obstacles in world frame into constraint list.
  //    call_back.acquire_robotpose();

  //    call_back.create_obstacle_list(call_back.obstacle_list, call_back.robot2World_Transform);

  //    obstacles.updateStateAndObstacle(Eigen::Vector3d(0,0,0), call_back.obstacle_list);


  //write a service client thing for this.
  //goal(std::Nan, std::Nan, std::Nan, std::Nan, std::Nan);

  /*****************************  END Initialization of MPC Objects   *****************************/



  /*******************   BEGIN LOOP *********************/

  while(ros::ok())
  {

    // Get robot pose
    lqr_obsAvoid.acquire_robotpose();

    // extract state from LqrObsAvoid object for clarity
    Eigen::Vector3d currentState( lqr_obsAvoid.pose.x, lqr_obsAvoid.pose.y, lqr_obsAvoid.pose.heading );

    // PRINT FOR DEBUG
    //cout<<"current state: "<<currentState<<endl;
    //for(int trajind=0;trajind<numPoints;trajind++)
    //{
    //cout<<"pathTrajectory.target_output:  "<<pathTrajectory[trajind].target_output<<endl;
    //cout<<"pathTrajectory.target_control:  "<<pathTrajectory[trajind].target_control<<endl;
    //}
    //cout<<"constraintList: "<<&constraintList<<endl;
    //cout<<"soft_ConstraintList: "<<&soft_ConstraintList<<endl;


    //call_back.create_world_obstacle_list();
#ifdef PRINT_DEBUG
    //cout << "world obs size: " << call_back.world_obstacle_list.size() << endl;
#endif

    // UPDATE TRAJECTORY DEFINITION
    pathTrajectory.target_control = Eigen::Vector2d(lqr_obsAvoid.target_velocity,0);
    //Calculating control gains
    Eigen::Vector3d error;

    pathTrajectory.target_output.head(2) = currentState.head(2);
    pathTrajectory.p_dynamic_system->solve(pathTrajectory.target_output, pathTrajectory.target_control, pathTrajectory.t_span, true);

    error = pathTrajectory.target_output - currentState;

    Eigen::MatrixXd Ad = pathTrajectory.p_dynamic_system->Jx_sens()+0.5*(pathTrajectory.p_dynamic_system->Jxx_sens()*error).rd_slice(0);
    Eigen::MatrixXd Bd = pathTrajectory.p_dynamic_system->Ju_sens()+(pathTrajectory.p_dynamic_system->Jxu_sens()*error).rd_slice(0);

    ARE_Solution sln = dare(Ad,Bd,pathTrajectory.Q,pathTrajectory.R);

    // obstacle avoidance
    Eigen::Vector3d Cx = lqr_obsAvoid.obstacles.x_sens(0,currentState,pathTrajectory.target_control,SOFT);

    error -= Cx;

    Eigen::Vector2d u = pathTrajectory.target_control+sln.K*error;

    //Velocity Saturation
    double w_max = .25;//M_PI/2.0;
    if( std::abs(u(1)) > w_max )
      u *= w_max/std::abs(u(1));

    //Control Commands
    double linear_velocity = u(0);//assigns vel in x dir
    double angular_velocity = u(1);//assigns ang vel about z dir

#ifdef SOLVED_VELOCITY
    cout<<"soln linear vel: "<<linear_velocity<<endl;
    cout<<"soln angular vel: "<<angular_velocity<<endl;
#endif

    //-----------------------------------------------------------------------------------------------


    //Define Velocity Command
    geometry_msgs::Twist velocity_command;//twist message to send velocity command.

    //set velocities to be sent to the husky
    velocity_command.linear.x =  linear_velocity;
    velocity_command.angular.z =  angular_velocity;

    //publish velocity command to husky.
    if(lqr_obsAvoid.flag)
      velocity_cmdPublisher.publish(velocity_command);

    lqr_obsAvoid.flag = false;
    // Call LqrObsAvoid functions for the velocity command published data
    ros::spinOnce();

    //flag check
    //cout << "Flag: " << lqr_obsAvoid.flag << endl;

    //wait for next iteration
    rate.sleep();
  }//ros::ok
  //shuts down the node gracefully (ctrl-c)

  geometry_msgs::Twist velocity_command;//twist message to send velocity command.

  //set velocities to be sent to the husky
  velocity_command.linear.x =  0;
  velocity_command.angular.z =  0;

  //publish velocity command to husky.
  velocity_cmdPublisher.publish(velocity_command);

  // Call LqrObsAvoid functions for the velocity command published data
  ros::spinOnce();

  ros::shutdown();
}//main

