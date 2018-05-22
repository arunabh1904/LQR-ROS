cmake_minimum_required(VERSION 2.8.3)
project(lqr_obsavoid)

##find_package(catkin REQUIRED)

##Declare catkin and requried dependencies
##cmake_modules for Eigen
##cv_bridge for OpenCV
find_package(catkin REQUIRED 
 COMPONENTS roscpp rospy std_msgs geometry_msgs tf2 tf sensor_msgs velodyne_msgs nav_msgs cmake_modules
)

##Declare Eigen, PCL and OpenCV as required packages
find_package(Eigen REQUIRED)
find_package(PCL 1.7 REQUIRED)
find_package(OpenCV REQUIRED)

##Include directories
## (Include all paths to the headers and include catkin, Eigen, PCL and OpenCV libraries by function variables)
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}  
)


###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
##  INCLUDE_DIRS
##  /home/arunabh/CodeBase/
##  /home/arunabh/MPC/MPC_CPP_Development/
##  ${catkin_INCLUDE_DIRS}
##  ${Eigen_INCLUDE_DIRS}
##  ${PCL_INCLUDE_DIRS}
##  ${OpenCV_INCLUDE_DIRS}
  
  ##LIBRARIES geobot_navigation
##  CATKIN_DEPENDS 
##  roscpp  std_msgs  geometry_msgs  tf2  tf  sensor_msgs  velodyne_msgs  nav_msgs 

##  DEPENDS system_lib
##  Eigen
)


## Declare all C++ executables 
## (Declare all the cpp files for every class that are required to compile the code)
add_executable(lqr_obsavoid src/lqr_obsavoid.cpp
src/husky_collision_constraint.cpp
src/husky_system.cpp
src/superellipsoid_keepout.cpp
src/are_solver.cpp
src/Timer.cpp
src/integrationmethod.cpp
src/ivp_solver.cpp
src/precisionCalendar.cpp
src/PseudoInverse.cpp
)


## Link Libraries
## Link libraries to your project that are required to make the project run
target_link_libraries(lqr_obsavoid ${catkin_LIBRARIES} ${OpenCV_LIBS} ${PCL_LIBRARIES})



