#include<ros/ros.h>
#include<ros/console.h>
//#include<mavros_msgs/SwarmCommands.h>
#include <eigen_conversions/eigen_msg.h>
#include<geometry_msgs/PoseStamped.h>
#include<mavros_msgs/CommandBool.h>
#include<mavros_msgs/SetMode.h>
#include<mavros_msgs/State.h>
#include<mavros_msgs/GlobalPositionTarget.h>
#include<mavros_msgs/HomePosition.h>
#include <mavros/mavros_plugin.h>
#include <sensor_msgs/NavSatFix.h>
#include <sensor_msgs/NavSatStatus.h>
#include "std_msgs/String.h"
#include <geographic_msgs/GeoPointStamped.h>
#include <geographic_msgs/GeoPoseStamped.h>
#include <geometry_msgs/TransformStamped.h>
#include <GeographicLib/Geocentric.hpp>
#include<string>
#include <math.h>
#include <cmath>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>

#include <random>
#include <iostream>
#include <memory>
#include <functional>
#include "boost/bind.hpp"

int N = 4;
// mavros_msgs::State current_state;
std::vector<mavros_msgs::State> current_state_list(N);
void state_cb(const mavros_msgs::State::ConstPtr& msg, mavros_msgs::State current_state){
  current_state=*msg;
}

// geometry_msgs::PoseStamped current_pos;
std::vector<geometry_msgs::PoseStamped> current_pos_list(N);
void current_pos_cb(const geometry_msgs::PoseStamped::ConstPtr& msg1, geometry_msgs::PoseStamped current_pos){
  current_pos=*msg1;
}

// geometry_msgs::PoseStamped offset;
std::vector<geometry_msgs::PoseStamped> offset_list(N);
bool offset_check_flag = false;

void ENUoff_cb(const geometry_msgs::PoseStamped::ConstPtr& msg1, geometry_msgs::PoseStamped offset){
  offset=*msg1;
  offset_check_flag = true;
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "circle_traj_node");
  ros::NodeHandle nh;
  // tf::TransformBroadcaster broadcaster;
  // tf::Transform quad_body_frame(tf::Transform::getIdentity());
  ros::Subscriber state_sb = nh.subscribe<mavros_msgs::State>("/quad1/mavros/state", 10, boost::bind(state_cb,boost::placeholders::_1,current_state_list[0]));
  std::vector<tf::Transform> body_frame_quad(N);
  for (int i = 0; i < N; i++)
  {
      body_frame_quad[i]=tf::Transform::getIdentity();
  }

}
