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
const float PI = 3.141567;

mavros_msgs::State current_state;
void state_cb(const mavros_msgs::State::ConstPtr& msg){
  current_state=*msg;
}

geometry_msgs::PoseStamped current_pos;
void current_pos_cb(const geometry_msgs::PoseStamped::ConstPtr& msg1){
  current_pos=*msg1;
}

geometry_msgs::PoseStamped offset;
bool offset_check_flag = false;

void ENUoff_cb(const geometry_msgs::PoseStamped::ConstPtr& msg1){
  offset=*msg1;
  offset_check_flag = true;
}


int main(int argc, char **argv)
{
  //std::string id="1";
  ros::init(argc, argv, "hover_node");
  ros::NodeHandle nh;
  tf::TransformBroadcaster broadcaster;
  tf::Transform quad_body_frame(tf::Transform::getIdentity());
  
  int id;
  nh.param<int>("id", id, 0); //id of the quadrotor
  
//M-Air Dimensions (The local common frame follows ENU convection) with the origin close to the small door near the pavilion 

float X_max=21, X_min=-1, Y_max=-2.5, Y_min=-36;



//Subscriber and Publisher Block

  ros::Subscriber state_sb = nh.subscribe<mavros_msgs::State>("mavros/state", 10, state_cb);
  ros::Subscriber curr_pos = nh.subscribe<geometry_msgs::PoseStamped>("gstation_position",20, current_pos_cb);
  ros::Subscriber curr_off_sb = nh.subscribe<geometry_msgs::PoseStamped>("local_ENU_offset", 10, ENUoff_cb);
 
  ros::Publisher setpoint_pub = nh.advertise<geometry_msgs::PoseStamped>("desired_setpoint", 10);

  ros::Rate rate(20.0);

  // get the parameters for the circular trajectory from the user otherwise set some default values
  double circle_center_x ,circle_center_y, radius;
  
  // Make the drones hover for hover_time seconds
    double hover_time, hover_altitude;

    nh.param<double>("hover_time", hover_time, 25);
    nh.param<double>("hover_atitude", hover_altitude, 1.5);

  ros::Time last_request = ros::Time::now(),last_point;
  int runFlag=0;

// Set the current position in the gstation ENU frame, which is the local_ENU_offset, with a given hover_altitude as the setpoint
  geometry_msgs::PoseStamped setpoint;
  setpoint.header.frame_id = "world";

  while(offset_check_flag==false)
  {
  ros::spinOnce();
  }  

  ROS_INFO_STREAM("Received local_ENU_offset message for Quadrotor "<<id);

  setpoint.pose.position.x = offset.pose.position.x; 
  setpoint.pose.position.y = offset.pose.position.y;
  setpoint.pose.position.z = offset.pose.position.z + hover_altitude ; //hover at the current location at 1.5 m heigher than the current location
  

  std::string quad_frame_string;
  quad_frame_string = "quad";
  quad_frame_string.append(std::to_string(id));

  //Send hover commands
  double t0 = ros::Time::now().toSec();
  double t = ros::Time::now().toSec();
   ROS_INFO_STREAM("Sending hover position commands to Quadrotor "<<id<<" for "<<hover_time<<" seconds");
   while ((t-t0)<hover_time)
  {  	    	
  	setpoint_pub.publish(setpoint); 

    quad_body_frame.setOrigin(tf::Vector3(setpoint.pose.position.x, setpoint.pose.position.y, setpoint.pose.position.z));
    broadcaster.sendTransform(tf::StampedTransform(quad_body_frame, ros::Time::now(), "world", quad_frame_string)); 
  	
  	ros::spinOnce();
    rate.sleep();          
    t = ros::Time::now().toSec();
  }
  ROS_INFO_STREAM("Hover completed for Quadrotor "<<id);
  nh.setParam("flag_hover_completed",1);
 
  //ros::shutdown();
}