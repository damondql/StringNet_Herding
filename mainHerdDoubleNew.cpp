#include "AllParametersExperiment.cpp"
// #include "symDerivative.cpp"
#include "findCommGraphAndFormDist.cpp"
#include "defInitDesirePos.cpp"
#include "controlAttacker4.cpp"
#include "controlDefender5.cpp"
#include "defDesiredOpenForm.cpp"
#include "controlFiniteTimeTrajTracking.cpp"
#include "inhull.cpp"
#include "controlDefenderFormation4.cpp"
#include "defDesiredClosedForm.cpp"
#include "modifiedDIDynamics.cpp"
#include "findCoordOnPath.cpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <vector>
#include <list>
#include <deque>
#include <iterator>
#include <chrono>
#include "boost/bind.hpp"
#include "boost/function.hpp"


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

using namespace std;
using namespace arma;

// mavros_msgs::State current_state;
std::vector<mavros_msgs::State> current_state_list;
void state_cb(const mavros_msgs::State::ConstPtr& msg, int i){
  current_state_list[i]=*msg;
}

// geometry_msgs::PoseStamped current_pos;
std::vector<geometry_msgs::PoseStamped> current_pos_list;
void current_pos_cb(const geometry_msgs::PoseStamped::ConstPtr& msg1, int i){
  current_pos_list[i]=*msg1;
}

// geometry_msgs::PoseStamped offset;
std::vector<geometry_msgs::PoseStamped> offset_list;
bool offset_check_flag = false;

void ENUoff_cb(const geometry_msgs::PoseStamped::ConstPtr& msg1, int i){
  offset_list[i]=*msg1;
  offset_check_flag = true;
}


int flagExp = 1;
int flagGazeboExp = 0;


double t = 0;
auto Time = 0;
auto t_max = 2000;
auto endGame = 0;
double Niter;
mat vA_dot;
mat X;
mat XA;
//Allocate memory
void AllocateMemory(int NA, int ND) {
    vA_dot.resize(2,NA);
    XA.resize(4,NA);
    Niter = floor(t_max/dt);
    X.set_size(4*(ND+NA+1), Niter+1);
    X.fill(0);
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < NA; j++)
        {
            XA(i,j) = XA0(i,j);
        }
        
    }
    for (size_t i = 0; i < XD0.n_rows; i++)
    {
        for (size_t j = 0; j < XD0.n_cols; j++)
        {
            XD(i,j) = XD0(i,j);
        }
        
    }
    
}

mat Cov_YA(4,4,fill::zeros);
mat YA;
mat SD;
mat SD_arr;
mat rAcm;
mat vAcm;
mat rDcm;
mat vDcm;
double flagNDeven;
void measurements(int NA, int ND) {
    Cov_YA = {{pow(0.9/3,2),0,0,0},
                  {0,pow(0.9/3,2),0,0},
                  {0,0,pow(0.12,2),0},
                  {0,0,0,pow(0.12,2)}};
    // YA = XA + arma::mvnrnd(arma::zeros(4,1), Cov_YA, NA);
    // YA.load("../../../../../Downloads/swarm_matlab/controlD/YA.txt");
    YA = XA;
    // XD0.load("../../../../../Downloads/swarm_matlab/controlD/XD0.txt");
    SD = arma::zeros(ND,1);
    SD_arr = arma::zeros(ND,Niter+1);
    rAcm = arma::sum(YA,1)/NA;
    rAcm = rAcm(2,0, arma::size(2,1));
    rDcm = arma::sum(XD,1)/ND;
    vDcm = rDcm(2,0,arma::size(2,1));
    rDcm = rDcm(0,0, arma::size(2,1));
    flagNDeven = int(ND) % 2;
}


//Initialization of vectors and matrices
//Initial control
mat uA;
mat RPA_arr,RAD_arr;
mat Psi, Psi_prime, Psi_des, Psi_des_dot, Psi_des_ddot, Psi_prime_dot, Psi_prime_ddot;
mat vA_dot_arr, FA_dot_arr, XD_Des, U, times, minEAO, minEDO, minRDD, minRAD, minRAA, arr_norm_F_AO;
mat AOut;
double NAout = 0;
double AIcount = 0; //counter for Attacker interceptions;
double safe_flag=0;
double RPA, epsilon, RAD, RA0_des, Rii0, Rii1;
double RDF_open = 0;
double R_DD_string, RDF_closed;
// Initial string chain
int flagHerd = 0;
int flagDForm = 0;
mat WDString;
CommGraph attacker_graph;
CommGraph defender_graph_close;
CommGraph defender_graph_open;
void initial_contorl(int NA, int ND) {
    uA = vA;
    RPA_arr = arma::zeros(1, Niter+1);
    RAD_arr=arma::zeros<mat>(1,Niter+1);
    Psi=arma::zeros<mat>(1,Niter+1);
    Psi_prime=arma::zeros<mat>(1,Niter+1);
    Psi_des=arma::zeros<mat>(1,Niter+1);
    Psi_des_dot=arma::zeros<mat>(1,Niter+1);
    Psi_des_ddot=arma::zeros<mat>(1,Niter+1);
    Psi_prime=arma::zeros<mat>(1,Niter+1);
    Psi_prime_dot=arma::zeros<mat>(1,Niter+1);
    Psi_prime_ddot=arma::zeros<mat>(1,Niter+1);
    vA_dot_arr=arma::zeros<mat>(2*NA,Niter+1);
    FA_dot_arr=arma::zeros<mat>(2*NA,Niter+1);
    XD_Des=arma::zeros<mat>(4*ND,Niter+1);
    U=arma::zeros<mat>(2*(ND+NA+1),Niter+1);
    times=arma::zeros<mat>(1,Niter+1);
    minEAO=arma::zeros<mat>(1,Niter+1);
    minEDO=minEAO;
    minRDD=minEAO;
    minRAD=minEAO;
    minRAA=minRAD;
    arr_norm_F_AO=minEAO;
    std::complex<double> comp = {rA(0,0)-rP[0],rA(1,0) - rP[1]};
    RPA=sqrt(std::norm(comp));
    epsilon=M_PI/100;
    RAD=RAD_max;
    RA0_des=0;
    Rii0 = RA0*sqrt(2*(1-cos(2*M_PI/NA)));
    Rii1=Rii0*sqrt(2*(1-cos(180-2*M_PI/NA)));
    attacker_graph = findCommGraphAndFormDist(NA,1,RA);
    defender_graph_close = findCommGraphAndFormDist(ND+1, 2, rho_sn);
    // double RDF_open = 0;
    if (flagExp == 1) {
        RDF_open = 4.4;
    } else if (flagGazeboExp == 1)
    {
        RDF_open = 1.5 * rho_sn;
    } else {
        RDF_open = rho_sn + 20;
    }
    defender_graph_open = findCommGraphAndFormDist(ND+1, 4, RDF_open);
    R_DD_string = 1.5 * defender_graph_close.Rij_tilde(0,1);
    RDF_closed = 1.1 * rho_sn;
    WDString = arma::zeros<mat>(ND, ND);
    // cout  << "RDF_closed: " << RDF_closed << endl;
    // cout << "RDF_open" <<RDF_open << endl;
}

DesiredPos motionP_result;

void getMotionPlan(int NA, int ND){
    motionP_result =  defInitDesiredPos(YA, XD0, NA, ND, RDF_open, v_maxA[0], 105);
    motionP_result.dDf.phi += M_PI;
    motionP_result.mP.startTime += ones<mat>(ND,1) * 25;
    // cout << "finish motionP_result" << endl;
    // XD.print("XD:");
    arma::mat tempM;
    tempM.reshape(XD.n_rows, XD.n_cols+1);
    tempM.submat(0,0,XD.n_rows-1,XD.n_cols-1) = XD;
    tempM.submat(0,tempM.n_cols-1,motionP_result.dDf.rDFc0.n_rows-1, tempM.n_cols-1) = motionP_result.dDf.rDFc0;
    tempM.submat(motionP_result.dDf.rDFc0.n_rows,tempM.n_cols-1, tempM.n_rows-1,tempM.n_cols-1) = zeros<mat>(2,1);
    XD = tempM;
    // XD.print("XD:");
    tempM.resize(X.n_rows,1);
    tempM.submat(0,0,XA.n_cols*XA.n_rows-1,0) = XA.as_col();
    tempM.submat(XA.n_cols*XA.n_rows,0,tempM.n_rows-1,0) = XD.as_col();
    // tempM.print("tempM:");
    X.col(0) = tempM;
    // X.submat(0,0,19,9).print("X sub:");

}

vec indDef;
void calDistance(int ND){
    //Defenders indices for the formation
    indDef = regspace(1,ND+1);
    int na = motionP_result.mP.assign.n_elem;
    int nid = indDef.n_elem;
    vec tempV(indDef.n_elem);
    int indx_num = na + regspace(na+1,nid).n_elem;
    vec indx(indx_num);
    indx.subvec(0,na-1) = motionP_result.mP.assign;
    indx.subvec(na,indx_num-1) = regspace(na+1,nid);
    
    for (int i = 0; i < tempV.n_elem; i++)
    {
        tempV(indx(i)-1) = indDef(i);
    }
    indDef = tempV;
    // cout << "finish first for loop" << endl;

    //minimum distance between the defenders
    vec RDjDl(ND);
    vec arr_minRDD(ND-1);
    vec arr_minRAD(ND-1);
    // cout << "start sencond for loop" << endl;
    if (ND >1)
    {
        for (int j = 0; j < ND-1; j++)
        {
            for (int i = j+1 ; i < ND; i++)
            {
                RDjDl(i) = norm(rD.col(j) - rD.col(i));
            }
            arr_minRDD(j) = RDjDl.subvec(j+1,RDjDl.n_rows-1).min();
            mat tempM(rA.n_rows, rA.n_cols, fill::zeros);
            for (int i = 0; i < rA.n_cols; i++)
            {
                tempM.col(i) = rD.col(j) - rA.col(i);
            }
            
            arr_minRAD(j) = norm(tempM);
        }
        minRDD(0,0) = arr_minRDD.min();
        minRAD(0,0) = arr_minRAD.min();
        
    }
    

}

double dpsiAi0=0.02;
double vA_ref_max=4;
extern double dpsi_var;
int flagDefReachOpen=0;
int flagDefReachClosed=0;
int flagDefConnect=0;
int flagAttInSight=0;
mat defReachCount;
mat flagAttInObs;
mat FlagDefInObs;
mat betaAv0;
mat betaDv0;
mat mDv0;
mat cDv0;
mat mAv0;
mat cAv0;
mat speedA0;
mat assignment; // needs to be initialized later
int distTol=3;
int countDDes=0;

int flagGather=1;
int flagSeek=0;
int flagEnclose=0;
int flagAttackerStayTogether=1;

void control_loop(int NA, int ND, double z_h,
                  std::vector<ros::Subscriber> state_sb_list,
                  std::vector<ros::Subscriber> curr_pos_list,
                  std::vector<ros::Subscriber> curr_off_sb_list,
                  std::vector<ros::Publisher> setpoint_pub_list,
                  ros::Rate rate,
                  tf::TransformBroadcaster broadcaster,
                  std::deque<tf::Transform> body_frame_quad){
    assignment = regspace(1,1,ND);
    defReachCount.resize(ND,1);
    flagAttInObs.resize(NA,NO);
    FlagDefInObs.resize(ND,NO);
    betaAv0.resize(NA,NO);
    betaDv0.resize(ND,NO);
    mDv0=betaDv0;
    cDv0=mDv0;
    mAv0=betaAv0;
    cAv0=mAv0;
    speedA0=betaAv0;
    mat RefTraj(4*ND, Niter);
    mat SigmaProdD_arr;
    mat rAcm_arr;
    mat vAcm_arr;
    mat sigmaProd;
    mat minRAAProjS(Niter,1);
    int out_ti;
    int bound = Niter;
    mat uD;
    mat XD_des;
    mat XD_des_dot;
    mat uDFc_trans;
    mat XDF_des;
    cube WDString_mat(ND,ND, Niter);
    geometry_msgs::PoseStamped setpoint;
    setpoint.header.frame_id = "world";
    auto start = std::chrono::high_resolution_clock::now();
    for (int ti = 0; ti < bound; ti++)
    {
        mat Psi(NA,1, fill::zeros);
        mat Psi_dot;
        mat Psi_ddot;
        Psi_dot = Psi;
        Psi_ddot = Psi;
        
        int countAinS = 0;
        for (int ii = 0; ii < NA; ii++)
        {
            if (arma::norm(rA.col(ii) - rS) < rho_S)
            {
                countAinS++;
            }
        }

        int countDinS = 0;
        for (int jj = 0; jj < ND; jj++)
        {
            if (arma::norm(rD.col(jj) - rS) < rho_S)
            {
                countDinS++;
            }
        }
        if (countAinS > NA-1 && countDinS > ND-1)
        {
            break;
        }
        //////////////////////////////////////////
        /////////   safe flag protion     ////////
        //////////////////////////////////////////

        mat xd0(XD.n_rows, ND);
        for (int i = 0; i < ND; i++)
        {
            xd0.col(i) = XD.col(i);
        }
        XD_Des.col(ti+1) = xd0.as_col();
        
        vec pairDA(ND);
        pairDA = regspace(1,1,ND);
        double kref = 0.1;

        mat refTraj(4,ND, fill::ones);
        RefTraj.col(ti) = refTraj.as_col();


        vec XA_lead_goal(4,fill::zeros);
        XA_lead_goal.subvec(0,1) = rP;
        mat XA_goal(XA_lead_goal.n_rows, NA);
        mat XA_goal_dot(XA_lead_goal.n_rows, NA);
        XA_goal.col(0) = XA_lead_goal;
        XA_goal_dot.col(0) = XA_lead_goal;

        // if(ti == bound -1)
        // {   
        //     cout << "controlAttacker input: " << endl;
        //     cout << "ti: " << ti << endl;
        //     XA.print("XA: ");
        //     XA_goal.print("XA_goal: ");
        //     XA_goal_dot.print("XA_goal_dot: ");
        //     cout << "flagEnClose: " << flagEnclose << endl;
        //     cout << "flagHerd: " << flagHerd << endl;
        //     XD.print("XD: ");
        //     attacker_graph.W.print("W: ");
        //     WDString.print("WDString: ");
        //     cout << endl;
        // }

        control_attacker_t control_A_result = controlAttacker4(XA,XA_goal,XA_goal_dot,flagEnclose, flagHerd, XD,attacker_graph.W,WDString,NA,ND, ti, bound);
        uA = control_A_result.uA;
        // if(ti == bound -1)
        // {
        //     cout << "ti: " << ti << endl;
        //     control_A_result.uA.print("uA:");
        //     control_A_result.uA0.print("uA0: ");
        //     control_A_result.F_A.print("F_A: ");
        //     control_A_result.F_A_dot.print("F_A_dot: ");
        //     cout << "R_AO_min: " << control_A_result.R_AO_min << endl;
        //     cout << "R_AAProjS_min: " << control_A_result.R_AAProjS_min << endl;
        //     control_A_result.vA_des.print("vA_des");
        //     control_A_result.vA_des_dot.print("vA_des_dot: ");
        //     cout << endl;
        // }
        
        sigmaProd = control_A_result.SigmaProdD;
        // cout << "controlAttacker4:" << endl;
        // sigmaProd.print("sigmaProd:");
        SigmaProdD_arr.resize(control_A_result.SigmaProdD.n_rows,ti+1);
        SigmaProdD_arr.col(ti) = control_A_result.SigmaProdD;
        // cout << "SigmaProd_arr" << endl;

        rAcm = arma::sum(XA.submat(0,0, 1,XA.n_cols-1),1)/NA;
        vAcm = arma::sum(XA.submat(2,0, 3,XA.n_cols-1),1)/NA;
        
        WDString = zeros<mat>(ND,ND);
        
        for (int j = 1; j <= ND-1; j++)
        {
            uvec j11 = find(motionP_result.mP.assign == j);
            uvec j22 = find(motionP_result.mP.assign == j+1);
            int j1 = j11(0);
            int j2 = j22(0);
            if (arma::norm(XD.submat(0,j1,1,j1) - XD.submat(0,j2,1,j2)) <= R_DD_string)
            {
                WDString(j1,j2) = 1;
                WDString(j2,j1) = 1;
            } else {
                WDString(j1,j2) = 0;
                WDString(j2,j1) = 0;
            } 
        }
        
        double ti_2,ti_3, ti_e, ti_g;
        if (flagGather == 1 && flagSeek != 1 && flagEnclose != 1 && flagHerd !=1)
        {
            
            XD_des = motionP_result.dDf.XD_des0;
            XD_des_dot = motionP_result.dDf.XD_des_dot0;
            uD = controlDefender5(XD, SD, regspace(1,1,ND+1), motionP_result.mP.assign, XD_des, XD_des_dot, motionP_result.mP, times(ti), ND);
            for (int j = 0; j < ND; j++)
            {
                if (arma::norm(XD.submat(0,indDef(j)-1, 1,indDef(j)-1)- XD_des.submat(0,j,1,j)) < 1e-3) {
                    defReachCount(j) = 1;
                    if (accu(defReachCount) >= ND)
                    {
                        flagDefReachOpen = 1;
                        flagSeek = 1;
                        flagGather = 0;
                    }
                    
                }
            }
            ti_g = ti;
            // uD.print("uD: ");
        } else if (flagGather!=1 && flagSeek==1 && flagEnclose!=1 && flagHerd!=1)
        {
            // cout << "enter Seek phase when ti = " << ti << endl;
            uvec j11 = find(motionP_result.mP.assign == 1);
            int j1 = j11(0);
            uvec jND_v = find(motionP_result.mP.assign == ND);
            int jND = jND_v(0);
            double phi_ddot;
            // if ( ti == bound-1)
            // {
            //     cout << "" << endl;
            //     cout << "defDesiredOpenForm input: " << endl;
            //     XD.col(ND).print("XD: ");
            //     cout << "RDF_open" << RDF_open << endl;
            //     YA.print("YA:");
            //     cout << "phi: " << motionP_result.dDf.phi << endl;
            //     cout << "phi_dot: " << motionP_result.dDf.phi_dot << endl;
            //     cout << 

            // }
            defDesiredOpenForm(XD.col(ND), RDF_open, YA, motionP_result.dDf.phi, motionP_result.dDf.phi_dot, NA, ND, &XD_des, &XD_des_dot, &phi_ddot, &uDFc_trans, &flagAttInSight);
            // if ( ti == bound-1)
            // {
            //     cout << "" << endl;
            //     cout << "defDesiredOpenForm output: " << endl;
            //     XD_des.print("XD_des: ");
            //     XD_des_dot.print("XD_des_dot:");
            //     cout << "phi_ddot" << phi_ddot << endl;
            //     uDFc_trans.print("uDFc_trans:");
            //     cout << "flagAttInSight: " << flagAttInSight << endl;

            // }
            motionP_result.dDf.phi += dt * motionP_result.dDf.phi_dot;
            if (motionP_result.dDf.phi > 2*M_PI)
            {
                motionP_result.dDf.phi -= 2*M_PI;
            }
            motionP_result.dDf.phi_dot += dt * phi_ddot;
            uD = controlFiniteTimeTrajTracking(XD, indDef, XD_des, XD_des_dot , uDFc_trans, YA, ND,1);
            double ti_2 = ti;

            // if ( ti == bound-1)
            // {
                
            //     cout << "arma::norm(rAcm - rDcm) = "<< arma::norm(rAcm - rDcm) <<endl;
            //     cout << "flagAttInSight:" << flagAttInSight << endl; 
            // }
            
            if (arma::norm(rAcm - rDcm) < 2*rho_sn  && flagAttInSight) {
                flagEnclose = 1;
                flagSeek = 0;
            }
            double ti_s = ti;


        } else if (flagGather != 1 && flagSeek != 1 && flagEnclose == 1 && flagHerd != 1) {
            // cout << "enter enclose formation at ti = " << ti << endl;
            if (flagDForm != 1)
            {
                flagDForm = 1;
            }
            mat rD_des;
            for (int j = 0; j < ND; j++)
            {
                // cout << "j: " << j << endl;
                double RD0 = RDF_closed;
                double thetaD;
                thetaD = motionP_result.dDf.phi + 2*M_PI *(j+1)/ND - M_PI / ND;
                mat tempM1(2,1);
                tempM1(0,0) = cos(thetaD);
                tempM1(1,0) = sin(thetaD);
                // tempM1.print("tempM1: ");
                rD_des.insert_cols(j, rAcm+RD0*tempM1);
                // rD_des.print("rD_des: ");
                mat tempM2;
                tempM2 = join_cols(rD_des.col(j), vAcm);
                
                // tempM2.print("tempM2: ");
                // XD_des.print("XD_des: ");
                XD_des.col(j) = tempM2.as_col();
                // XD_des.print("XD_des: ");
                XD_des_dot.col(j) = zeros<vec>(4);
                // XD_des_dot.print("XD_des_dot: ");
                rSD_goal.resize(rS.n_rows, j+1);
                // rSD_goal.print("rSD_goal: ");
                rSD_goal.col(j) = rS + RD0 * tempM1;
                // rSD_goal.print("rSD_goal: ");

            }
            
            // XD_des.print("XD_des: ");
            // XD_des_dot.print("XD_des_dot: ");
            // rSD_goal.print("rSD_goal: ");
            // Check if the terminal defenders should be connected or not
            // first check if all attackers are inside the convex hull of the
            // defenders or not;
            int flagAttackInHull = 1;
            for (int i = 0; i < NA; i++)
            {
                if(!poly_contain(XD.submat(0,0,1,ND-1), rA.col(i)))
                {
                    flagAttackInHull = 0;
                    break;
                }
            }
            // if (ti == bound-1)
            // {
            //     cout << "flagAttackInHull" << flagAttackInHull << endl;
            // }
            
            uvec j1_v = find(motionP_result.mP.assign==1);
            int j1 = j1_v(0);
            uvec jND_v = find(motionP_result.mP.assign == ND);
            int jND = jND_v(0);
            if ((arma::norm(XD.submat(0,j1,1,j1)-XD_des.submat(0,0,1,0))< bd && arma::norm(XD.submat(0,jND,1,jND)-XD_des.submat(0,ND-1,1,ND-1))<bd) ||  (arma::norm(XD.submat(0,j1,1,j1)-XD.submat(0,jND,1,jND))<=0.95*R_DD_string && flagAttackInHull))
            {
                WDString(j1,jND) = 1;
                WDString(jND,j1) = 1;
            } else 
            {
                WDString(j1,jND) = 0;
                WDString(jND,j1) = 0;
            }
            
            //Check if all the defenders are connected to each other and the
            //StringNet is formed
            int countDefConnect = 0;
            for (int j = 0; j < ND; j++)
            {
                if (j >= ND-1)
                {
                    if (WDString(indDef(ND-1)-1, indDef(0)-1) == 1) 
                    {
                        countDefConnect++;
                    }
                    
                } else
                {
                    if (WDString(indDef(j)-1, indDef(j+1)-1) == 1)
                    {
                        countDefConnect++;
                    }
                }
            }
            if (countDefConnect == ND)
            {

                flagDefConnect = 1;
                flagHerd = 1;
                flagEnclose = 0;
                rDcm = arma::sum(XD.submat(0,0,1,XD.n_cols-1),1) / ND;
                // double thetaD1 = atan2(XD(1,indDef(0)-1) - rDcm(2), XD(0,indDef(0)-1) - rDcm(0));
                // for (int j = 0; j < ND; j++)
                // {
                //     double thetaD = thetaD1+2*M_PI*(j)/(ND);
                // }
                
            }
            
            // if (ti == bound -1)
            // {
            //     cout << "controlDefenderFormation4 inputs: " << endl;
            //     XD.print("XD: ");
            //     indDef.print("indDef: ");
            //     motionP_result.mP.assign.print("assign: ");
            //     XD_des.print("XD_des: ");
            //     XD_des_dot.print("XD_des_dot: ");
            //     uDFc_trans.print("uDFc_trans: ");
            //     YA.print("YA: ");
            // }
            
            uD = controlDefenderFormation4(XD, indDef, motionP_result.mP.assign, XD_des, XD_des_dot, uDFc_trans, YA, NA, ND, 1);
            
            if (ti == bound - 1)
            {
                uD.print("uD: ");
            }
            
            ti_3 = ti;
            ti_e = ti;

            
            XDF_des = XD_des;
            
            X.submat(4*(NA+ND+1)-4,ti,4*(NA+ND+1)-1,ti) = join_cols(rAcm, vAcm);

        } else if (flagDefReachClosed != 1) 
        {
            // cout << "enter loop flagDefReachClosed at ti = " << ti << endl;
            uvec j1_v = find(motionP_result.mP.assign==1);
            int j1 = j1_v(0);
            uvec jND_v = find(motionP_result.mP.assign == ND);
            int jND = jND_v(0);
            if (arma::norm(XD.submat(0,j1,1,j1)-XD.submat(0,jND,1,jND))<=R_DD_string)
            {
                WDString(j1,jND) = 1;
                WDString(jND,j1) = 1;
            } else 
            {
                WDString(j1,jND) = 0;
                WDString(jND,j1) = 0;
            }
            flagDefReachClosed = 1;
            flagHerd = 1;
            for (int j = 0; j < ND; j++)
            {
                if (arma::norm(XD.submat(0,indDef(j)-1,1,indDef(j)-1)-XDF_des.submat(0,j,1,j))<2 && arma::norm(XD.submat(2,j,3,j))<1e-1 && flagDefReachClosed!=1)
                {
                    defReachCount(j) = 1;
                    if(accu(defReachCount) >= ND) {
                        flagDefReachClosed = 1;
                        flagHerd = 1;
                        break;
                    }
                }
            }
            
            uD = zeros<mat>(2,ND+1);
            for (int j = 0; j < ND; j++)
            {
                uD.col(indDef(j)-1) = -0.1*(XD.submat(0,indDef(j)-1,1,indDef(j)-1)-XDF_des.submat(0,j,1,j))-0.1*(XD.submat(2,indDef(j)-1,3,indDef(j)-1));
                double infNorm_uD = max(abs(uD(0,indDef(j)-1)),abs(uD(1,indDef(j)-1)));
                if (infNorm_uD > u_maxD(indDef(j)-1))
                {
                    uD.col(indDef(j)-1) = uD.col(indDef(j)-1)*u_maxD(indDef(j)-1)/infNorm_uD;
                }
            }
            uD.col(ND) = C_d*XD.submat(2,ND,3,ND);
            XD_des = XDF_des;
            ti_3 = ti;
            ti_e = ti;

        }else if (flagGather!=1 && flagSeek!=1 && flagEnclose!=1 && flagHerd==1)
        {
            uvec j1_v = find(motionP_result.mP.assign==1);
            int j1 = j1_v(0);
            uvec jND_v = find(motionP_result.mP.assign == ND);
            int jND = jND_v(0);
            if (arma::norm(XD.submat(0,j1,1,j1)-XD.submat(0,jND,1,jND))<=R_DD_string)
            {
                WDString(j1,jND) = 1;
                WDString(jND,j1) = 1;
            } else 
            {
                WDString(j1,jND) = 0;
                WDString(jND,j1) = 0;
            }

            if (ti - ti_3 < 100)
            {
                XD.col(ND) = join_cols(rAcm,vAcm);
            }
            mat XDFc_des = XD.col(ND);
            double delta_t = dt*(ti-ti_3);
            defDesiredClosedForm(XDFc_des, RDF_closed, motionP_result.dDf.phi, YA, NA, ND, flagNDeven, delta_t, &XD_des, &XD_des_dot, &uDFc_trans);
            uD = controlFiniteTimeTrajTracking(XD,indDef, XD_des, XD_des_dot, uDFc_trans, YA, ND, 0 );
        }
        XD_Des.col(ti) =  XD_des.as_col();
        WDString_mat.slice(ti) = WDString;
        // XD_des.print("XD_Des:");
        // If any attacker has reached inside the protected area stop their
        // control action
        for (int ii = 0; ii < NA; ii++)
        {
            if (arma::norm(rA.col(ii)-rP) < rho_P)
            {
                uA.col(ii) = zeros<mat>(2,1);
            }
        }
        
        U.col(ti) = join_cols(uA.as_col(), uD.as_col());
        // U.col(ti).print("U col(ti)");
        // X.col(ti).print("X col(ti) before modified");
        X.col(ti+1) = modifiedDIDynamics(X.col(ti), U.col(ti), dt, C_d);
        // X.col(ti+1).print("X col(ti+1) after modifed");
        // cout << "11111111111111" << endl;
        XA = reshape(X.submat(0,ti+1,4*NA-1,ti+1),4,NA);
        // YA = XA + arma::mvnrnd(arma::zeros(4,1), Cov_YA, NA);
        YA = XA;
        rA = XA.submat(0,0,1,XA.n_cols-1);
        vA = XA.submat(2,0,3,XA.n_cols-1);
        //  cout << "222222222222222222222" << endl;
        XD = reshape(X.submat(4*NA,ti+1,4*(NA+ND+1)-1,ti+1),4,ND+1);
        mat XDp = reshape(X.submat(4*NA,ti,4*(NA+ND+1)-1,ti),4,ND+1);
        //  cout << "3333333333333333333333" << endl;
        rD = XD.submat(0,0,1,XD.n_cols-1);
        mat vD = XD.submat(2,0,3,XD.n_cols-1);
        
        for (int j = 0; j < NA+ND; j++)
        {
            
            setpoint.pose.position.x = X(4*j, ti+1);
            setpoint.pose.position.y = X(4*j+1, ti+1);
            setpoint.pose.position.z = z_h;
            setpoint_pub_list[j].publish(setpoint);
            // cout<<"message published"<<endl;
            
            body_frame_quad[j].setOrigin(tf::Vector3(setpoint.pose.position.x, setpoint.pose.position.y, setpoint.pose.position.z));
            // cout<<"message broadcasted 1 "<<endl;
            broadcaster.sendTransform(tf::StampedTransform(body_frame_quad[j], ros::Time::now(), "world", "quad" + to_string(j+1)));
            // cout<<"message broadcasted 2"<<endl;

            //std:advance(it,1);
        }
        cout<<"messages broadcasted at ti = "<<ti<<endl;
        ros::spinOnce();
        // cout<<"Ros spin at ti = "<<ti<<endl; 
        rate.sleep();



        // Saturate the velocity if beyond the maximum
        for (int i = 0; i < NA; i++)
        {
            double norm_vA = arma::norm(vA.col(i));
            if(norm_vA > v_maxA(i))
            {
                vA.col(i) = vA.col(i) * v_maxA(i) / norm_vA;
            }
        }
        for (int j = 0; j < ND; j++)
        {
            double norm_vD = arma::norm(vD.col(j));
            if (norm_vD > v_maxD(j))
            {
                vD.col(j) = vD.col(j)*v_maxD(j)/norm_vD;
            }
        }
        XA = join_cols(rA, vA);
        XD = join_cols(rD, vD);

        rDcm = arma::sum(XD.submat(0,0,1,ND-1),1) /ND;
        vDcm = arma::sum(XD.submat(2,0,3,ND-1),1) / ND;
        // if (ti == bound-1)
        // {
        //     cout << "ti: " << ti << endl;
        //     XD.print("XD: ");
        //     XA.print("XA: ");
        //     cout << " " << endl;
        //     rAcm.print("rAcm: ");
        //     rDcm.print("rDcm: "); 
        // }
        
        

        for (int j = 0; j < ND; j++)
        {
            SD(j) = SD_arr(j,ti) + arma::norm(XD.submat(0,j,1,j) - XDp.submat(0,j,1,j));
        }
        SD_arr.col(ti+1) = SD;
        
        arma::mat arr_minRDD;
        mat arr_minRAD(NA,ND);
        mat arr_minRAA(1,NA);
        mat RDjDl(1,ND);
        for (int j = 0; j < ND; j++)
        {
            if (j < ND-1)
            {
                for (int l = j+1; l < ND; l++)
                {
                    RDjDl.col(l) = arma::norm(rD.col(j)- rD.col(l));
                }
                arr_minRDD.resize(1,j+1);
                mat tempM;
                tempM = RDjDl.submat(0,j+1,0,RDjDl.n_cols-1);
                arr_minRDD.col(j) = tempM.min();
                
            }
            for (int i = 0; i < NA; i++)
            {
                arr_minRAD(i,j) = arma::norm(rD.col(j) - rA.col(i));
            }
            
        }

        if(!arr_minRDD.is_empty())
        {
            minRDD(ti+1) = arr_minRDD.min();
        }

        for (int i = 0; i < NA; i++)
        {
            arr_minRAD(i,ND-1) = arma::norm(rD.col(ND-1)- rA.col(i));
        }
        minRAD(ti+1) = arr_minRAD.min();

        //  Find critical distances for the attackers
        mat RAiAl(1,NA);
        for (int i = 0; i < NA; i++)
        {
            if (i < NA-1)
            {
                for (int ii = i+1; ii < NA; ii++)
                {
                    RAiAl = arma::norm(rA.col(i) - rA.col(ii));
                }
                arr_minRAA.resize(1,i+1);
                mat tempM;
                tempM = RDjDl.submat(0,i+1,0,RDjDl.n_cols-1);
                arr_minRDD.col(i) = tempM.min();
            }
        }
        if(!arr_minRAA.is_empty())
        {
            minRAA(ti+1) = arr_minRAA.min();
        }
    
        if (NA > 1)
        {
            minEAO(ti)=(R_m_AA-rhoA_safe)/(control_A_result.R_AO_min);
        }

        minRAAProjS(ti)=(R_m_AD-rhoAD_safe)/control_A_result.R_AAProjS_min;
        
        t=t+dt;
        times(ti+1)=t;

        out_ti = ti;
        // control_A_result.uA.print("uA:");
        // control_A_result.uA0.print("uA0: ");
    }
    auto stop = std::chrono::high_resolution_clock::now();
    
    // U.col(out_ti).print("U col out_ti: ");
    // U.col(out_ti-1).print("U col ti-1: ");
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "out ti : " << out_ti << endl;
    cout << "go through the control loop for " << out_ti+1 << " times " << endl;
    cout << "in total takes: " << duration.count() << " microseconds" << endl;
    cout << "each loop takes: " << duration.count()/(out_ti+1) << " microseconds in average" << endl;
    U.col(out_ti) = U.col(out_ti-1);
    RefTraj.col(out_ti) = RefTraj.col(out_ti-1);
    SigmaProdD_arr.col(out_ti) = sigmaProd;
    minEDO(out_ti) = minEDO(out_ti-1);
    minEAO(out_ti) = minEAO(out_ti-1);
    minRAAProjS(out_ti)=minRAAProjS(out_ti-1);


    // rA.print("rA: ");
    // rD.print("rD: ");
    // delete the unnecessry elements ti<Niter

    times = times.submat(0,0,times.n_rows-1,out_ti + 1);
    X = X.submat(0,0,X.n_rows-1, out_ti+1);
    // X.col(out_ti).print("X col out_ti: ");
    // X.col(out_ti+1).print("col out_ti+1: ");
    U = U.submat(0,0,U.n_rows-1,out_ti+1);
    XD_Des = XD_Des.submat(0,0,XD_Des.n_rows-1, out_ti+1);
    vA_dot_arr = vA_dot_arr.submat(0,0,vA_dot_arr.n_rows-1, out_ti+1);
    FA_dot_arr = FA_dot_arr.submat(0,0,FA_dot_arr.n_rows-1, out_ti+1);
    minRDD = minRDD.submat(0,0,minRDD.n_rows-1, out_ti+1);
    minRAD = minRAD.submat(0,0,minRAD.n_rows-1, out_ti+1);
    minRAA = minRAA.submat(0,0,minRAA.n_rows-1, out_ti+1);
    minEDO = minEDO.submat(0,0,minEDO.n_rows-1, out_ti+1);
    minEAO = minEAO.submat(0,0,minEAO.n_rows-1, out_ti+1);
    SigmaProdD_arr = SigmaProdD_arr.submat(0,0,SigmaProdD_arr.n_rows-1, out_ti);
    // times.print("times");
    // X.save("../../../../../Downloads/swarm_matlab/cppResult/X.csv", csv_ascii);
    // rP.save("../../../../../Downloads/swarm_matlab/cppResult/rP.csv", csv_ascii);
    // rS.save("../../../../../Downloads/swarm_matlab/cppResult/rS.csv", csv_ascii);
    minRDD.save("/home/dasclab/Desktop/minRDD.csv", csv_ascii);
    // std::ofstream myfile;
    // myfile.open ("rho_S.csv");
    // myfile << rho_S;
    // myfile.close();
    // myfile.open("rho_Acon.csv");
    // myfile << rho_Acon;
    // myfile.close();
    // cout << "rho_S: " << rho_S<< endl;
    // cout <<"rho_Acon: " << rho_Acon << endl;
    // cout << "rho_P: " << rho_P << endl;
    // times.save("../../../../../Downloads/swarm_matlab/cppResult/times.csv", csv_ascii);
    // WDString_mat.save("../../../../../Downloads/swarm_matlab/cppResult/WDString_mat.csv");
    
}






int main(int argc, char **argv) {
    ros::init(argc, argv, "single_swarm_stringnet_herding_node");
    ros::NodeHandle nh;
    tf::TransformBroadcaster broadcaster;
    tf::Transform quad_body_frame(tf::Transform::getIdentity());
    std::deque<tf::Transform> body_frame_quad;

    
    //M-Air Dimensions (The local common frame follows ENU convection) with the origin close to the small door near the pavilion 
    float X_max=21, X_min=-1, Y_max=-2.5, Y_min=-36;
    ros::Rate rate(20.0);


    std::vector<ros::Subscriber> state_sb_list;
    std::vector<ros::Subscriber> curr_pos_list;
    std::vector<ros::Subscriber> curr_off_sb_list;
    std::vector<ros::Publisher> setpoint_pub_list;

    int NA;
    int ND;
    double z_h;
    int flag_initial_position; // the flag to set initial position of the code:
                               // 0-the value in cpp code, 1-current positon from gstation topic, 2-manually setting from launch file
                               // setting format: [x1, y1, x2, y2,...xN, yN]

    nh.param<int>("num_of_attackers",NA,1);
    nh.param<int>("num_of_defenders",ND,3);
    nh.param<double>("z_h", z_h, 2.5);
    nh.param<int>("setting_initial_position", flag_initial_position,0);
    std::vector<int> quad_ids(NA+ND);
    std::vector<int> quad_ids0(NA+ND);
    std::vector<double> initial_pos_from_launch(NA+ND);
    std::vector<double> initial_pos_from_launch_default(NA+ND); //which is the same as the code position ==> flag_initial_positon == 0

    for (int i = 0; i < NA+ND; i++)
    {
        body_frame_quad.push_back(tf::Transform::getIdentity());
    }
    initial_pos_from_launch_default = {6.1548,-33.0553,
                                11.9907, 9.3724, 
                                16.1838, -10.7580,
                                -7.1278, -5.8411};
    nh.param("initial_pos", initial_pos_from_launch, initial_pos_from_launch_default);
    std::cout<<"positions read" << std::endl;
    for (int i = 0; i < ND+NA; i++)
    {
        quad_ids0.push_back(i);
    }
    std::cout<<"positions stored" << std::endl;


    nh.param("quad_ids",quad_ids,quad_ids0);
    for (int i = 0; i < ND+NA; i++)
    {
        string state = "mavros/state";
        string quad = "quad";
        string num = to_string(quad_ids[i]);
        string quad_state = "/"+quad+num+"/"+state;
        //cout<<"subscribers creating"<<endl;
        ros::Subscriber state_sb = nh.subscribe<mavros_msgs::State>(quad_state, 10, boost::bind(state_cb,boost::placeholders::_1,i));
        //cout<<"subscribers created"<<endl;

        string gstation = "gstation_position";
        string quad_gstation = "/"+quad+num+"/"+gstation;
        ros::Subscriber curr_pos = nh.subscribe<geometry_msgs::PoseStamped>(quad_gstation,20, boost::bind(current_pos_cb,boost::placeholders::_1, i));

        string offset = "local_ENU_offset";
        string quad_offset = "/"+quad+num+"/"+offset;
        ros::Subscriber curr_off_sb = nh.subscribe<geometry_msgs::PoseStamped>(quad_offset, 10, boost::bind(ENUoff_cb,boost::placeholders::_1, i));

        string setpoint = "desired_setpoint";
        string quad_setpoint = "/"+quad+num+"/"+setpoint;
        ros::Publisher setpoint_pub = nh.advertise<geometry_msgs::PoseStamped>(quad_setpoint, 10);

        state_sb_list.push_back(state_sb);
        curr_pos_list.push_back(curr_pos);
        curr_off_sb_list.push_back(curr_off_sb);
        setpoint_pub_list.push_back(setpoint_pub);
    }

    cout<<"Subscribers  & publishers created"<<endl;
    for (int i = 0; i < quad_ids.size(); i++)
    {
        cout << "quad: " << quad_ids[i] << endl;
    }
    
    calAllparametersExperiment(NA,ND, flag_initial_position, current_pos_list, initial_pos_from_launch);
    AllocateMemory(NA,ND);
    measurements(NA,ND);
    initial_contorl(NA, ND);
    // cout << "complete initial control" << endl;
    getMotionPlan(NA,ND);
    // cout << "complete getMotionPlan" << endl;
    calDistance(ND);
    // cout << "complete calDistance" << endl;
    control_loop(NA,ND,z_h,
                 state_sb_list,
                 curr_pos_list,
                 curr_off_sb_list,
                 setpoint_pub_list,
                 rate,
                 broadcaster,
                 body_frame_quad);
    cout << "R_m_DD" << R_m_DD <<endl;
    ros::shutdown();
}