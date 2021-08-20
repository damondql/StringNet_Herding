#pragma once
#include <iostream>
#include <armadillo>
#include <fstream>

#include "findCoordOnPath.cpp"
#include "motionPlanForDefOpenForm.cpp"
using namespace std;
using namespace arma;
// struct tanG_prime_elem
// {
//     mat G;
//     mat G_pathType;
//     cube rVO;
//     cube rVO2;
//     mat gTO_all;
//     mat rTO_all;
//     mat obsld_all;
//     mat obsVertId_all;
//     mat obsVertPos_all;
// };

// struct interSec_elem
// {
//     int Flag;
//     int Pbar1;
//     int Pbar2;
//     int flag0;
// };

// struct pathVel_elem
// {
//     mat v;
//     double v_bar;
//     double s_bar1;
//     double s_bar2;
//     mat T;
//     double T_bar1;
//     double T_bar2;
//     double u_maxD;
//     double v_maxD;
//     double v_maxDC;
// };

// struct motionPlan
// {
//     mat assign;
//     double assignCost;
//     double maxOptT;
//     std::vector<path_elem> path;
//     std::vector<tanG_prime_elem> tanG_prime;
//     std::vector<::vector<interSec_elem>> interSec;
//     std::vector<pathVel_elem> pathVel;
//     mat optT;
//     mat leadTime;
//     mat startTime;
//     double tCompMIQP;
//     double tCompMILP;
//     double tComp;
//     mat XD;
//     mat XD_des;
// };


struct defDesForm
{
    mat rDFc0;
    mat XD_des0;
    mat XD_des_dot0;
    double phi;
    double phi_dot;
};

struct DesiredPos {
    defDesForm dDf;
    motionPlan mP;
    int flagGatheringPossible;
};

DesiredPos defInitDesiredPos(mat XA0, mat XD, int NA, int ND, double RD0, double rhoA_con, double v_maxA, double DeltaT_s) {
    DesiredPos result;
    // cout << "enter function  " << endl;
    colvec rAcm = arma::sum(XA0.submat(0,0,1,XA0.n_cols-1),1) / NA;
    colvec vAcm = arma::sum(XA0.submat(2,0,3,XA0.n_cols-1),1) / NA;
    double thetaAcm0=atan2(rAcm(1)-rP(1),rAcm(0)-rP(0))+20*M_PI/180;
    if (thetaAcm0 < 0) {
        thetaAcm0 = thetaAcm0 + 2 * M_PI;
    }
    // cout <<"111111111111" <<endl;
    rAcm.print("rAcm: ");
    rP.print("rP: ");
    path_elem path_Acm = findShortestPath(rAcm, rP);
    path_Acm.S.print("path_Acm, S: ");
    double Qa = path_Acm.S(1);
    // cout << "Qa:" << Qa << endl;
    double q1 = 0;
    double q2 = Qa - rho_P;
    cout << "rho_P: " << rho_P << endl;
    cout <<  "q2: " << q2 << endl;
    double f = 10;
    int i = 0;
    double tComp = 0;
    double df = 10000;
    mat rD_des(2,ND);
    mat XD_des0(4,ND,fill::zeros);
    mat XD_des_dot0(4,ND,fill::zeros);
    while (fabs(f) > 3 && i<500)
    {
        cout << "enter while loop" << endl;
        i = i + 1;
        double q = (q1+q2)/2;
        vec q_v(1);
        q_v(0) = q;
        CoorOnPath Acm0 = findCoordOnPath(q_v,path_Acm);
        // Acm0.rp.print("rDFc0");
        // Acm0.thetap.print("thetaAcm0");
        // double thetaD = Acm0.thetap(0) - M_PI / 2;
        double theta1 = Acm0.thetap(0) - M_PI/2;
        double theta2 = Acm0.thetap(0) - M_PI/2 + M_PI;
        double dR = 2*RD0/(ND-1);
        for (int j = 0; j < ND; j++) {
            vec tempV1(2);
            tempV1(0) = cos(theta1);
            tempV1(1) = sin(theta1);
            vec tempV2(2);
            tempV2(0) = cos(theta2);
            tempV2(1) = sin(theta2);
            XD_des0.submat(0,j,1,j) = Acm0.rp + RD0 * tempV1 + j * dR * tempV2;
            XD_des0.submat(2,j,3,j) = zeros<mat>(2,1);
            XD_des_dot0.col(j) = zeros<mat>(4,1);
        }
        // Acm0.rp.print("rDFc0:");
        // rD_des.print("rD_des:");
        // XD_des0.print("XD_des0:");
        result.mP = motionPlanForDefOpenForm(XD, XD_des0, ND);
        // result.mP.assign.print("assign: ");
        // cout << "assignCost: " <<result.mP.assignCost << endl;
        // cout << "maxOptT: " << result.mP.maxOptT << endl;
        tComp += result.mP.tComp;
        // df = f - (result.mP.maxOptT - q/v_maxA + DeltaT_s);
        f = result.mP.maxOptT - (q-rhoA_con)/v_maxA +DeltaT_s;
        if (f<0)
        {
            q2 = q;
        } else {
            q1 = q;
        }
        result.dDf.rDFc0 = Acm0.rp;
        result.dDf.phi = Acm0.thetap(0);
        result.dDf.phi_dot = 0;
    }
    result.mP.tComp = tComp;
    result.dDf.XD_des0 = XD_des0;
    result.dDf.XD_des_dot0 = XD_des_dot0;
    
    if(i >= 500)
    {
        result.flagGatheringPossible = 0;
    } else
    {
        result.flagGatheringPossible = 1;
    }


    return result;
    
}

// int main(){
//     // arma::vec rAcm = {5.9807,-33.2064};
    
//     rP = {0,200};
//     // path_elem p = findShortestPath(rAcm, rP);
//     // p.rV.print("rV: ");
//     // p.S.print("S: ");
//     // double q = 16.0949;
//     // CoorOnPath A = findCoordOnPath(q, p);
//     // A.rp.print("rDFc0: ");
//     // A.thetap.print("thetaAcom0: ");
//     double DeltaT_S = 50;
//     int ND = 18;
//     int NA = ND;
//     calAllparametersExperiment(NA,ND);
//     cout << "rho_safe: " << rho_safe << endl; 
//     // double RD0 = 4.4;
//     // rho_Acon = 3.78533660770541;
//     // rho_P = 1;
//     // rho_safe = 2.01179739054263;
//     // rho_sn = 5.79713399824804;
//     // v_maxA[0] = 0.22;
//     // u_maxD = {0.0982320000000000, 0.0982320000000000, 0.0982320000000000};
//     // v_maxD = {0.495560288965934,0.49556028896593,0.495560288965934};
//     // v_maxDC = {0.392380551560576,0.392380551560576,0.392380551560576};
//     XD.load("../../../../../Downloads/multi_swarm/defV/XD.txt");
//     XA0.load("../../../../../Downloads/multi_swarm/defV/XA0.txt");
//     XA0.print("XA0: ");
//     XD.print("XD: ");

//     std::ifstream fin("../../../../../Downloads/multi_swarm/defV/RD0.txt");
//     double RD0;
//     fin >> RD0;
//     fin.close();

//     fin.open("../../../../../Downloads/multi_swarm/defV/rhoA_con.txt");
//     double rhoA_con;
//     fin >> rhoA_con;
//     fin.close();

//     fin.open("../../../../../Downloads/multi_swarm/defV/v_maxA.txt");
//     double v_maxA_;
//     fin >> v_maxA_;
//     fin.close();

//     DesiredPos a = defInitDesiredPos(XA0, XD, NA, ND, RD0, rhoA_con , v_maxA_, DeltaT_S);
//     a.mP.assign.print("motion plan assign: ");
//     cout << "motion plan assign cost: " << a.mP.assignCost << endl;
//     cout << "motion plan maxOptT: " << a.mP.maxOptT << endl;
//     cout << "defDesForm phi " <<a.dDf.phi << endl;
//     cout << "tcomp: " << a.mP.tComp << endl;
// }