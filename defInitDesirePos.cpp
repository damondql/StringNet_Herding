#pragma once
#include <iostream>
#include <armadillo>
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


// struct defDesForm
// {
//     mat rDFc0;
//     mat XD_des0;
//     mat XD_des_dot0;
//     double phi;
//     double phi_dot;
// };

// struct DesiredPos {
//     defDesForm dDf;
//     motionPlan mP;
// };

void defInitDesiredPos(mat XA0, mat XD, int Na, int ND, double RD0, double v_maxA, double DeltaT_s) {
    cout << "enter function  " << endl;
    colvec rAcm = arma::sum(XA0.submat(0,0,1,XA0.n_cols-1),1) / Na;
    colvec vAcm = arma::sum(XA0.submat(2,0,3,XA0.n_cols-1),1) / Na;
    double thetaAcm0=atan2(rAcm(1)-rP(1),rAcm(0)-rP(0))+20*M_PI/180;
    if (thetaAcm0 < 0) {
        thetaAcm0 = thetaAcm0 + 2 * M_PI;
    }
    cout <<"111111111111" <<endl;
    path_elem path_Acm = findShortestPath(rAcm, rP);
    double Qa = path_Acm.S(1);
    cout << "Qa:" << Qa << endl;
    double q1 = 0;
    double q2 = Qa - rho_P;

    double f = 10000;
    int i = 0;
    double tComp = 0;
    double df = 10000;
    mat rD_des(2,ND);
    mat XD_des0(4,ND,fill::zeros);
    mat XD_des_dot0(4,ND,fill::zeros);
    while (fabs(f) > 3 && fabs(df) > 1)
    {
        i = i + 1;
        double q = (q1+q2)/2;
        CoorOnPath Acm0 = findCoordOnPath(q,path_Acm);
        double thetaD = Acm0.thetap(0) - M_PI / 2;
        for (int j = 0; j < ND; j++) {
            vec tempv(2);
            tempv(0) = cos(thetaD);
            tempv(1) = sin(thetaD);
            rD_des.submat(0,j,1,j) = Acm0.rp + 0.5 * (ND-1) * RD0 * tempv - RD0 * tempv*(j);
            XD_des0.submat(0,j,1,j) = rD_des.submat(0,j,1,j);
        }
    }
    rD_des.print("rD_des:");
    XD_des0.print("XD_des0:");
    
    
}

int main(){
    arma::vec rAcm = {5.9807,-33.2064};
    arma::vec rP = {14,-1};
    path_elem p = findShortestPath(rAcm, rP);
    p.rV.print("rV: ");
    p.S.print("S: ");
    double q = 16.0949;
    CoorOnPath A = findCoordOnPath(q, p);
    A.rp.print("rDFc0: ");
    A.thetap.print("thetaAcom0: ");
}