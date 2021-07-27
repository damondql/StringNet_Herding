#pragma once
#include "AllParametersExperiment.hpp"
#include <armadillo>
#include <math.h>
#include "helperFunction.cpp"
using namespace std;
using namespace arma;

struct control_attacker_t
{
    mat uA;
    mat uA0;
    double R_AO_min;
    double R_AAProjS_min;
    mat vA_des;
    mat vA_des_dot;
    mat FlagAttInObs;
    mat BetaAv0;
    mat SpeedA0;
    mat mAv0;
    mat cAv0;
    double SigmaProdD;
    mat FA;
    mat FA_dot;
};



void controlAttacker4(mat XA, mat XA_goal, mat XA_goal_dot,
                      mat A_lead,mat FlagAttInObs,mat flagEnclose,int flagHerd,
                      mat BetaAv0, mat SpeedA0,mat mAv0,mat cAv0, mat XD, mat WA,
                      mat WDString,int Na, int ND){
    double tol = 0.5;
    int NO = 0;

    mat rD, vD;
    if (ND > 0)
    {
        rD.resize(2,ND);
        vD.resize(2,ND);
        rD = XD.submat(0,0,1,ND-1);
        vD = XD.submat(2,0,3,ND-1);
    }
    
    arma::colvec rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1)/Na;
    arma::colvec vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1)/Na;
    
}

int main() {
//     double R_bar = 10.4177;
//     double R_underbar = 7.4177;
//     mat a = sigma_parameters(R_underbar, R_bar);
//     a.print("a:");

//     double alphav_12 = 0.5;
//     double kr_12 = 0.5;
//     double kv_12 = 0.5;
//     double R_m_12 = 5.4177;
//     double rho_sens_1 = 20;
//     double R_tilde = 20.4176960858139;
//     mat X1(4,1);
//     X1(0,0) = 6.15480000000000;
//     X1(1,0) = -33.0553000000000;
//     X1(2,0) = 0;
//     X1(3,0) = 0;
//     mat X2(4,2);
//     X2 = {{9.37240000000000,16.1838000000000},
//          {-7.12780000000000,-5.84110000000000},
//          {0,0},
//          {0,0}};
//     vec re = potentialControl(X1,X2, rho_sens_1, a.as_col(),R_m_12, R_underbar, R_bar,R_tilde, kr_12,kv_12,alphav_12);

//     vec r = {6.15494381116162,
//         -33.0547123910252};
//     vec r1 = {9.37240000000000,
// -7.12780000000000};
//     vec r2 = {11.9907000000000,
// -10.7580000000000};
//     vec aa = projectionOnLine(r,r1,r2);
//     aa.print("aa:");
    return 0;
}