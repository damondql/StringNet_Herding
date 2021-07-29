#pragma once
#include "AllParametersExperiment.cpp"
#include <armadillo>
#include <math.h>
#include "helperFunction.cpp"
// using namespace std;
using namespace arma;

struct control_attacker_t
{
    arma::mat uA;
    mat uA0;
    double R_AO_min;
    double R_AAProjS_min;
    mat vA_des;
    mat vA_des_dot;
    mat SigmaProdD;
    mat F_A;
    mat F_A_dot;
};



control_attacker_t controlAttacker4(mat XA, mat XA_goal, mat XA_goal_dot,
                      int flagEnclose,int flagHerd,
                      mat XD, mat WA,
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
    
    mat SigmaProdD(NA,1,fill::ones);
    double R_AO_min = INFINITY;
    double R_AAProjS_min = INFINITY;

    mat F_A(2,NA,fill::zeros);
    mat F_A_dot(2,NA,fill::zeros);
    mat F_AD(2,NA,fill::zeros);
    mat F_AD_dot(2,NA,fill::zeros);
    mat uA;
    mat uA0;
    mat vA_dot;
    mat vA_Des;
    mat vA_Des_dot;
    // cout << "before for loop" << endl;
    for (int i = 0; i < NA; i++)
    {
        mat WDString_temp = WDString;
        mat rA = XA.submat(0,i,1,i);
        mat vA = XA.submat(2,i,3,i);

        if (arma::norm(vA) > v_maxA(i))
        {
            vA = vA * v_maxA(i) / arma::norm(vA);
        }
        
        mat rA_goal = XA_goal.submat(0,i,1,i);
        mat vA_goal = XA_goal.submat(2,i,3,i);

        double digmaProd = 1;
        double sigmaBarProd = 1;
        double sigmaSum = 0;
        double sigmaBarSum = 0;
        mat Sigma = zeros<mat>(NO+ND+NA,1);
        mat Sigma_dot = Sigma;

        vec uAOv = zeros<vec>(2);
        vec uAOr = uAOv;

        vec uAFv = zeros<vec>(2);
        vec uAFr = uAFv;

        uvec find_result = arma::find(WA.row(i) == 1);

        //check for nearby defenders
        // cout << "111111111111111" << endl;
        double sigmaSumD = 0;
        double sigmaProdD = 1;
        double minRAD = INFINITY;
        vec uADv(2,fill::zeros);
        vec uADr(2,fill::zeros);

        int countAPS = 0;
        vec uAD_pot(2);
        // cout << "22222222222222222222" << endl;
        if (i > NA - NA_sep)
        {
            // cout << "enter if " << endl;
            double R_m = 15*R_m_AD;
            double R_underbar = R_m+20;
            double R_bar = R_m + 25;
            vec potentialControl_result = potentialControl(XA.col(i),XD.submat(0,0,XD.n_rows-1, ND-1),
                                                           2*rho_c_A,sigma_parameters(R_underbar,R_bar),
                                                           R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);
            uAD_pot = potentialControl_result.subvec(0,1);
            minRAD = potentialControl_result(3);
        } else {
            // cout << "enter else" << endl;
            mat sig_para(4,1);
            sig_para(0,0) = A_A_D;
            sig_para(1,0) = B_A_D;
            sig_para(2,0) = C_A_D;
            sig_para(3,0) = D_A_D;
            // XD.submat(0,0,XD.n_rows-1, ND-1).print("potential input:");
            vec potentialControl_result = potentialControl(XA.col(i),XD.submat(0,0,XD.n_rows-1, ND-1),
                                                           rho_c_A,sig_para,
                                                           R_m_AD,R_bar_AD, R_u_AD, Rij0(0),kADr,kADv, alphaADv);
            // potentialControl_result.print("finish potential Control:");
            uAD_pot = potentialControl_result.subvec(0,1);
            minRAD = potentialControl_result(2);
        }
        // uAD_pot.print("uAD_pot: ");
        // cout << "minRAD: " << minRAD << endl;
        // cout << "3333333333333333333" << endl;
        double R_underbar = R_bar_AD;
        double R_bar = R_u_AD;
        double R_m = R_m_AD;
        vec rA_temp;
        vec vA_temp;
        if(flagHerd == 1)
        {
            rA_temp = rAcm;
            vA_temp = vAcm;
            R_underbar = R_bar_AD+rho_Acon;
            R_bar = R_u_AD + rho_Acon;
            R_m = R_m_AD + rho_Acon;
            double distMin = INFINITY;
            for (int j = 0; j < ND; j++)
            {
                double dist = arma::norm(rAcm - XD.submat(0,j,1,j));
                if (dist < distMin)
                {
                    distMin = dist;
                }
                
            }
            mat rv_Acm(rAcm.n_rows+vAcm.n_rows, rAcm.n_cols);
            rv_Acm.submat(0,0,rAcm.n_rows-1,0) = rAcm;
            rv_Acm.submat(rAcm.n_rows,0,rv_Acm.n_rows-1,0) = vAcm;
            vec potentialControl_result = potentialControl(rv_Acm, XD.submat(0,0,XD.n_rows-1,ND-1),
                                                            rho_c_A,sigma_parameters(R_underbar,R_bar),
                                                            R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);
            uAD_pot = potentialControl_result.subvec(0,1);
        } else {
            rA_temp = rA;
            vA_temp = vA;
        }
        // rA_temp.print("rA_temp: ");
        // vA_temp.print("vA_temp: ");
        // cout << "444444444444444444444" << endl;
        mat rAProjS;
        mat vAProjS;
        vec uAAProj(2,fill::zeros);
        for (int j = 0; j < ND; j++)
        {
            find_result = arma::find(WDString_temp.row(j) == 1);
            if (!find_result.is_empty())
            {
                for (int ii = 0; ii < find_result.n_elem; ii++)
                {
                    int jj = find_result(ii);
                    countAPS++;
                    // cout << "555555555555555555555555" << endl;
                    vec projection_resutl = projectionOnLine(rA_temp, XD.submat(0,jj,1,jj), XD.submat(0,j,1,j));
                    rAProjS.insert_cols(countAPS-1, projection_resutl.subvec(0,1));
                    mat rTDD = XD.submat(0,j,1,j) - XD.submat(0,jj,1,jj);
                    rTDD = rTDD / arma::norm(rTDD);
                    mat tempM;
                    tempM = rTDD.t() * vA_temp;
                    // cout << "666666666666666666666666" << endl;
                    if (tempM(0,0) < 0)
                    {
                        rTDD = -rTDD;
                    }
                    mat vAProjS0, vAProjS1;
                    // rTDD.print("rTDD: ");
                    // vA_temp.print("vA_temp: ");
                    tempM = rTDD.t() * vA_temp;
                    vAProjS0 = tempM(0,0) * rTDD;
                    // cout<< "6.1" << endl;
                    vAProjS1 = XD.submat(2,j,3,j) + 
                               (XD.submat(2,jj,3,jj) - XD.submat(2,jj,3,jj)) 
                               * arma::norm(rAProjS.submat(0,countAPS-1,1,countAPS-1) - XD.submat(0,j,1,j))
                               * arma::norm(XD.submat(0,jj,1,jj) - XD.submat(0,j,1,j));
                    // cout << "777777777777777777777777" << endl;
                    mat cal_result;
                    cal_result = vAProjS1 + vAProjS0;
                    vAProjS.insert_cols(countAPS-1,cal_result);
                    WDString_temp(j,jj) = 0;
                    WDString_temp(jj,j) = 0;
                }
                
            }
            
        }
        
        // cout << "88888888888888888888888" << endl;
        if(flagHerd == 1) {
            if (countAPS > 2)
            {
                vec dist;
                dist.resize(countAPS);
                for (int js = 0; js < countAPS; js++)
                {
                    dist(js) = arma::norm(rAcm - rAProjS.submat(0,js,1,js));
                }
                uvec ind = sort_index(dist);
                mat copy_m;
                copy_m = rAProjS;
                rAProjS.reset();
                rAProjS.insert_cols(rAProjS.n_cols,copy_m.col(ind(0)));
                rAProjS.insert_cols(rAProjS.n_cols,copy_m.col(ind(1)));
                copy_m = vAProjS;
                vAProjS.insert_cols(vAProjS.n_cols,copy_m.col(ind(0)));
                vAProjS.insert_cols(vAProjS.n_cols,copy_m.col(ind(1)));
                countAPS = 2;
            }
        }
        // cout << "999999999999999999999999999" << endl;
        for (int js = 0; js < countAPS-1; js++)
        {
            mat rAvA_temp;
            rAvA_temp = arma::join_cols(rA_temp, vA_temp);
            mat rAvA_ProjS_temp;
            // cout << "0.1" << endl;
            rAvA_ProjS_temp = rAProjS.submat(0,js,1,js);
            rAvA_ProjS_temp = arma::join_cols(rAvA_ProjS_temp, vAProjS.col(js));
            // cout << "0.2" << endl;
            vec potentialControl_result = potentialControl(rAvA_temp, rAvA_ProjS_temp, 
                                                           2*rho_c_A,sigma_parameters(R_underbar,R_bar),
                                                           R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);

        }
        // cout << "000000000000000000000000000000" << endl;
        if (arma::norm(rAcm) < 1.2 * rho_Acon)
        {
            double thetaAAcm = atan2(rA(1)-rAcm(1), rA(0)- rAcm(0));
            vec temp_v = {cos(thetaAAcm), sin(thetaAAcm)};
            mat rAProjC = rAcm+rho_Acon * temp_v;
            vec rTP(2);
            rTP(0) = cos(thetaAAcm + M_PI/2);
            rTP(1) = sin(thetaAAcm + M_PI/2);
            mat cal_result = rTP.t() * vA;
            if (cal_result(0,0) < 0)
            {
                rTP = -rTP;
            }
            mat vAProjC = rTP.t()* vA * rTP + vAcm;
            double Rik0 = Rik00(0);
            double Ri_ik = arma::norm(rA - rAProjC);
            if(Ri_ik - R_m_AA > tol) {
                mat nabla_ri_ViP = kAOr2*(rA-rAProjC)/Ri_ik/abs(Ri_ik-R_m_AA)*(pow((Ri_ik-R_m_AA),2)-pow(Rik0,2))/(pow((Ri_ik-R_m_AA),2)+pow(Rik0,2));
            } else {
                mat nabla_ri_ViP = - kAOr2*(rA-rAProjC)*largeP;
            }
        }
        // cout << "11111111111111111111111111111111" << endl;
        mat duA_goal;
        if(i == 0) 
        {
            duA_goal = XA_goal_dot.submat(2,i,3,i) - (kAPr*(rA-rA_goal)+kAPv*(vA-vA_goal));
        } else
        {
            duA_goal = XA_goal_dot.submat(2,i,3,i) - (10*kAPr*(rA-rA_goal)+10*kAPv*(vA-vA_goal));
        }
        double norm_duA_goal = arma::norm(duA_goal);
        // cout << "norm_duA_goal: " << norm_duA_goal << endl;
        if (norm_duA_goal > 1e-10 && flagHerd != 1)
        {
            mat cal_result;
            cal_result = duA_goal+uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*arma::norm(vA)*vA;
            uA.insert_cols(uA.n_cols,cal_result);
        } else 
        {
            mat cal_result;
            cal_result = uAFv+uAFr+uAOv+uAOr+uAD_pot+uAAProj+C_d*arma::norm(vA)*vA;
            uA.insert_cols(uA.n_cols,cal_result);
        }
        // cout << "22222222222222222222222222222222222" << endl;
        uA0.insert_cols(uA0.n_cols,uA.col(i));

        if(flagHerd != 1 && i ==1)
        {
            R_m=1.5*R_m_AD;
            R_underbar=R_m+2;
            R_bar=R_m+5;  
            mat temp_M;
            vec num;
            num = regspace(2,1,ND);
            temp_M.resize(XD.n_rows, num.n_elem);
            for (int ii = 0; ii < num.n_elem; ii++)
            {
                temp_M.col(ii) = XD.col(num(ii));
            }
            vec potentialControl_result = potentialControl(XA.col(i), temp_M, 
                                                           2*rho_c_A,sigma_parameters(R_underbar,R_bar),
                                                           R_m,R_underbar,R_bar, R_bar+10,kADr,kADv, alphaADv);
            vec uAD_pot2;
            uAD_pot2 = potentialControl_result.subvec(0,1);
            uA.col(i) = uA.col(i) + uAD_pot2;
        }

        double norm_uA = arma::norm(uA.col(i));
        // cout << "333333333333333333333333333333" << endl;
        double uMaxA;
        if (flagEnclose == 1)
        {
            uMaxA = 0.9 * u_maxA(i);
        } else 
        {
            uMaxA = u_maxA(i);
        }
        
        if (norm_uA > uMaxA) {
            uA.col(i) = uA.col(i) * uMaxA / norm_uA;
        }

        vA_dot.insert_cols(i,uA.col(i));
        vA_Des.insert_cols(i,zeros<vec>(2));
        vA_Des_dot.insert_cols(i,zeros<vec>(2));
        SigmaProdD(i,0) = sigmaProdD;
        // cout << "44444444444444444444444444444444" << endl;
        // uA.print("uA: ");
        // uA0.print("uA0:");
        // cout << "R_AO_min: "<< R_AO_min << endl;
        // cout << "R_AAProjS_min" << R_AAProjS_min << endl;
        // vA_Des.print("vA_Des: ");
        // vA_Des_dot.print("vA_Des_dot:");
        // SigmaProdD.print("SigmaProdD");
        // F_A.print("F_A: ");
        // F_A_dot.print("F_A_dot: ");
        
    }
    control_attacker_t control_result;
    control_result.uA = uA;
    control_result.uA0 = uA0;
    control_result.vA_des = vA_Des;
    control_result.vA_des_dot = vA_Des_dot;
    control_result.R_AO_min = R_AO_min;
    control_result.R_AAProjS_min = R_AAProjS_min;
    control_result.F_A = F_A;
    control_result.F_A_dot = F_A_dot;
    return control_result;

}

// int main() {
// //     double R_bar = 10.4177;
// //     double R_underbar = 7.4177;
// //     mat a = sigma_parameters(R_underbar, R_bar);
// //     a.print("a:");

// //     double alphav_12 = 0.5;
// //     double kr_12 = 0.5;
// //     double kv_12 = 0.5;
// //     double R_m_12 = 5.4177;
// //     double rho_sens_1 = 20;
// //     double R_tilde = 20.4176960858139;
// //     mat X1(4,1);
// //     X1(0,0) = 6.15480000000000;
// //     X1(1,0) = -33.0553000000000;
// //     X1(2,0) = 0;
// //     X1(3,0) = 0;
// //     mat X2(4,2);
// //     X2 = {{9.37240000000000,16.1838000000000},
// //          {-7.12780000000000,-5.84110000000000},
// //          {0,0},
// //          {0,0}};
// //     vec re = potentialControl(X1,X2, rho_sens_1, a.as_col(),R_m_12, R_underbar, R_bar,R_tilde, kr_12,kv_12,alphav_12);

// //     vec r = {6.15494381116162,
// //         -33.0547123910252};
// //     vec r1 = {9.37240000000000,
// // -7.12780000000000};
// //     vec r2 = {11.9907000000000,
// // -10.7580000000000};
// //     vec aa = projectionOnLine(r,r1,r2);
// //     aa.print("aa:");
//     calAllparametersExperiment();
//     mat XA;
//     mat XA_goal;
//     mat XA_goal_dot;
//     int flagEnclose;
//     int flagHerd;
//     mat BetaAv0;
//     mat SpeedA0;
//     mat mAv0;
//     mat cAv0;
//     mat XD;
//     mat WA;
//     mat WDString;
//     int ND = 3;
    
//     XA.load("../../../../../Downloads/swarm_matlab/controlA/XA.txt");
//     XA.print("XA:");
//     XA_goal.load("../../../../../Downloads/swarm_matlab/controlA/XA_goal.txt");
//     XA_goal.print("XA_goal:");
//     XA_goal_dot.load("../../../../../Downloads/swarm_matlab/controlA/XA_goal_dot.txt");
//     XA_goal_dot.print("XA_goal_dot:");
//     XD.load("../../../../../Downloads/swarm_matlab/controlA/XD.txt");
//     XD.print("XD:");
//     WA.load("../../../../../Downloads/swarm_matlab/controlA/WA.txt");
//     WA.print("WA:");
//     WDString.load("../../../../../Downloads/swarm_matlab/controlA/WDString.txt");
//     WDString.print("WDString:");
//     control_attacker_t a =  controlAttacker4(XA,XA_goal, XA_goal_dot, 0,0,XD,WA,WDString,NA,ND);
//     a.uA.print("uA:");
//     a.uA0.print("uA0:");

//     return 0;
// }