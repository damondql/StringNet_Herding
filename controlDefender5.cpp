#pragma once
#include "AllParameters.hpp"
#include <armadillo>
#include "defInitDesirePos.cpp"

using namespace std;
using namespace arma;

mat controlDefender5(mat XD, mat SD, vec indD, vec assign, mat XD_des, mat XD_des_dot, vec indDes, mat uD, motionPlan motionPlan_result, double t, int ND){
    // mat temp_M;
    // for (int i = 0; i < indDef.n_elem; i++)
    // {
    //     temp_M.insert_cols(i,XD.col(indDef(i)-1));
    // }
    // XD = temp_M;
    // temp_M.reset();
    // for (int i = 0; i < motionPlan_result.assign.n_elem; i++)
    // {
    //     temp_M.insert_cols(i, XD_des.col(motionPlan_result.assign(i)-1));
    // }
    // XD_des = temp_M;
    mat Rji0;
    Rji0 = 0.6 * R_bar_AD * ones<mat>(1,ND);

    mat rDcm,vDcm;
    rDcm = arma::sum(XD.submat(0,0,1,ND-1),1)/ND;
    vDcm = arma::sum(XD.submat(2,0,3,ND-1),1)/ND;
    
    mat rD,vD;
    rD = XD.submat(0,0,1,XD.n_cols-1);
    vD = XD.submat(2,0,3,XD.n_cols-1);

    uD= zeros<mat>(2,ND);
    for (int jj = 0; jj < ND; jj++)
    {
        int j = indD(jj);
        // cout << "enter for loop" << endl;
        // cout << "j: " << j << endl;
        if (t >= motionPlan_result.startTime(j))
        {
            // cout << "t is greater than startTime" << endl;
            double s_ddot = 0;
            double s_dot = arma::norm(vD.col(j));
            vec S = motionPlan_result.path[j].S;
            double P = motionPlan_result.path[j].P;
            mat rV = motionPlan_result.path[j].rV;
            mat rVC = motionPlan_result.path[j].rVc;
            int segType = motionPlan_result.path[j].segType;
            double v_bar = motionPlan_result.pathVel[j].v_bar;
            double s_bar1 = motionPlan_result.pathVel[j].s_bar1;
            double s_bar2 = motionPlan_result.pathVel[j].s_bar2;
            uvec ind0 = find(S <= SD(j));
            int indS = ind0(ind0.n_elem-1);
            // ind0.print("ind0");
            // cout << "indS: " << indS << endl;
            if (indS <= 0) //if indS <= length(P) which is 1 since obs free
            {
                // cout << "cal ud in if loop" << endl;
                double S_prev = S(indS);
                if (indS % 2 == 0)
                {
                    if (SD(j) < s_bar1 + S_prev)
                    {
                        s_ddot=u_maxD(j);
                    } else if (SD(j) > s_bar2+S_prev) 
                    {
                        s_ddot=-(u_maxD(j));
                    } else 
                    {
                        s_ddot=u_maxD(j);
                    }
                    // cout << "s_ddot" << s_ddot << endl;
                    // cout << "calculating ud in if loop" << endl;
                    uD.col(j) = s_ddot*(rV.col(indS+1)-rV.col(indS))/norm(rV.col(indS+1)-rV.col(indS));
                }


            } else {
                // cout << "cal ud in else loop"<< endl;
                // rD.print("rD: ");
                // XD_des.print("XD_des: ");
                // vD.print("vD: ");
                uD.col(j) = -(rD.submat(0,j,1,j) - XD_des.submat(0,jj,1,jj)) - vD.submat(0,j,1,j);
            }
            

        } else {
            // cout << "t is less than startTime, enter else loop" << endl;
            uD.col(j) = zeros<vec>(2);
        }
        
    }
    return uD;
    
}
