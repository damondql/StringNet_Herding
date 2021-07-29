#pragma once
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

mat sigma_parameters(double R_underbar, double R_bar){
    mat sigma_params(1,4);
    double dR_cube=pow((R_bar-R_underbar),3);
    sigma_params(0,0)=2/dR_cube;
    sigma_params(0,1)=-3*(R_bar+R_underbar)/dR_cube;
    sigma_params(0,2)=6*R_bar*R_underbar/dR_cube;
    sigma_params(0,3)=pow(R_bar,2)*(R_bar-3*R_underbar)/dR_cube;
    return sigma_params;
}


//Control for Collision avoidance of agent 1 from agents 2 using
//potential function approach
//return a 3x1 vector where the first 2 elements are u_pot_12
// and the 3rd element is minDist_12
vec potentialControl(mat X1,mat X2, double rho_sens_1, mat sigma_params_12,
                     double R_m_12, double R_underbar_12, double R_bar_12, 
                     double R_tilde_12, double kr_12, double kv_12, double alphav_12){
    // cout << "enter potentialControl" << endl;
    mat r2(2,X2.n_cols);
    r2 = X2.submat(0,0,1,X2.n_cols-1);
    mat v2(2,X2.n_cols);
    v2 = X2.submat(2,0,3,X2.n_cols-1);
    mat r1(2,X1.n_cols);
    r1 = X1.submat(0,0,1,X1.n_cols-1);
    mat v1(2,X2.n_cols);
    v1 = X1.submat(2,0,3,X1.n_cols-1);
    // v1.print("in potentialControl v1:");
    vec uv_12(2,fill::zeros);
    vec ur_12(2,fill::zeros);
    double minDist_12 = INFINITY;
    double tol = 0.1;
    double largeNum = -1/tol*(pow((tol),2)-pow(R_tilde_12,2))/(pow((tol),2)+pow(R_tilde_12,2));
    double sigma = 0;
    for (int j = 0; j < r2.n_cols; j++)
    {
        // cout << "00000000000000000" << endl;
        vec r12j = r2.col(j) - r1;
        double R_12j = arma::norm(r12j);
        if(R_12j <= rho_sens_1) {
            if(R_12j < R_bar_12) {
                // cout << "111111111111111111111111111" << endl;
                if(R_12j < R_bar_12) {
                    sigma = 1;
                } else if (R_12j > R_underbar_12 && R_12j<R_bar_12) {
                    sigma = sigma_params_12(0)*pow(R_12j,3)+sigma_params_12(1)*pow(R_12j,2)+sigma_params_12(2)*R_12j+sigma_params_12(3);
                }
                vec nabla_r1_V12;
                // cout << "22222222222222222222" << endl;
                if( (R_12j - R_m_12) > tol){
                    nabla_r1_V12 = kr_12 * (r1-r2.col(j)) / R_12j/abs(R_12j-R_m_12) * (pow((R_12j-R_m_12),2)-pow(R_tilde_12,2))/(pow((R_12j-R_m_12),2)+pow(R_tilde_12,2));
                } else {
                    nabla_r1_V12 = -kr_12*(r1-r2.col(j)) / R_12j*largeNum;
                }
                double norm_dv12 = arma::norm(v1 - v2.col(j));
                mat a = v1 - v2.col(j);
                a = a.t();
                mat b  = r1 - r2.col(j);
                mat c = a * b;
                vec dv;
                // cout << "33333333333333333333333333" << endl;
                if (norm_dv12 > 1e-16 && c(0,0) < 0) {
                    dv = kv_12 * (v1 - v2.col(j)) * pow(norm_dv12,alphav_12-1);
                } else {
                    dv = zeros<vec>(2);
                }
                uv_12 = uv_12 - sigma * dv;
                ur_12 = ur_12 - sigma * nabla_r1_V12;
            }
        }
        if (minDist_12 > R_12j)
        {
            minDist_12 = R_12j;
        }
        
    }

    vec result(3);
    result.subvec(0,1) = uv_12 + ur_12;
    result(2) = minDist_12;
    
    return result;
}

vec projectionOnLine(vec r, vec r1, vec r2) {
    double dx = r2(0) - r1(0);
    double dy = r2(1) - r1(1);
    double mVV = dy/dx;
    double cVV = r1(1) - mVV*r1(0); 
    vec rProj(2);
    rProj(0) = (mVV*r(1)+r(0)-mVV*cVV)/(1+pow(mVV,2));
    rProj(1) = mVV*rProj(0)+cVV;
    double lambdaP;
    if (mVV < 1e16)
    {
        lambdaP = (rProj(0)-r1(0))/dx;
    } else {
        lambdaP=(rProj(1)-r1(1))/dy;
    }
    if (lambdaP<0){
        rProj=r1;
    } else if (lambdaP>1) {
        rProj=r2;
    }
    vec result(3);
    result.subvec(0,1) = rProj;
    result(2) = lambdaP;
    return result;
}