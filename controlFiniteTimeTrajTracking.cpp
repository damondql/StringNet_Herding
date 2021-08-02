#pragma once
#include "AllParametersExperiment.hpp"
#include <armadillo>

using namespace arma;
using namespace std;

mat controlFiniteTimeTrajTracking(mat XD, mat indDef, mat XD_des, mat XD_des_dot, mat uDFc_trans, mat XA, int ND, int flagAvoidAcon) {
    mat tempM;
    mat uD(2,ND+1, fill::zeros);
    tempM = XD;
    for (int i = 0; i < indDef.n_elem; i++)
    {
        XD.col(i) = tempM.col(indDef(i)-1);
    }
    XD.print("XD after switch:");
    // XD = tempM;
    int Na;
    Na = XA.n_cols;
    int NO = 0;

    mat Rji0;
    Rji0 = 0.6*R_bar_AD*arma::ones<mat>(1,ND);

    mat rDcm,vDcm,rA,vA,rAcm,vAcm;
    rDcm = arma::sum(XD.submat(0,0,1,XD.n_cols-1),1) / ND;
    vDcm = arma::sum(XD.submat(2,0,3,XD.n_cols-1),1) / ND;
    
    rA = XA.submat(0,0,1,XA.n_cols-1);
    vA = XA.submat(2,0,3,XA.n_cols-1);

    rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1) / Na;
    vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1) / Na;
    vec Sigma;
    vec Sigma_dot;
    // cout << "111111111111111" << endl;
    for (int j = 0; j < ND; j++)
    {
        mat rD = XD.submat(0,j,1,j);
        mat rD_des = XD_des.submat(0,j,1,j);
        mat vD = XD.submat(2,j,3,j);
        mat vD_des = vAcm;
        // rD.print("rD");
        // rD_des.print("rD_des");
        // vD.print("vD");
        // vD_des.print("vD_des:");

        mat dr0 = rD-rD_des + 1/(kDFr2*(2-alphaDFv))*(vD-vD_des) * pow(arma::norm(vD-vD_des),(1-alphaDFv));
        // dr0.print("dr0: ");
        mat dr,dv;
        if (arma::norm(dr0) > 1e-5)
        {
            dr = -kDFr2*(dr0)*pow(arma::norm(dr0),(alphaDFr-1));
        } else {
            dr = zeros<mat>(2,1);
        }

        if (arma::norm(vD-vD_des) > 1e-5)
        {
            dv=-kDFr2*(vD-vD_des)*pow(arma::norm(vD-vD_des),(alphaDFv-1));
        } else {
            dv = zeros<mat>(2,1);
        }
        // dr.print("dr: ");
        // dv.print("dv: ");
        mat uDDv(2,1,fill::zeros);
        mat uDDr(2,1,fill::zeros);
        
        int count = 0;
        for (int l = 0; l < ND; l++)
        {
            double sigma;
            double sigma_dot;
            if (l != j)
            {
                count++;
                mat rDjDl = XD.submat(0,l,1,l) - rD;
                double R_DjDl = arma::norm(rDjDl);
                mat R_DjDl_dot = -rDjDl.t() * vD / R_DjDl;
                // R_DjDl_dot.print("R_DjDl_dot");
                if (R_DjDl < R_u_DD)
                {
                    if (R_DjDl < R_bar_DD)
                    {
                        sigma = 1;
                        sigma_dot = 0;
                    } else if (R_DjDl > R_bar_DD && R_DjDl < R_u_DD)
                    {
                        sigma = A_D_D*pow(R_DjDl,3)+B_D_D*pow(R_DjDl,2)+C_D_D*R_DjDl+D_D_D;
                        sigma_dot = (3*A_D_D*pow(R_DjDl,2)+B_D_D*R_DjDl+C_D_D)*R_DjDl_dot(0);
                    }
                    // cout << "sigma " << sigma << endl;
                    // cout << "sigma_dot " << sigma_dot << endl;
                    double Rjj = Rjj0(j);
                    // cout << "Rjj" << Rjj << endl;
                    mat nabla_rj_Vjj= kDDr*(-rDjDl)/R_DjDl/abs(R_DjDl-R_m_DD) * (pow((R_DjDl-R_m_DD),2)-pow(Rjj,2))/(pow((R_DjDl-R_m_DD),2)+pow(Rjj,2));
                    // cout << "kDDr: " << kDDr <<endl;
                    // cout << "R_m_DD: " << R_m_DD << endl;
                    // cout << "abs(R_DjDl-R_m_DD)" <<abs(R_DjDl-R_m_DD) << endl;
                    // nabla_rj_Vjj.print("nabla_rj_Vjj: ");
                    dv = kDDv*(vD-XD.submat(2,l,3,l)) * pow(arma::norm(vD-XD.submat(2,l,3,l)),(alphaDDv-1));
                    uDDv = uDDv-sigma*dv;
                    uDDr=uDDr-sigma*nabla_rj_Vjj;
                    // uDDr.print("uDDr: ");
                    // uDDv.print("uDDv: ");
                } else {
                    sigma = 0;
                    sigma_dot = 0;
                }
                
                Sigma.resize(count);
                Sigma_dot.resize(count);
                Sigma(count-1) = sigma;
                Sigma_dot(count-1) = sigma_dot;
            }
            
        }
        // Sigma.print("Sigma: ");
        // Sigma_dot.print("Sigma_dot:");

        mat nabla_rj_VjP(2,1,fill::zeros);
        mat dvP(2,1,fill::zeros);
        // cout << "22222222222222222222222" <<endl;
        if (flagAvoidAcon)
        {
            double  thetaDAcm=atan2(rD(1)-rAcm(1),rD(0)-rAcm(0));
            tempM.resize(2,1);
            tempM(0,0) = cos(thetaDAcm);
            tempM(1,0) = sin(thetaDAcm);
            mat rDProj = rAcm+rho_Acon*tempM;
            // rDProj.print("rDProj: ");
            mat rTP(2,1);
            rTP(0,0) = cos(thetaDAcm+M_PI/2);
            rTP(1,0) = sin(thetaDAcm+M_PI/2);
            tempM.reset();
            tempM = rTP.t() * vD;
            if (tempM(0,0) < 0 )
            {
                rTP = -rTP;
            }
            // rTP.print("rTP: ");
            tempM.reset();
            tempM = rTP.t()*vD;
            mat vDProj = tempM(0,0)*rTP + vAcm;
            // vDProj.print("vDProj: ");
            double Rjk0 = Rjk00(0);
            double Rj_jk = arma::norm(rD - rDProj);
            if (Rj_jk - R_m_DD > 1e-1)
            {
                nabla_rj_VjP = kDOr*(rD-rDProj)/Rj_jk/abs(Rj_jk-R_m_DD)*(pow((Rj_jk-R_m_DD),2)-pow(Rjk0,2))/(pow((Rj_jk-R_m_DD),2)+pow(Rjk0,2));
            } else {
                nabla_rj_VjP= -kDOr*(rD-rDProj)*largeP;
            }
            if (arma::norm(vD-vDProj) != 0) {
                dvP = -kDOv2*(vD-vDProj)*pow(arma::norm(vD-vDProj),(alphaDOv-1));
            } else {
                dvP = zeros<mat>(2,1);
            }
            // nabla_rj_VjP.print("nabla_rj_Vjp: ");
            // dvP.print("dvP:");
        }
        
        mat uDOv(2,1,fill::zeros);
        mat uDOr(2,1,fill::zeros);
        // cout << "33333333333333333333333" << endl;
        mat uD1,uD2;
        double norm_uD1, norm_uD2;
        uD1 = dr+uDDv+uDDr+dvP-nabla_rj_VjP+uDOv+uDOr;
        norm_uD1=arma::norm(uD1);
        uD2=dv+C_d*XD.submat(2,j,3,j)*arma::norm(XD.submat(2,j,3,j));
        norm_uD2=arma::norm(uD2);
        cout << "norm_uD1" << norm_uD1 << endl;
        cout << "norm_uD2" << norm_uD2 << endl;
        mat uD1_hat, uD2_hat;
        if (norm_uD1 > 1e-15)
        {
            uD1_hat = uD1/norm_uD1;
        } else {
            uD1_hat = zeros<mat>(2,1);
        }

        if (norm_uD2 > 1e-15)
        {
            uD2_hat = uD2/norm_uD2;
        } else {
            uD2_hat = zeros<mat>(2,1);
        }
        uD1_hat.print("uD1_hat: ");
        uD2_hat.print("uD2_hat: ");
        cout << "std::min(umd1,norm_uD1)" << std::min(umd1,norm_uD1) << endl;
        cout << "std::min(umd2,norm_uD2)" << std::min(umd2,norm_uD2) << endl;
        tempM.reset();
        tempM = XD_des_dot.submat(2,j,3,j)+std::min(umd1,norm_uD1)*uD1_hat+std::min(umd2,norm_uD2)*uD2_hat;;
        tempM.print("tempM:");
        uD.print("uD before cal:");
        uD.col(j) = XD_des_dot.submat(2,j,3,j)+std::min(umd1,norm_uD1)*uD1_hat+std::min(umd2,norm_uD2)*uD2_hat;
        uD.col(j).print("uD col j:");
        uD.print("uD after cal: ");
        

    }
    uD.col(ND) = uDFc_trans;
    uD.print("uD after all cal before switch:");
    tempM.reset();
    tempM = uD;
    for (int i = 0; i < uD.n_cols; i++)
    {
        uD.col(indDef(i)-1) = tempM.col(i);
    }

    


    return uD;
}