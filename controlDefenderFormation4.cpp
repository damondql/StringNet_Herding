#pragma once
#include "AllParametersExperiment.hpp"
#include "helperFunction.cpp"
#include <armadillo>

using namespace arma;

mat controlDefenderFormation4(mat XD, mat indDef, mat assign, mat XD_des, mat XD_des_dot, mat uDFc_trans, mat XA, int Na, int ND, int flagAvoidAcon){
    int na = assign.n_elem;
    int nid = indDef.n_elem;
    mat uD(2,ND+1,fill::zeros);
    mat tempM;
    tempM = XD;
    for (int i = 0; i < indDef.n_elem; i++)
    {
        XD.col(i) = tempM.col(indDef(i)-1);
    }
    mat rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1)/NA;
    mat vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1)/NA;
    for (int j = 0; j < ND; j++)
    {
        mat rD = XD.submat(0,j,1,j);
        mat vD = XD.submat(2,j,3,j);

        mat uAcon(2,1,fill::zeros);
        if (flagAvoidAcon)
        {
            double thetaDAcm = atan2(rD(1)-rAcm(1),rD(0)-rAcm(0));
            tempM.resize(2,1);
            tempM(0,0) = cos(thetaDAcm);
            tempM(1,0) = sin(thetaDAcm);
            mat rDProj = rAcm+rho_Acon*tempM;
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
            double R_underbar = 4;
            double R_bar = 3;
            double R_m = 1.5;
            mat rvDproj;
            rvDproj = arma::join_cols(rDProj, vDProj);
            // rvDproj.print("rvDproj");
            vec potentialControl_result = potentialControl(0.1,XD.col(j),rvDproj,20,sigma_parameters(R_underbar,R_bar),R_m,R_underbar,R_bar, R_bar+10,kDOr,kDOv2, alphaDOv);
            uAcon = potentialControl_result.subvec(0,1);
        }

        mat uDOv(2,1,fill::zeros);
        mat uDOr(2,1,fill::zeros);
        mat uD1,uD2;
        double norm_uD1, norm_uD2;
        uD1=-kDDesr*(rD-XD_des.submat(0,j,1,j))-kDDesv*(vD-XD_des.submat(2,j,3,j)) + uAcon+uDOv+uDOr;
        norm_uD1=arma::norm(uD1);
        // cout << "norm_uD1: " << norm_uD1 << endl;
        uD2=-kDDesv*(vD-XD_des.submat(2,j,3,j))+C_d*vD*arma::norm(vD);
        norm_uD2=arma::norm(uD2);
        // cout << "norm_uD2: " << norm_uD2 << endl;
        uD.col(j)=XD_des_dot.submat(2,j,3,j) + std::min(umd_e1,norm_uD1)*uD1/norm_uD1 + std::min(umd_e2,norm_uD2)*uD2/norm_uD2;

    }
    uD.col(ND) = uDFc_trans;
    // uD.print("uD after all cal before switch:");
    tempM.reset();
    tempM = uD;
    for (int i = 0; i < uD.n_cols; i++)
    {
        uD.col(indDef(i)-1) = tempM.col(i);
    }

    return uD;
}