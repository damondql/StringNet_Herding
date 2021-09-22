#pragma once
#include "AllParameters.cpp"
#include <armadillo>
using namespace arma;
// using namespace std;


void defDesiredOpenForm(field<mat> XDFc, int clusterDNum, mat *XD_des, mat *XD_des_dot, vec RDF_OPEN, mat rho_SN, std::vector<vec> indDes, int NDinCluster, int NClusterD, mat XA,
                        mat rAcm, mat vAcm, mat rhoA_con, int NClusterAD, int clusterADNum, double phi, double phi_dot, int NA, int ND, int flagAvoidObs){
                        // mat *XD_des, mat *XD_des_dot, double *phi_ddot,
                        // mat *uDFc_trans, int *flagAttInSight) {
                        // ,mat XA, double phi, double phi_dot,int NA,int ND, mat *XD_des, mat *XD_des_dot, double *phi_ddot, mat *uDFc_trans, int *flagAttInSight
    double chi = 0.55;
    double RDF0 = RDF_OPEN(clusterDNum);
    double phi_ddot_max=(0.9-chi)*u_maxD(1)/RDF0;

    mat rAcm,vAcm;
    rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1)/NA;
    vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1)/NA;

    mat rDFc,vDFc;
    rDFc = XDFc(clusterDNum).submat(0,0,1,0);
    vDFc = XDFc(clusterDNum).submat(2,0,3,0);

    mat drc,dvc;
    drc = rAcm.col(clusterADNum) - rDFc;
    dvc = vAcm.col(clusterADNum) - vDFc;

    double phi_des = atan2(-drc(1), -drc(0));
    if(phi_des < 0) {
        phi_des += 2*M_PI;
    }
    double phi_des_dot=(-drc(0)*dvc(1)+drc(1)*dvc(0))/pow(arma::norm(drc),2);
    double dphi = (phi - phi_des);

    if (dphi < -M_PI) {
        dphi += 2*M_PI;
    }
    int flagAttInSight1 = 0;
    // drc.print("drc: ");
    // cout << "arma::norm(drc): " << arma::norm(drc) << endl;
    // cout << "asin(RDF0): " << asin(RDF0) << endl;
    // cout <<"asin(RDF0) / arma::norm(drc): " << asin(RDF0) / arma::norm(drc) << endl; 
    // cout << "fabs(asin(RDF0) / arma::norm(drc)): " << fabs(asin(RDF0) / arma::norm(drc)) << endl;
    // cout << "dphi: " << dphi << endl;
    if (dphi<abs(asin((RDF0)/arma::norm(rDFc - rAcm.col(clusterADNum)))))
    {
        flagAttInSight1 = 1;
    } else {
        flagAttInSight1 = 0;
    }
    // cout <<"in defDesiredOpenForm flagAttInSight1: " << flagAttInSight1 << endl;
    double phi_ddot0=-kDFphi*dphi-kDFphid*(phi_dot-phi_des_dot);
    double phi_ddot=arma::sign(phi_ddot0)* std::min(fabs(phi_ddot0),phi_ddot_max);

    // Transiational part
    mat uDFc0;
    if (arma::norm(rDFc - rAcm.col(clusterADNum)) < 3*rho_SN(clusterADNum))
    {
        uDFc0 = -kDFphir*(rDFc-rAcm.col(clusterADNum))-kDFphiv*(vDFc-vAcm.col(clusterADNum));
    } else
    {
        uDFc0=-kDFphir*(rDFc-rAcm.col(clusterADNum));
    }
    
    mat uDFcOv = zeros<mat>(2,1);
    mat uDFcOr = uDFcOv;
    // Avoid other attackers' clusters
    mat nabla_rj_VjAc = zeros<mat>(2,1);
    mat dvAc = zeros<mat>(2,1);

    for (int c = 0; c < NClusterAD; c++)
    {
        if (c != clusterADNum)
        {
            double thetaDFcAcm = atan2(rDFc(1) - rAcm(1,c), rDFc(0)-rAcm(0,c));
            mat tempM(2,1);
            tempM(0,0) = cos(thetaDFcAcm);
            tempM(1,0) = sin(thetaDFcAcm);
            mat rDFcProj = rAcm.col(c) + rhoA_con(c) * tempM;
            mat rTP(2,1);
            rTP(0,0) = cos(thetaDFcAcm + M_PI/2);
            rTP(1,0) = sin(thetaDFcAcm + M_PI/2);
            tempM = rTP.t()*vDFc;
            if (tempM(0,0) < 0)
            {
                rTP = - rTP;
            }
            mat vDFcProj = rTP.t()*vDFc*rTP + vAcm.col(c);
            double Rjk0 = Rjk00(0);
            double Rj_jk = arma::norm(rDFc - rDFcProj);
            double R_underbar = R_m_DD + RDF0;
            if (Rj_jk - R_underbar > 1e-1)
            {
                nabla_rj_VjAc = nabla_rj_VjAc + kDOr*(rDFc-rDFcProj)/Rj_jk/abs(Rj_jk- R_underbar)*(pow((Rj_jk- R_underbar),2)-pow(Rjk0,2))/(pow((Rj_jk- R_underbar),2)+pow(Rjk0,2));
            } else
            {
                nabla_rj_VjAc = nabla_rj_VjAc - kDOr*(rDFc-rDFcProj)*largeP;
            }
            if (arma::norm(vDFc - vDFcProj) != 0)
            {
                dvAc = dvAc - kDOv2*(vDFc-vDFcProj)*pow(arma::norm(vDFc-vDFcProj),(alphaDOv-1));
            } else 
            {
                dvAc = dvAc + 0;
            }
        }
        
    }
    
    mat nabla_rj_VjDc = zeros<mat>(2,1);
    mat dvDc = zeros<mat>(2,1);

    for (int c = 0; c < NClusterD; c++)
    {
        if(c != clusterDNum)
        {
            mat drDrD = XDFc(clusterDNum).submat(0,0,1,0) - XDFc(c).submat(0,0,1,0);
            mat dvDvD = XDFc(clusterDNum).submat(2,0,3,0) - XDFc(c).submat(2,0,3,0);
            double Rjk0 = Rjk00(0);
            double Rj_j = arma::norm(drDrD);
            double R_underbar = R_m_DD + M_PI/2 * (RDF_OPEN(clusterDNum)) + RDF_OPEN(c);
            if(Rj_j - R_underbar > 1e-1)
            {
                nabla_rj_VjDc = nabla_rj_VjDc + kDOr*(drDrD)/Rj_j/abs(Rj_j- R_underbar)*(pow((Rj_j- R_underbar),2)-pow(Rjk0,2))/(pow((Rj_j- R_underbar),2)+pow(Rjk0,2));
            } else
            {
                nabla_rj_VjDc = nabla_rj_VjDc - kDOr*(drDrD)*largeP;
            }

            if(arma::norm(dvDvD) != 0)
            {
                dvDc = dvDc - kDOv2*(dvDvD)*pow(amra::norm(dvDvD),(alphaDOv-1));
            } else
            {
                dvDc = dvDc + 0;
            }
        }
    }

    mat uDFc_trans = uDFc0+uDFcOr+uDFcOv+dvAc-nabla_rj_VjAc+dvDc-nabla_rj_VjDc;
    
    double infNorm_uDFc_trans=max(fabs(uDFc_trans(0)),fabs(uDFc_trans(1)));
    if(infNorm_uDFc_trans > chi*u_maxD(0))
    {
        uDFc_trans = uDFc_trans*chi*u_maxD(0)/infNorm_uDFc_trans;
    }

    double phi1 = phi + M_PI/2;
    double phi2 = phi + 3*M_PI/2;
    double dR = 2*RDF0/(NDinCluster-1);
    // mat uDFcr0 = -kDFphir * drc;
    // mat uDFcOv(2,1,fill::zeros);
    // mat uDFcOr(2,1,fill::zeros);

    // mat uDFc_trans1 = uDFcr0 + uDFcOr + uDFcOv;
    // double norm_uDFc_trans1 = arma::norm(uDFc_trans1);
    // mat uDFc_trans2 = -kDFphiv*(dvc)+C_d*vDFc*arma::norm(vDFc);
    // double norm_uDFc_trans2=arma::norm(uDFc_trans2);

    // mat uDFc_trans0;
    // uDFc_trans0 =std::min(umdf_s1,norm_uDFc_trans1)*uDFc_trans1/norm_uDFc_trans1 + std::min(umdf_s2,norm_uDFc_trans2)*uDFc_trans2/norm_uDFc_trans2;

    // // rDFc.print("rDFc: ");
    // // vDFc.print("vDFc: ");
    // mat uD;
    // mat XD_des1;
    // mat XD_des_dot1;
    // for (int j = 0; j < ND; j++)
    // {
    //     double phi_j = phi + M_PI/2;
    //     mat M1(2,1);
    //     mat M2(2,1);
    //     M1(0,0) = -sin(phi_j);
    //     M1(1,0) = cos(phi_j);
    //     M2(0,0) = cos(phi_j);
    //     M2(1,0) = sin(phi_j);
    //     mat uD_rot = RDF0*phi_ddot1*M1-pow(RDF0*phi_dot,2)*M2;    
    //     uD.insert_cols(j, uDFc_trans0+uD_rot);
    //     vec XD_des_v(4);
    //     // M2.print("M2:");
    //     XD_des_v.subvec(0,1) = rDFc+(0.5*(ND-1)*RDF0-(j)*RDF0)*M2;
    //     // XD_des_v.print("XD_des_v,1,2:");
    //     XD_des_v.subvec(2,3) = vDFc+(0.5*(ND-1)*RDF0-(j)*RDF0)*phi_dot*M1;
    //     // XD_des_v.print("XD_des_v,3,4:");
    //     vec XD_des_dot_v(4);
    //     XD_des_dot_v.subvec(0,1) = XD_des_v.subvec(2,3);
    //     XD_des_dot_v.subvec(2,3) = uD.col(j);
    //     XD_des1.insert_cols(j,XD_des_v);
    //     XD_des_dot1.insert_cols(j, XD_des_dot_v);
    // }
    // // uD.print("uD:");
    // // XD_des1.print("XD_des:");
    // // XD_des_dot1.print("XD_des_dot:");
    // // cout << "phi_dot" << phi_dot << endl;
    // // cout << "phi_ddot" << phi_ddot0 << endl;
    // // uDFc_trans.print("uDFc_trans:");
    // *XD_des = XD_des1;
    // *XD_des_dot = XD_des_dot1;
    // *phi_ddot = phi_ddot1;
    // *uDFc_trans = uDFc_trans0;
    // *flagAttInSight = flagAttInSight1;
    // *XD_des = XD_des1;

}

