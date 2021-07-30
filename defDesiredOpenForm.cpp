#pragma once
#include "AllParametersExperiment.hpp"
#include <armadillo>
using namespace arma;
// using namespace std;


void defDesiredOpenForm(mat XDFc,double RDF0,mat XA, double phi, double phi_dot,int Na,int ND, mat *XD_des, mat *XD_des_dot, double *phi_ddot, mat *uDFc_trans, int *flagAttInSight){
                        // mat *XD_des, mat *XD_des_dot, double *phi_ddot,
                        // mat *uDFc_trans, int *flagAttInSight) {
    double chi = 0.55;
    double phi_ddot_max=(0.9-chi)*u_maxD(1)/RDF0;

    mat rAcm,vAcm;
    rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1)/NA;
    vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1)/NA;

    mat rDFc,vDFc;
    rDFc = XDFc.submat(0,0,1,0);
    vDFc = XDFc.submat(2,0,3,0);

    mat drc,dvc;
    drc = rDFc - rAcm;
    dvc = vDFc - vAcm;

    double phi_des = atan2(-drc(1), -drc(0));
    if(phi_des < 0) {
        phi_des += 2*M_PI;
    }
    double phi_des_dot=(-drc(0)*dvc(1)+drc(1)*dvc(0))/pow(arma::norm(drc),2);
    double dphi = (phi - phi_des);

    if (dphi < -M_PI) {
        dphi += 2*M_PI;
    }
    // int flagAttInSight = 0;
    if (dphi < fabs(asin(RDF0) / arma::norm(drc)))
    {
        *flagAttInSight = 1;
    } else {
        *flagAttInSight = 0;
    }
    
    double phi_ddot0=-kDFphi*dphi-kDFphid*(phi_dot-phi_des_dot);
    double phi_ddot1=arma::sign(phi_ddot0)* std::min(fabs(phi_ddot0),phi_ddot_max);

    mat uDFcr0 = -kDFphir * drc;
    mat uDFcOv(2,1,fill::zeros);
    mat uDFcOr(2,1,fill::zeros);

    mat uDFc_trans1 = uDFcr0 + uDFcOr + uDFcOv;
    double norm_uDFc_trans1 = arma::norm(uDFc_trans1);
    mat uDFc_trans2 = -kDFphiv*(dvc)+C_d*vDFc*arma::norm(vDFc);
    double norm_uDFc_trans2=arma::norm(uDFc_trans2);

    mat uDFc_trans0;
    uDFc_trans0 =std::min(umdf_s1,norm_uDFc_trans1)*uDFc_trans1/norm_uDFc_trans1 + std::min(umdf_s2,norm_uDFc_trans2)*uDFc_trans2/norm_uDFc_trans2;

    // rDFc.print("rDFc: ");
    // vDFc.print("vDFc: ");
    mat uD;
    mat XD_des1;
    mat XD_des_dot1;
    for (int j = 0; j < ND; j++)
    {
        double phi_j = phi + M_PI/2;
        mat M1(2,1);
        mat M2(2,1);
        M1(0,0) = -sin(phi_j);
        M1(1,0) = cos(phi_j);
        M2(0,0) = cos(phi_j);
        M2(1,0) = sin(phi_j);
        mat uD_rot = RDF0*phi_ddot1*M1-pow(RDF0*phi_dot,2)*M2;    
        uD.insert_cols(j, uDFc_trans0+uD_rot);
        vec XD_des_v(4);
        // M2.print("M2:");
        XD_des_v.subvec(0,1) = rDFc+(0.5*(ND-1)*RDF0-(j)*RDF0)*M2;
        // XD_des_v.print("XD_des_v,1,2:");
        XD_des_v.subvec(2,3) = vDFc+(0.5*(ND-1)*RDF0-(j)*RDF0)*phi_dot*M1;
        // XD_des_v.print("XD_des_v,3,4:");
        vec XD_des_dot_v(4);
        XD_des_dot_v.subvec(0,1) = XD_des_v.subvec(2,3);
        XD_des_dot_v.subvec(2,3) = uD.col(j);
        XD_des1.insert_cols(j,XD_des_v);
        XD_des_dot1.insert_cols(j, XD_des_dot_v);
    }
    // uD.print("uD:");
    // XD_des1.print("XD_des:");
    // XD_des_dot1.print("XD_des_dot:");
    // cout << "phi_dot" << phi_dot << endl;
    // cout << "phi_ddot" << phi_ddot0 << endl;
    // uDFc_trans.print("uDFc_trans:");
    *XD_des = XD_des1;
    *XD_des_dot = XD_des_dot1;
    *phi_ddot = phi_ddot1;
    *uDFc_trans = uDFc_trans0;
    // *XD_des = XD_des1;

}

