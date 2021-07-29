#pragma once
#include "AllParametersExperiment.hpp"
#include <armadillo>
using namespace arma;
// using namespace std;


void defDesiredOpenForm(mat XDFc,double RDF0,mat XA, double phi, double phi_dot,int Na,int ND) {
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
    int flagAttInSight = 0;
    if (dphi < fabs(asin(RDF0) / arma::norm(drc)))
    {
        flagAttInSight = 1;
    } else {
        flagAttInSight = 0;
    }
    
    double phi_ddot0=-kDFphi*dphi-kDFphid*(phi_dot-phi_des_dot);
    double phi_ddot=arma::sign(phi_ddot0)* std::min(fabs(phi_ddot0),phi_ddot_max);

}