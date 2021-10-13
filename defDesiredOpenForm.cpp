#pragma once
#include "AllParameters.cpp"
#include <armadillo>
using namespace arma;
using namespace std;

struct OpenForm
{
    double phi_dot;
    double phi_ddot;
    mat uDFc_trans;
    int flagAttInSight;
};


OpenForm defDesiredOpenForm(std::vector<mat> XDFc, int clusterDNum, mat *XD_des, mat *XD_des_dot, vec RDF_OPEN, mat rho_SN, std::vector<vec> indDes, int NDinCluster, int NClusterD, mat XA,
                        mat rAcm, mat vAcm, mat rhoA_con, int NClusterAD, int clusterADNum, double phi, double phi_dot, int NA, int ND, int flagAvoidObs){
                        // mat *XD_des, mat *XD_des_dot, double *phi_ddot,
                        // mat *uDFc_trans, int *flagAttInSight) {
                        // ,mat XA, double phi, double phi_dot,int NA,int ND, mat *XD_des, mat *XD_des_dot, double *phi_ddot, mat *uDFc_trans, int *flagAttInSight
    double chi = 0.55;
    double RDF0 = RDF_OPEN(clusterDNum-1);
    double phi_ddot_max= (0.9-chi)*u_maxD(0)/RDF0;
    // cout << "ND: " << ND << endl;
    mat uD(2, ND, fill::zeros);
    // uD.print("uD: ");
    // mat rAcm,vAcm;
    // rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1)/NA;
    // vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1)/NA;
    // XDFc.print("XDFc: ");
    mat rDFc,vDFc;
    rDFc = XDFc[clusterDNum-1].submat(0,0,1,0);
    vDFc = XDFc[clusterDNum-1].submat(2,0,3,0);
    // rDFc.print("rDFc: ");
    // vDFc.print("vDFc: ");
    mat drc,dvc;
    drc = rAcm.col(clusterADNum-1) - rDFc;
    dvc = vAcm.col(clusterADNum-1) - vDFc;
    // drc.print("drc: ");
    // dvc.print("dvc: ");

    double phi_des = atan2(drc(1), drc(0));
    if(phi_des < 0) {
        phi_des += 2*M_PI;
    }
    // cout << "phi_des: " << phi_des << endl;
    double phi_des_dot=(drc(0)*dvc(1)-drc(1)*dvc(0))/pow(arma::norm(drc),2);
    // cout << "phi_des_dot: " << phi_des_dot << endl;
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
    // cout << "phi_ddot0: " << phi_ddot0 << endl;
    // cout << "phi_ddot: " << phi_ddot << endl;
    // Transiational part
    mat uDFc0;
    if (arma::norm(rDFc - rAcm.col(clusterADNum-1)) < 3*rho_SN(clusterADNum-1))
    {
        uDFc0 = -kDFphir*(rDFc-rAcm.col(clusterADNum-1))-kDFphiv*(vDFc-vAcm.col(clusterADNum-1));
    } else
    {
        uDFc0=-kDFphir*(rDFc-rAcm.col(clusterADNum-1));
    }
    // uDFc0.print("uDFc0: ");
    mat uDFcOv = zeros<mat>(2,1);
    mat uDFcOr = uDFcOv;
    // Avoid other attackers' clusters
    mat nabla_rj_VjAc = zeros<mat>(2,1);
    mat dvAc = zeros<mat>(2,1);

    for (int c = 0; c < NClusterAD; c++)
    {
        if (c != clusterADNum-1)
        {
            double thetaDFcAcm = atan2(rDFc(1) - rAcm(1,c), rDFc(0)-rAcm(0,c));
            // cout << "thetaDFcAcm" << thetaDFcAcm << endl;
            mat tempM(2,1);
            tempM(0,0) = cos(thetaDFcAcm);
            tempM(1,0) = sin(thetaDFcAcm);
            mat rDFcProj = rAcm.col(c) + rhoA_con(c) * tempM;
            // rDFcProj.print("rDFcProj: ");
            mat rTP(2,1);
            rTP(0,0) = cos(thetaDFcAcm + M_PI/2);
            rTP(1,0) = sin(thetaDFcAcm + M_PI/2);
            tempM = rTP.t()*vDFc;
            if (tempM(0,0) < 0)
            {
                rTP = - rTP;
            }
            // rTP.print("rTP: ");
            tempM = rTP.t()*vDFc;
            mat vDFcProj = tempM(0,0) * rTP + vAcm.col(c);
            // vDFcProj.print("vDFcProj: ");
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
    // dvAc.print("dvAc: ");
    mat nabla_rj_VjDc = zeros<mat>(2,1);
    mat dvDc = zeros<mat>(2,1);

    for (int c = 0; c < NClusterD; c++)
    {
        if(c != clusterDNum-1)
        {
            mat drDrD = XDFc[clusterDNum-1].submat(0,0,1,0) - XDFc[c].submat(0,0,1,0);
            mat dvDvD = XDFc[clusterDNum-1].submat(2,0,3,0) - XDFc[c].submat(2,0,3,0);
            double Rjk0 = Rjk00(0);
            double Rj_j = arma::norm(drDrD);
            double R_underbar = R_m_DD + M_PI/2 * (RDF_OPEN(clusterDNum-1)) + RDF_OPEN(c);
            if(Rj_j - R_underbar > 1e-1)
            {
                nabla_rj_VjDc = nabla_rj_VjDc + kDOr*(drDrD)/Rj_j/abs(Rj_j- R_underbar)*(pow((Rj_j- R_underbar),2)-pow(Rjk0,2))/(pow((Rj_j- R_underbar),2)+pow(Rjk0,2));
            } else
            {
                nabla_rj_VjDc = nabla_rj_VjDc - kDOr*(drDrD)*largeP;
            }

            if(arma::norm(dvDvD) != 0)
            {
                dvDc = dvDc - kDOv2*(dvDvD)*pow(arma::norm(dvDvD),(alphaDOv-1));
            } else
            {
                dvDc = dvDc + 0;
            }
        }
    }
    // dvDc.print("dvDc: ");
    mat uDFc_trans = uDFc0+uDFcOr+uDFcOv+dvAc-nabla_rj_VjAc+dvDc-nabla_rj_VjDc;
    // uDFc_trans.print("uDFc_trans: ");
    double infNorm_uDFc_trans=max(fabs(uDFc_trans(0)),fabs(uDFc_trans(1)));
    if(infNorm_uDFc_trans > chi*u_maxD(0))
    {
        uDFc_trans = uDFc_trans*chi*u_maxD(0)/infNorm_uDFc_trans;
    }
    // uDFc_trans.print("uDFc_trans: ");
    double phi1 = phi + M_PI/2;
    double phi2 = phi + 3*M_PI/2;
    double dR = 2*RDF0/(NDinCluster-1);

    for(int jj = 0; jj < NDinCluster; jj++)
    {
        int j = indDes[clusterDNum-1](jj);
        mat tempM1(2, 1); //[cos(phi1), sin(phi1)]'
        mat tempM2(2, 1); //[cos(phi2), sin(phi2)]'
        mat tempM3(2, 1); //[-sin(phi1), cos(phi1)]'
        mat tempM4(2, 1); //[-sin(phi2), cos(phi2)]
        tempM1(0, 0) = cos(phi1);
        tempM1(1, 0) = sin(phi1);
        tempM2(0, 0) = cos(phi2);
        tempM2(1, 0) = sin(phi2);
        tempM3(0, 0) = -sin(phi1);
        tempM3(1, 0) = cos(phi1);
        tempM4(0, 0) = -sin(phi2);
        tempM4(1, 0) = cos(phi2);
        XD_des->submat(0,j,1,j) = rDFc + RDF0 * tempM1 + jj * dR * tempM2;
        XD_des->submat(2,j,3,j) = vDFc + RDF0 * phi_dot * tempM3 + jj * dR * phi_dot * tempM4;
        mat uD_rot = RDF0 * phi_ddot * tempM3 + jj * dR * phi_ddot * tempM4 - RDF0 * pow(phi_dot,2) * tempM1 - jj * dR * pow(phi_dot,2) * tempM2;
        // uD_rot.print("uD_rot");
        // uD.print("uD: ");
        // cout << "j: " << j << endl;
        uD.col(j) = uDFc_trans + uD_rot;
        // cout << "1111111111111" << endl;
        XD_des_dot->submat(2,j,3,j) = uD.col(j);
        // cout << "2222222222222" << endl;
        XD_des_dot->submat(0,j,1,j) = XD_des->submat(2,j,3,j);
        // cout << "3333333333333" << endl;
    }

    OpenForm result;
    result.flagAttInSight = flagAttInSight1;
    result.phi_dot = phi_dot;
    result.phi_ddot = phi_ddot;
    result.uDFc_trans = uDFc_trans;
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
    return result;
}

// int main(){
//     int NA = 18;
//     int ND = 18;

//     calAllparametersExperiment(NA,ND);
//     std::vector<mat> XDFc(3); 
//     int clusterDNum; 
//     mat XD_des;
//     mat XD_des_dot;
//     vec RDF_OPEN; 
//     mat RDF_OPEN_mat;
//     mat rho_SN;
//     std::vector<vec> indDes;
//     int NDinCluster;
//     int NClusterD;
//     mat XA;

//     mat rAcm;
//     mat vAcm;
//     mat rhoA_con;
//     int NClusterAD;
//     int clusterADNum;
//     double phi;
//     double phi_dot;
//     // int NA;
//     // int ND;
//     int flagAvoidObs;

//     mat tempM(4,1, fill::zeros);
//     tempM(0,0) = 385.153071102075;
//     tempM(1,0) =-176.086389266759;
//     XDFc[0] = tempM;
//     tempM(0,0) = 443.375146160435;
//     tempM(1,0) = 37.3666215895702;
//     XDFc[1] = tempM;
//     tempM(0,0) = 512.183053047588;
//     tempM(1,0) = 289.629270783414;
//     XDFc[2] = tempM;
//     // XDFc.print("XDFc: ");

//     clusterDNum = 1;
//     XD_des.load("/home/damon/Downloads/multi_swarm/defDopen/XD_des.txt");
//     // XD_des.print("XD_des: ");
//     XD_des.load("/home/damon/Downloads/multi_swarm/defDopen/XD_des.txt");
//     XD_des_dot.load("/home/damon/Downloads/multi_swarm/defDopen/XD_des_dot.txt");
//     RDF_OPEN_mat.load("/home/damon/Downloads/multi_swarm/defDopen/RDF_open.txt");
//     RDF_OPEN = RDF_OPEN_mat.as_col();
//     rho_SN.load("/home/damon/Downloads/multi_swarm/defDopen/rho_sn.txt");
//     RDF_OPEN.print("RDF_OPEN: ");
//     rho_SN.print("rho_SN");
//     vec tempV;
//     tempV = {14,15,16,17,18};
//     tempV = tempV - 1;
//     indDes.push_back(tempV);
//     tempV = {8,9,10,11,12,13};
//     tempV = tempV - 1;
//     indDes.push_back(tempV);
//     tempV = {1,2,3,4,5,6,7};
//     tempV = tempV - 1;
//     indDes.push_back(tempV);

//     NDinCluster = 5;
//     NClusterD = 3;
//     XA.load("/home/damon/Downloads/multi_swarm/defDopen/XA.txt");
//     rAcm.load("/home/damon/Downloads/multi_swarm/defDopen/rAcm.txt");
//     vAcm.load("/home/damon/Downloads/multi_swarm/defDopen/vAcm.txt");
//     rhoA_con.load("/home/damon/Downloads/multi_swarm/defDopen/rhoA_con.txt");
//     NClusterAD = 3;
//     clusterADNum = 1;
//     phi = -0.266285265448926;
//     phi_dot = 0;
//     NA = 18;
//     ND = 18;
//     flagAvoidObs = 0;
//     OpenForm result = defDesiredOpenForm(XDFc, clusterDNum, &XD_des, &XD_des_dot, RDF_OPEN, rho_SN, indDes, NDinCluster, NClusterD, XA, rAcm, vAcm, rhoA_con, NClusterAD, clusterADNum, phi, phi_dot, NA, ND, flagAvoidObs);
//     XD_des.print("XD_des");
//     XD_des_dot.print("XD_des_dot: ");
//     cout << "phi_ddot: " << result.phi_ddot << endl;
//     return 1;
// }