#include "AllParameters.hpp"
#include <armadillo>

using namespace arma;

// void defDesiredClosedForm(mat XDFc, double RDF0,double phi0, mat XA, int NA, int ND, int flagNDeven, double delta_t, mat *XD_des, mat *XD_des_dot, mat *uDFc_trans){
//     mat rDFc, vDFc;
//     rDFc = XDFc.submat(0,0,1,0);
//     vDFc = XDFc.submat(2,0,3,0);

//     mat uDFc0;
//     if (arma::norm(rDFc - rS) < 3*rho_sn)
//     {
//         uDFc0 = -kDFphir*(rDFc-rS)-kDFphiv*(vDFc);
//     } else {
//         uDFc0 = -kDFphir*(rDFc-rS);
//     }

//     mat uDFcOv(2,1,fill::zeros);
//     mat uDFcOr(2,1,fill::zeros);

//     double R_DFcP = arma::norm(rDFc - rP) - rho_P;
//     double theta = atan2(rDFc(1)-rP(1),rDFc(0)-rP(0));
//     mat rDFcProjP(2,1);
//     rDFcProjP(0,0) = rho_P*cos(theta);
//     rDFcProjP(1,0) = rho_P*sin(theta);
//     mat rTP(2,1);
//     rTP(0,0) = cos(theta+M_PI/2);
//     rTP(1,0) = sin(theta+M_PI/2);
//     mat tempM;
//     tempM = rTP.t()*vDFc;
//     if (tempM(0,0) < 0)
//     {
//         rTP = - rTP;
//     }
//     tempM = rTP.t()*vDFc;
//     mat vDFcProjP = tempM(0,0)*rTP;
//     mat uDFc_trans1=uDFc0+uDFcOr+uDFcOv;
//     double chi2=0.7*u_maxA(0)/u_maxD(0);
//     double chi1=0.0;
//     double delta_t_max1=0;
//     double delta_t_max2=20;
//     double delta_t_max=delta_t_max2-delta_t_max1;
//     double chi;
//     if (delta_t < delta_t_max1)
//     {
//         chi = chi1;
//     } else if (delta_t_max1<delta_t && delta_t<delta_t_max2)
//     {
//         chi=chi1+(chi2-chi1)*(delta_t-delta_t_max1)/delta_t_max;
//     } else 
//     {
//         chi=chi2;
//     }
//     double infNorm_uDFc_trans=std::max(abs(uDFc_trans1(0)),abs(uDFc_trans1(1)));
//     double norm_uDFc_trans=arma::norm(uDFc_trans1);
//     if (norm_uDFc_trans>chi*u_maxD(0))
//     {
//         uDFc_trans1=uDFc_trans1*chi*u_maxD(0)/norm_uDFc_trans;
//     }
//     mat XD_des_1(4,ND);
//     mat XD_des_dot_1(4,ND);

//     for (int j = 0; j < ND; j++)
//     {
//         double phi_j;
//         phi_j = phi0 + 2*M_PI*(j+1)/ND - M_PI/ND;
        
//         tempM.resize(2,1);
//         tempM(0,0) = cos(phi_j);
//         tempM(1,0) = sin(phi_j);
//         XD_des_1.submat(0,j,1,j) = rDFc + RDF0*tempM;
//         XD_des_1.submat(2,j,3,j) = vDFc;
//         XD_des_dot_1.submat(2,j,3,j) = uDFc_trans1;
//         XD_des_dot_1.submat(0,j,1,j) = XD_des_1.submat(2,j,3,j);
//     }
    

//     *XD_des = XD_des_1;
//     *XD_des_dot = XD_des_dot_1;
//     *uDFc_trans = uDFc_trans1;

// }

void defDesiredClosedForm(std::vector<mat> XDFc, int clusterDNum, mat rS, mat RDF_CLOSED, std::vector<arma::vec> indDes, int NDinCluster, int NClusterD,double phi0, mat XA, mat rAcm, mat vAcm, mat rhoA_con, int NClusterAD,
                          int clusterADNum,  int NA, int ND, int flagNDeven, int flagAvoidObs,  double delta_t, mat &XD_des, mat &XD_des_dot, mat &uDFc_trans){
    double RDF0 = RDF_CLOSED(clusterDNum-1);
    
    mat rDFc, vDFc;
    rDFc = XDFc[clusterDNum-1].submat(0,0,1,0);
    vDFc = XDFc[clusterDNum-1].submat(2,0,3,0);

    mat uDFc0;
    if (arma::norm(rDFc - rS) < 3*rho_sn)
    {
        uDFc0 = -kDFphir*(rDFc-rS)-kDFphiv*(vDFc);
    } else {
        uDFc0 = -kDFphir*(rDFc-rS);
    }
    
    mat uDFcOv(2,1,fill::zeros);
    mat uDFcOr(2,1,fill::zeros);

    mat nabla_rj_VjAc(2,1,fill::zeros);
    mat dvAc(2,1,fill::zeros);

    for(int c = 0; c < NClusterAD; c++)
    {
        if(c != clusterADNum - 1)
        {
            double thetaDFcAcm = atan2(rDFc(1) - rAcm(1,c), rDFc(0) - rAcm(0,c));
            mat tempM(2,1);
            tempM(0) = cos(thetaDFcAcm);
            tempM(1) = sin(thetaDFcAcm);
            mat rDFcProj = rAcm.col(c) + rhoA_con(c) * tempM;
            mat rTP(2,1);
            rTP(0) = cos(thetaDFcAcm + M_PI / 2);
            rTP(1) = sin(thetaDFcAcm + M_PI / 2);
            mat result = rTP.t() * vDFc;
            if(result(0,0) < 0)
            {
                rTP = -rTP;
            }
            mat vDFcProj = result(0,0) * rTP + vAcm.col(c);
            double Rjk0 = Rjk00(0);
            double Rj_jk = arma::norm(rDFc - rDFcProj);
            double R_underbar = R_m_DD + RDF0;
            if (Rj_jk - R_underbar > 1e-1)
            {
                nabla_rj_VjAc = nabla_rj_VjAc + kDOr*(rDFc-rDFcProj)/Rj_jk/fabs(Rj_jk- R_underbar)*(pow((Rj_jk- R_underbar),2)-pow(Rjk0,2))/(pow((Rj_jk- R_underbar),2)+pow(Rjk0,2));
            } else
            {
                nabla_rj_VjAc = nabla_rj_VjAc - kDOr*(rDFc - rDFcProj) * largeP;
            }

            if(arma::norm(vDFc - vDFcProj) != 0)
            {
                dvAc =  dvAc - kDOv2*(vDFc-vDFcProj)*pow(arma::norm(vDFc-vDFcProj),(alphaDOv-1));
            }
            


        }
    }

    mat nabla_rj_VjDc(2,1,fill::zeros);
    mat dvDc(2,1,fill::zeros);

    for(int c = 0; c < NClusterD; c++)
    {
        if(c != clusterDNum)
        {
            mat drDrD = XDFc[clusterDNum].submat(0,0,1,0) - XDFc[c].submat(0,0,1,0);
            mat dvDvD = XDFc[clusterDNum].submat(2,0,3,0) - XDFc[c].submat(2,0,3,0);
            double Rjk0 = Rjk00(0);
            double Rj_j = arma::norm(drDrD);
            double R_underbar = R_m_DD + M_PI/2 * (RDF_CLOSED(clusterDNum)) + RDF_CLOSED(c);
            if(Rj_j - R_underbar > 1e-1)
            {
                nabla_rj_VjDc = nabla_rj_VjDc + kDOr * (drDrD) / Rj_j/ fabs(Rj_j - R_underbar) * (pow((Rj_j- R_underbar),2)-pow(Rjk0,2))/(pow((Rj_j- R_underbar),2)+pow(Rjk0,2));
            } else
            {
                nabla_rj_VjDc = nabla_rj_VjDc - kDOr*(drDrD)*largeP;
            }

            if(arma::norm(dvDvD) != 0)
            {
                dvDc = dvDc - kDOv2*(dvDvD)*pow(arma::norm(dvDvD),(alphaDOv-1));
            }
        }
    }

    uDFc_trans=uDFc0+uDFcOr+uDFcOv+dvAc-nabla_rj_VjAc+dvDc-nabla_rj_VjDc ;

    double chi2 = 0.5 * u_maxA(0) / u_maxD(0);
    double chi1 = 0.0;
    double delta_t_max1 = 10;
    double delta_t_max2 = 30;
    double delta_t_max = delta_t_max2 - delta_t_max1;
    double chi;
    if(delta_t < delta_t_max1)
    {
        chi = chi1;
    } else if (delta_t_max1 < delta_t && delta_t < delta_t_max2)
    {
        chi = chi1+(chi2-chi1)*(delta_t-delta_t_max1)/delta_t_max;
    } else
    {
        chi = chi2;
    }
    

    double infNorm_uDFc_trans = max(fabs(uDFc_trans(0)), fabs(uDFc_trans(1)));
    double norm_uDFc_trans = arma::norm(uDFc_trans);
    if(norm_uDFc_trans > chi * u_maxD(0))
    {
        uDFc_trans = uDFc_trans*chi*u_maxD(0)/norm_uDFc_trans;
    }

    for(int jj = 0; jj < NDinCluster; jj++)
    {
        int j = indDes[clusterDNum](jj);
        double phi_j;
        if(flagNDeven != 1)
        {
            phi_j = phi0 + 2 * M_PI * (jj+1) / NDinCluster - M_PI / NDinCluster;
        } else
        {
            phi_j = phi0 + 2 * M_PI * (jj+1) / NDinCluster - M_PI / NDinCluster;
        }
        mat tempM(2,1);
        tempM(0) = cos(phi_j);
        tempM(1) = sin(phi_j);
        XD_des.submat(0,j,1,j) = rDFc + RDF0 * tempM;
        XD_des.submat(2,j,3,j) = vDFc;
        XD_des_dot.submat(2,j,3,j) = uDFc_trans;
        XD_des_dot.submat(0,j,1,j) = XD_des_dot.submat(2,j,3,j);
    }
    // double R_DFcP = arma::norm(rDFc - rP) - rho_P;
    // double theta = atan2(rDFc(1)-rP(1),rDFc(0)-rP(0));
    // mat rDFcProjP(2,1);
    // rDFcProjP(0,0) = rho_P*cos(theta);
    // rDFcProjP(1,0) = rho_P*sin(theta);
    // mat rTP(2,1);
    // rTP(0,0) = cos(theta+M_PI/2);
    // rTP(1,0) = sin(theta+M_PI/2);
    // mat tempM;
    // tempM = rTP.t()*vDFc;
    // if (tempM(0,0) < 0)
    // {
    //     rTP = - rTP;
    // }
    // tempM = rTP.t()*vDFc;
    // mat vDFcProjP = tempM(0,0)*rTP;
    // mat uDFc_trans1=uDFc0+uDFcOr+uDFcOv;
    // double chi2=0.7*u_maxA(0)/u_maxD(0);
    // double chi1=0.0;
    // double delta_t_max1=0;
    // double delta_t_max2=20;
    // double delta_t_max=delta_t_max2-delta_t_max1;
    // double chi;
    // if (delta_t < delta_t_max1)
    // {
    //     chi = chi1;
    // } else if (delta_t_max1<delta_t && delta_t<delta_t_max2)
    // {
    //     chi=chi1+(chi2-chi1)*(delta_t-delta_t_max1)/delta_t_max;
    // } else 
    // {
    //     chi=chi2;
    // }
    // double infNorm_uDFc_trans=std::max(abs(uDFc_trans1(0)),abs(uDFc_trans1(1)));
    // double norm_uDFc_trans=arma::norm(uDFc_trans1);
    // if (norm_uDFc_trans>chi*u_maxD(0))
    // {
    //     uDFc_trans1=uDFc_trans1*chi*u_maxD(0)/norm_uDFc_trans;
    // }
    // mat XD_des_1(4,ND);
    // mat XD_des_dot_1(4,ND);

    // for (int j = 0; j < ND; j++)
    // {
    //     double phi_j;
    //     phi_j = phi0 + 2*M_PI*(j+1)/ND - M_PI/ND;
        
    //     tempM.resize(2,1);
    //     tempM(0,0) = cos(phi_j);
    //     tempM(1,0) = sin(phi_j);
    //     XD_des_1.submat(0,j,1,j) = rDFc + RDF0*tempM;
    //     XD_des_1.submat(2,j,3,j) = vDFc;
    //     XD_des_dot_1.submat(2,j,3,j) = uDFc_trans1;
    //     XD_des_dot_1.submat(0,j,1,j) = XD_des_1.submat(2,j,3,j);
    // }
    

    // *XD_des = XD_des_1;
    // *XD_des_dot = XD_des_dot_1;
    // *uDFc_trans = uDFc_trans1;

}