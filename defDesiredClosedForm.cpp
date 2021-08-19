#include "AllParameters.hpp"
#include <armadillo>

using namespace arma;

void defDesiredClosedForm(mat XDFc, double RDF0,double phi0, mat XA, int NA, int ND, int flagNDeven, double delta_t, mat *XD_des, mat *XD_des_dot, mat *uDFc_trans){
    mat rDFc, vDFc;
    rDFc = XDFc.submat(0,0,1,0);
    vDFc = XDFc.submat(2,0,3,0);

    mat uDFc0;
    if (arma::norm(rDFc - rS) < 3*rho_sn)
    {
        uDFc0 = -kDFphir*(rDFc-rS)-kDFphiv*(vDFc);
    } else {
        uDFc0 = -kDFphir*(rDFc-rS);
    }

    mat uDFcOv(2,1,fill::zeros);
    mat uDFcOr(2,1,fill::zeros);

    double R_DFcP = arma::norm(rDFc - rP) - rho_P;
    double theta = atan2(rDFc(1)-rP(1),rDFc(0)-rP(0));
    mat rDFcProjP(2,1);
    rDFcProjP(0,0) = rho_P*cos(theta);
    rDFcProjP(1,0) = rho_P*sin(theta);
    mat rTP(2,1);
    rTP(0,0) = cos(theta+M_PI/2);
    rTP(1,0) = sin(theta+M_PI/2);
    mat tempM;
    tempM = rTP.t()*vDFc;
    if (tempM(0,0) < 0)
    {
        rTP = - rTP;
    }
    tempM = rTP.t()*vDFc;
    mat vDFcProjP = tempM(0,0)*rTP;
    mat uDFc_trans1=uDFc0+uDFcOr+uDFcOv;
    double chi2=0.7*u_maxA(0)/u_maxD(0);
    double chi1=0.0;
    double delta_t_max1=0;
    double delta_t_max2=20;
    double delta_t_max=delta_t_max2-delta_t_max1;
    double chi;
    if (delta_t < delta_t_max1)
    {
        chi = chi1;
    } else if (delta_t_max1<delta_t && delta_t<delta_t_max2)
    {
        chi=chi1+(chi2-chi1)*(delta_t-delta_t_max1)/delta_t_max;
    } else 
    {
        chi=chi2;
    }
    double infNorm_uDFc_trans=std::max(abs(uDFc_trans1(0)),abs(uDFc_trans1(1)));
    double norm_uDFc_trans=arma::norm(uDFc_trans1);
    if (norm_uDFc_trans>chi*u_maxD(0))
    {
        uDFc_trans1=uDFc_trans1*chi*u_maxD(0)/norm_uDFc_trans;
    }
    mat XD_des_1(4,ND);
    mat XD_des_dot_1(4,ND);

    for (int j = 0; j < ND; j++)
    {
        double phi_j;
        phi_j = phi0 + 2*M_PI*(j+1)/ND - M_PI/ND;
        
        tempM.resize(2,1);
        tempM(0,0) = cos(phi_j);
        tempM(1,0) = sin(phi_j);
        XD_des_1.submat(0,j,1,j) = rDFc + RDF0*tempM;
        XD_des_1.submat(2,j,3,j) = vDFc;
        XD_des_dot_1.submat(2,j,3,j) = uDFc_trans1;
        XD_des_dot_1.submat(0,j,1,j) = XD_des_1.submat(2,j,3,j);
    }
    

    *XD_des = XD_des_1;
    *XD_des_dot = XD_des_dot_1;
    *uDFc_trans = uDFc_trans1;

}