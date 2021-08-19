#pragma once
#include "findShortestPath.cpp"
#include "findPathSpeeds.cpp"
#include "AllParameters.hpp"

#include <armadillo>

using namespace std;
using namespace arma;

vec findTimeOnPath(vec S, path_elem Path, pathVel_elem pathVel) {
    vec v = pathVel.v;
    double u_maxD_ = pathVel.u_maxD;
    double v_maxD_ = pathVel.v_maxD;
    double v_maxDC_ = pathVel.v_maxDC;
    double s_bar1 = pathVel.s_bar1;
    double s_bar2 = pathVel.s_bar2;
    double T_bar1 = pathVel.T_bar1;
    double T_bar2 = pathVel.T_bar2;
    double v_bar = pathVel.v_bar;

    double lambda0=sqrt(pow(C_d,2)*pow(rho_safe,2)+4);
    double kappa1=rho_safe*u_maxD_*lambda0;
    double kappa2=C_d*pow(rho_safe,2)*u_maxD_;

    double lambda1=sqrt(rho_safe*(lambda0-rho_safe*C_d));
    double lambda2=sqrt(rho_safe*(lambda0+rho_safe*C_d));
    double lambda3=lambda0/rho_safe*sqrt(u_maxD_/2);
    vec T;
    for (int i = 0; i < S.n_elem; i++)
    {
        arma::uvec ind0 = arma::find(Path.S < S(i));
        int segId;
        if(! ind0.is_empty()) {
            segId = ind0(ind0.n_elem-1);
        } else {
            segId = 1;
        }
        double S0 = Path.S(segId-1);
        double dS = S(i) - S0;
        double dT_bar1=T_bar1-pathVel.T(segId-1);
        double dT_bar2=T_bar2-pathVel.T(segId-1);
        double v2;
        double dT;
        double vtemp;
        if (segId % 2 == 1) {
            if(dS <= s_bar1)
            {
                vtemp=(u_maxD_-exp(-2*C_d*dS)*(u_maxD_-C_d*pow(v(segId-1),2)))/C_d;
                if (vtemp < 0) 
                {
                    v2 = imag(sqrt(vtemp));
                } else
                {
                    v2 = sqrt(vtemp);
                }
                dT = 1/sqrt(u_maxD_*C_d)*(atanh(v2/v_maxD_)-atanh(v(segId-1)/v_maxD_));

            } else if (dS > s_bar1 && dS < s_bar2) 
            {
                dT = dT_bar1+dS/v_bar;
            } else 
            {
                vtemp=(-u_maxD_+exp(-2*C_d*(dS-s_bar2))*(u_maxD_+C_d*pow(v_bar,2)))/C_d;
                if (vtemp < 0) 
                {
                    v2 = imag(sqrt(vtemp));
                } else
                {
                    v2 = sqrt(vtemp);
                }
                dT = dT_bar2+ 1/sqrt(u_maxD_*C_d)*(atan(v_bar/v_maxD_)-atan(v2)/v_maxD_);
            }
        } else {
            if (dS < s_bar1)
            {
                v2=0.5*sqrt(kappa1*tanh(dS*lambda0/rho_safe+atanh((kappa2+2*pow(v(segId-1),2))/kappa1))-kappa2);
                dT=1/lambda3*(atan(sqrt(2/u_maxD_)*v2/lambda1)/lambda1+atanh(sqrt(2/u_maxD_)*v2/lambda2)/lambda2)
                -1/lambda3*(atan(sqrt(2/u_maxD_)*v(segId-1)/lambda1)/lambda1+atanh(sqrt(2/u_maxD_)*v(segId-1)/lambda2)/lambda2);
            } else if (dS > s_bar1 && dS < s_bar2)
            {
                dT = dS/v_maxDC_;
            } else
            {
                v2=0.5*sqrt(kappa2-kappa1*tanh(-dS*lambda0/rho_safe+atanh((kappa2-2*pow(v(segId-1),2))/kappa1)));
                dT=1/lambda3*(atan(sqrt(2/u_maxD_)*v2/lambda2)/lambda2+atanh(sqrt(2/u_maxD_)*v2/lambda1)/lambda1)
                -1/lambda3*(atan(sqrt(2/u_maxD_)*v(segId-1)/lambda2)/lambda2+atanh(sqrt(2/u_maxD_)*v(segId-1)/lambda1)/lambda1) ;
            }
        }
        T.resize(i+1);
        T(i) = pathVel.T(segId-1) + dT;
    }

    return T;

    
}