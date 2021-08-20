#pragma once
#include "findShortestPath.cpp"
#include "findPathSpeeds.cpp"

mat findPosOnPath(vec T, path_elem path, pathVel_elem pathVel){
    mat S = zeros<mat>(1, T.n_elem);
    vec v = pathVel.v;
    double v_bar = pathVel.v_bar;
    double u_maxD_=pathVel.u_maxD;
    double v_maxD_=pathVel.v_maxD;
    double v_maxDC_=pathVel.v_maxDC;
    double s_bar1=pathVel.s_bar1;
    double s_bar2=pathVel.s_bar2;
    double T_bar1=pathVel.T_bar1;
    double T_bar2=pathVel.T_bar2;

    double lambda0=sqrt(pow(C_d,2)*pow(rho_safe,2)+4);
    double kappa1=rho_safe*u_maxD_*lambda0;
    double kappa2=C_d*pow(rho_safe,2)*u_maxD_;

    double lambda1=sqrt(rho_safe*(lambda0-rho_safe*C_d));
    double lambda10=sqrt(2/u_maxD_)/lambda1;
    double lambda2=sqrt(rho_safe*(lambda0+rho_safe*C_d));
    double lambda20=sqrt(2/u_maxD_)/lambda2;
    double lambda3=lambda0/rho_safe*sqrt(u_maxD_/2);
    double vT, dS;
    int segId;
    for (int i = 0; i < T.n_elem; i++)
    {
        uvec ind0 = arma::find(pathVel.T < T(i));
        
        if (!ind0.is_empty()) 
        {
            segId = ind0(ind0.n_elem-1);
        } else 
        {
            segId = 1;
        }
        double T0 = pathVel.T(segId-1);
        double dT = T(i) - T0;
        if (segId%2 == 1)
        {
            if(T(i) <= T_bar1){
                vT = sqrt(u_maxD_/C_d)*tanh(sqrt(u_maxD_*C_d)*dT+atanh(sqrt(C_d/u_maxD_)*v(segId-1)));
                dS=1/(2*C_d)*log((u_maxD_-C_d*pow(v(segId-1),2))/(u_maxD_-C_d*pow(vT,2)));
            } else if (T(i) > T_bar1 && T(i) < T_bar2)
            {
                dS=s_bar1+v_maxD_*(T(i)-T_bar1);
            } else
            {
                vT=sqrt(u_maxD_/C_d)*tan(sqrt(u_maxD_*C_d)*dT+atan(sqrt(C_d/u_maxD_)*v_bar));
                dS=s_bar2+1/(2*C_d)*log((u_maxD_+C_d*pow(v_bar,2))/(u_maxD_+pow(C_d*vT,2)));
            }
        } // skip the part if the segment is circular arc
        S(i) = path.S(segId-1) + dS;
    }
    return S;
}