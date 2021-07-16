#pragma once
#include "AllParametersExperiment.cpp"


using namespace std;
using namespace arma;

struct pathVel_elem
{
    vec v;
    double v_bar;
    double s_bar1;
    double s_bar2;
    mat T;
    double T_bar1;
    double T_bar2;
    double u_maxD;
    double v_maxD;
    double v_maxDC;
};

pathVel_elem findPathSpeeds(path_elem path, double v_maxDC, double u_maxD, double v_maxD) {
    pathVel_elem result;
    double lambda0 = sqrt(pow(C_d,2) * pow(rho_safe,2) + 4);
    double kappa1 = rho_safe * u_maxD * lambda0;
    double kappa2 = C_d * pow(rho_safe,2) * u_maxD;

    double lambda1 = sqrt(rho_safe * (lambda0 - rho_safe * C_d));
    double lambda2 = sqrt(rho_safe * (lambda0 + rho_safe * C_d));
    double lambda3 = lambda0 / rho_safe * sqrt(u_maxD/2);
    vec S = path.S;
    double P = path.P;
    int NS = S.n_elem - 1;
    if (NS = 1) {
        vec v = {0,0};
        int i = 1;
        double lambda = (u_maxD + C_d * pow(v(i),2)) / (u_maxD-C_d* pow(v(i-1),2)) * exp(P);
        double v_bar = sqrt((lambda-1) / (lambda+1) * u_maxD / C_d);
        if (fabs(v_bar - v_maxD) < 1e-10)
        {
            v_bar = (1-1e-10) * v_maxD;
        }
        result.s_bar1 = 1/(2*C_d)*log((u_maxD)/(u_maxD-C_d*pow(v_bar,2)));
        result.s_bar2 = P-1/(2*C_d)*log((u_maxD+C_d*pow(v_bar,2))/(u_maxD));
        result.v = v;
        result.v_bar = v_bar;


    } else {
        //////////////////////////////////////////////////
        //  Did not translate here since no obsticle    //
        //  findPathSpeeds.m Line 23-110                //
        //////////////////////////////////////////////////
    }
    vec T(2,fill:zeros);
    for (int i = 0; i < NS; i++)
    {
        result.T_bar1 = T(0)+1/sqrt(u_maxD*C_d)*(atanh(result.v_bar/v_maxD)-atanh(result.v(i)/v_maxD));
        result.T_bar2 = result.T_bar1 + (result.s_bar2-result.s_bar1)/result.v_bar;
        T(i+1) = result.T_bar2 + 1 / sqrt(u_maxD*C_d)*(atan((result.v_bar/v_maxD))-atan((result.v(i+1)/v_maxD)));
        //////////////////////////////////////////////////
        //  Did not translate here since no obsticle    //
        //  findPathSpeeds.m Line 153-175               //
        //////////////////////////////////////////////////
    }
    result.T = T;
    result.u_maxD = u_maxD;
    result.v_maxD = v_maxD;
    result.v_maxDC = v_maxDC;

}