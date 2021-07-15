#include "iostream"
#include <math.h>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex.h>
#include "symDerivative.hpp"

double Del,p,p_dot,p_ddot,pp,pp_dot,g,g_dot,g_ddot,FAO,FAO_dot,FAO_ddot;

double Delta_barS=1;
double F_bar_AO=0.9;

double dfdp(double FAO, double g, double p, double Del, double pp) {
    return  -FAO*cos(g)*(pow(tan(p),2.0)+1.0)-Del*cos(pp)*(pow(tan(p),2.0)+1.0);
}

double dfdpp(double p, double Del, double pp) {
    return Del*(cos(pp)+sin(pp)*tan(p));
}

double dfdg(double p, double FAO, double g) {
    return FAO*(cos(g)+sin(g)*tan(p));
}

double dfdFAO(double g, double p) {
    return sin(g)-cos(g)*tan(p);
}

double d2fdp2(double FAO, double g, double p, double Del, double pp) {
    return FAO*cos(g)*tan(p)*(pow(tan(p),2.0)+1.0)*-2.0-Del*cos(pp)*tan(p)*(pow(tan(p),2.0)+1.0)*2.0;
}

double d2fdpp2(double p, double Del, double pp) {
    return -Del*(sin(pp)-cos(pp)*tan(p));
}

double d2fdg2(double p, double FAO, double g) {
    return 0;
}

double psi_prime_dot0(double FAO, double g, double p, double Del, double pp) {
    return -(-p_dot*(FAO*cos(g)*(pow(tan(p),2.0)+1.0)+Del*cos(pp)*(pow(tan(p),2.0)+1.0))
           +FAO_dot*(sin(g)-cos(g)*tan(p))+FAO*g_dot*(cos(g)+sin(g)*tan(p)))
           /(Del*(cos(pp)+sin(pp)*tan(p)));
}

double psi_prime_ddot0(double FAO, double g, double p, double Del, double pp) {
    return (p_ddot*(FAO*cos(g)*(pow(tan(p),2.0)+1.0)+Del*cos(pp)*(pow(tan(p),2.0)+1.0))
    +(p_dot*p_dot)*(FAO*cos(g)*tan(p)*(pow(tan(p),2.0)+1.0)*2.0+Del*cos(pp)*tan(p)*
    (pow(tan(p),2.0)+1.0)*2.0)-FAO_ddot*(sin(g)-cos(g)*tan(p))-FAO*g_ddot*(cos(g)+sin(g)*tan(p))
    +FAO*(g_dot*g_dot)*(sin(g)-cos(g)*tan(p))+Del*(pp_dot*pp_dot)*(sin(pp)-cos(pp)*tan(p)))/(Del*(cos(pp)+sin(pp)*tan(p)));
}

double sig[1][2];
double sig_dot[1][2];
double sig_ddot[1][2];
void sig_generation(double sig1, double sig2,
                    double sig_dot1, double sig_dot2,
                    double sig_ddot1, double sig_ddot2) {
    sig[0][0] = sig1;
    sig[0][1] = sig2;
    sig_dot[0][0] = sig_dot1;
    sig_dot[0][1] = sig_dot2;
    sig_ddot[0][0] = sig_ddot1;
    sig_ddot[0][1] = sig_ddot2;
}

double sigProd(double sig[1][2]) {
    return(sig[0][0] - 1) * (sig[0][1] - 1);
}

double dsigProd[1][2];
double d2sigProd[1][2];
void sigProd_generation(double sig1, double sig2){
    dsigProd[0][0] = sig2-1.0;
    dsigProd[0][1] = sig1-1.0;
    d2sigProd[0][0] = 0;
    d2sigProd[0][1] = 0;
}

double sigProdA_dot(double sigA_dot1, double sigA_dot2, double sigA_dot3,
                    double sigA_dot4, double sigA_dot5, double sigA_dot6,
                    double sigA1, double sigA2, double sigA3,
                    double sigA4, double sigA5, double sigA6) {
    return sigA_dot1*(sigA2-1.0)*(sigA3-1.0)*(sigA4-1.0)*(sigA5-1.0)*(sigA6-1.0)
    +sigA_dot2*(sigA1-1.0)*(sigA3-1.0)*(sigA4-1.0)*(sigA5-1.0)*(sigA6-1.0)
    +sigA_dot3*(sigA1-1.0)*(sigA2-1.0)*(sigA4-1.0)*(sigA5-1.0)*(sigA6-1.0)
    +sigA_dot4*(sigA1-1.0)*(sigA2-1.0)*(sigA3-1.0)*(sigA5-1.0)*(sigA6-1.0)
    +sigA_dot5*(sigA1-1.0)*(sigA2-1.0)*(sigA3-1.0)*(sigA4-1.0)*(sigA6-1.0)
    +sigA_dot6*(sigA1-1.0)*(sigA2-1.0)*(sigA3-1.0)*(sigA4-1.0)*(sigA5-1.0);
}

double sigProdD_dot(double sigD_dot1, double sigD_dot2, double sigD_dot3, double sigD_dot4,
                    double sigD1, double sigD2, double sigD3, double sigD4) {
    return sigD_dot1*(sigD2-1.0)*(sigD3-1.0)*(sigD4-1.0)
    +sigD_dot2*(sigD1-1.0)*(sigD3-1.0)*(sigD4-1.0)
    +sigD_dot3*(sigD1-1.0)*(sigD2-1.0)*(sigD4-1.0)
    +sigD_dot4*(sigD1-1.0)*(sigD2-1.0)*(sigD3-1.0);
}

// double d2Beta_barF_dBeta2(double beta, double a0, double b0, double n0) {
//     double t2 = cos(beta);
//     double t3 = sin(beta);
//     double t4 = n0*2.0;
//     double t5 = t2*t2;
//     double t6 = t3*t3;
//     double t7 = pow(a0,t4);
//     double t8 = pow(b0,t4);
//     double t9 = t4-1.0;
//     double t10 = t9-1.0;
//     double t12 = pow(t2,t9);
//     double t15 = pow(t3,t9);
//     double t11 = t10-1.0;
//     double t13 = pow(t2,t10);
//     double t16 = pow(t3,t10);
//     double t18 = t8*t12;
//     double t19 = t7*t15;
//     double t14 = pow(t2,t11);
//     double t17 = pow(t3,t11);
//     double t20 = cimag(t18);
//     double t21 = creal(t18);
//     double t22 = cimag(t19);
//     double t23 = creal(t19);
//     double t25 = t2*t8*t9*t13;
//     double t26 = t2*t7*t9*t16;
//     double t27 = t3*t8*t9*t13;
//     double t28 = t3*t7*t9*t16;
//     double t24 = -t22;
//     double t29 = cimag(t26);
//     double t30 = cimag(t27);
//     double t31 = creal(t26);
//     double t32 = creal(t27);
//     double t34 = t5*t7*t9*t10*t17;
//     double t35 = t6*t8*t9*t10*t14;
//     double t36 = t20+t23;
//     double t33 = -t31;
//     double t37 = t36*t36;
//     double t38 = t21+t24;
//     double t39 = 1.0/t36;
//     double t42 = t29+t32;
//     double t40 = 1.0/t37;
//     double t41 = t38*t38;
//     double t43 = t30+t33;
//     double t46 = t39*t42;
//     double t44 = t37+t41;
//     double t47 = t38*t40*t43;
//     double t45 = 1.0/t44;
//     double t48 = -t47;
//     double t49 = t46+t48;
//     double t0 = -t37*t45*(t39*(cimag(t28)-cimag(t34)-creal(t25)+ creal(t35))
//                 +t38*t40*(cimag(t25)-cimag(t35)+creal(t28)-creal(t34))+t38*(t39*t39*t39)*(t43*t43)*2.0-t40*t42*t43*2.0)
//                 +t37*(t45*t45)*t49*(t36*t43*2.0+t38*t42*2.0)-t36*t43*t45*t49*2.0;
//     return t0;
// }