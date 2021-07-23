#pragma once
#include <armadillo>
#include <iostream>
#include <fstream>
#include <istream>
#include "AllParametersExperiment.cpp"
#include "lambdaInterSecCircLine.cpp"

using namespace std;
using namespace arma;

struct lineIntersec
{
    int Flag;
    mat Pos1;
    mat Pos2;
    vec S10;
    vec S20;
    double Pbar1;
    double Pbar2;
};

lineIntersec interSecLineLine(vec r1, vec r2, double mL1, double cL1, double theta1, double drx, double dry, double L1,
                              vec r21, vec r22, double mL2, double cL2, double theta2, double drx2, double dry2, double L2, double dtheta) {
    lineIntersec result;
    result.Pos1.resize(2,2);
    result.S10.resize(2);
    result.Pos2.resize(2,2);
    result.S20.resize(2);
    result.Flag = 0;
    vec rc(2);
    if (mL1 == mL2) {
        // cout << "enter if" << endl;
        double d  = fabs((cL1-cL2)/sqrt(1+pow(mL1,2)));
        if (d < 2 * rho_D)
        {
            result.Flag = 1;
            if (L1 < L2)
            {
                result.Pos1.col(0) = r1;
                result.Pos1.col(1) = r2;
                result.S10(0)=0;
                result.S10(1) = L1;
                result.Pbar1 = L1;
                result.Pos2(0,0) = (r1(0)+mL1*r1(1)-mL1*cL2)/(1+pow(mL1,2));
                result.Pos2(0,1) = (r2(0)+mL1*r2(1)-mL1*cL2)/(1+pow(mL1,2));
                if (mL2 < 1e16)
                {
                    result.Pos2(1,0) = mL2 * result.Pos2(0,0) + cL2;
                    result.Pos2(1,1) = mL2 * result.Pos2(0,1) + cL2;
                } else {
                    result.Pos2(1,0) = r1(1);
                    result.Pos2(1,1) = r2(1);
                }
                result.S20(0) = arma::norm(r21 - result.Pos2.col(0));
                result.S20(1) = result.S20(0) + L1;
                result.Pbar2 = L1;
            } else {
                result.Pos2.col(0) = r21;
                result.Pos2.col(1) = r22;
                result.S20(0) = 0;
                result.S20(1) = L2;
                result.Pbar2 = L2;
                result.Pos1(0,0) = (r21(0)+mL1*r21(1)-mL1*cL2)/(1+pow(mL1,2));
                result.Pos1(0,1) = (r22(0)+mL1*r22(1)-mL1*cL2)/(1+pow(mL1,2));
                if (mL1 < 1e16)
                {
                    result.Pos1(1,0) = mL1 * result.Pos1(0,0) + cL1;
                    result.Pos1(1,1) = mL1 * result.Pos1(0,1) + cL1;
                } else {
                    result.Pos1(1,0) = r21(1);
                    result.Pos1(1,1) = r22(2);
                }
                result.S10(0) = arma::norm(r1-result.Pos1.col(0));
                result.S10(1) = result.S10(0) + L2;
                result.Pbar1 = L2;
            }
        }
    } else {
        // cout << "enter else" << endl;
        double x_int = (cL2 - cL1) / (mL1 - mL2);
        double y_int;
        if (mL1 > 1e16)
        {
            y_int = mL2 * x_int + cL2;
        } else {
            y_int = mL1 * x_int + cL1;
        }
        double Pbar0 = max(rho_D*sqrt(2/(1-fabs(cos(dtheta)))),2*rho_D);
        double dx = Pbar0*fabs(cos(theta1));
        double dy = Pbar0*fabs(sin(theta1));
        double xbar11=x_int-arma::sign(drx)*dx;
        double ybar11=y_int-arma::sign(dry)*dy;
        double xbar12=x_int+arma::sign(drx)*dx;
        double ybar12=y_int+arma::sign(dry)*dy;

        dx=Pbar0*fabs(cos(theta2));
        dy=Pbar0*fabs(sin(theta2));
        double xbar21=x_int-arma::sign(drx2)*dx;
        double ybar21=y_int-arma::sign(dry2)*dy;
        double xbar22=x_int+arma::sign(drx2)*dx;
        double ybar22=y_int+arma::sign(dry2)*dy;
        // cout << "xbar11: " << xbar11 << endl;
        // cout << "ybar11: " << ybar11 << endl;
        // cout << "xbar12: " << xbar12 << endl;
        // cout << "ybar12: " << ybar12 << endl;
        // cout << "xbar21: " << xbar21 << endl;
        // cout << "ybar21: " << ybar21 << endl;
        // cout << "xbar22: " << xbar22 << endl;
        // cout << "ybar22: " << ybar22 << endl;
        double lambda1, lambda11, lambda12;
        double lambda2, lambda21, lambda22;
        if (drx != 0) {
            lambda1=(x_int-r1(0))/drx;//norm([x_int,y_int]'-r1)/L1;
            lambda11=(xbar11-r1(0))/drx;//norm([xbar11,ybar11]'-r1)/L1;
            lambda12=(xbar12-r1(0))/drx;//norm([xbar12,ybar12]'-r1)/L1;
        } else {
            lambda1=(y_int-r1(1))/dry;//norm([x_int,y_int]'-r1)/L1;
            lambda11=(ybar11-r1(1))/dry;//norm([xbar11,ybar11]'-r1)/L1;
            lambda12=(ybar12-r1(1))/dry;//norm([xbar12,ybar12]'-r1)/L1;
        }

        if (drx2 != 0) {
            lambda2=(x_int-r21(0))/drx2;//norm([x_int,y_int]'-r21)/L2;
            lambda21=(xbar21-r21(0))/drx2;//norm([xbar21,ybar21]'-r21)/L2;
            lambda22=(xbar22-r21(0))/drx2;//norm([xbar22,ybar22]'-r21)/L2;
        } else {
            lambda2=(y_int-r21(1))/dry2;//norm([x_int,y_int]'-r21)/L2;
            lambda21=(ybar21-r21(1))/dry2;//norm([xbar21,ybar21]'-r21)/L2;
            lambda22=(ybar22-r21(1))/dry2;//norm([xbar22,ybar22]'-r21)/L2;
        }
        if (lambda1 < 1e-6 && lambda1 > -1e-6)
        {
            lambda1 = 0;
        } else if (lambda1 < 1+1e-6 && lambda1 > 1-1e-6) {
            lambda1 = 1;
        }
        if (lambda2 < 1e-6 && lambda2 > -1e-6)
        {
            lambda2 = 0;
        } else if (lambda2 < 1+1e-6 && lambda2 > 1-1e-6) {
            lambda2 = 1;
        }
        
        // cout << "lambda1: " << lambda1 << endl;
        // cout << "lambda11: " << lambda11 << endl;
        // cout << "lambda12: " << lambda12 << endl;
        // cout << "lambda2: " << lambda2 << endl;
        // cout << "lambda21: " << lambda21 << endl;
        // cout << "lambda22: " << lambda22 << endl;
        // cout << "1-lamda1: " << 1-lambda1 << endl;
        vec xy_int = {x_int, y_int};
        double lam11,lam12,lam21,lam22;
        if (lambda2 >= 1.0)
        {
            // cout << "Process: lambda2 >=1 !!!" <<endl;
            if (lambda1 > 0.0 && lambda1 < 1.0)
            {
                // cout << "Process:0 < lambda1 < 1 !!!" <<endl;
                if (arma::norm(r22 - xy_int) < 2*rho_D/fabs(sin(dtheta)))
                {
                    result.Flag = 1;
                    lam22 = 1;
                    rc = r22;
                    vec lambda0 = lambdaInterSecCircLine(rc, 2*rho_D, r1, mL1, cL1, drx, dry);
                    mat a = r21 - xy_int;
                    a = a.t();
                    mat b = r2 - xy_int;
                    mat c = a * b;
                    if(c(0,0) > 0) {
                        lam11 = max(0.0, lambda0.min());
                        if (lambda12 > 1) {
                            lam12 = 1;
                            rc = r2;
                            vec lambda00 = lambdaInterSecCircLine(rc, 2*rho_D, r21, mL2, cL2, drx2, dry2);
                            if (lambda00.min() <=1 )
                            {
                                lam21 = max(0.0, lambda00.min());
                            } else {
                                lam21 = 1;
                                result.Flag =1;
                            }
                        } else {
                            lam12 = lambda12;
                            lam21 = max(0.0, lambda21);
                        }
                    } else {
                        lam12 = min(1.0, lambda0.min());
                        if (lambda11 < 0) {
                            lam11 = 0;
                            rc = r1;
                            vec lambda00 = lambdaInterSecCircLine(rc, 2*rho_D, r21, mL2, cL2, drx2, dry2);
                            if (lambda00.min() <= 1) {
                                lam21 = max(0.0, lambda00.min());
                            } else {
                                lam21 = 1;
                                result.Flag =1;
                            }
                        } else {
                            lam11 = lambda11;
                            lam21 = max(0.0, lambda21);
                        }
                    }
                    
                }
                
            } else if (lambda1 <= 0.0)
            {
                // cout << "Process: lambda1 <=0 !!!" <<endl;
                if (arma::norm(r22 - r1) < 2*rho_D)
                {
                    lam11 = 0;
                    lam22 = 1;
                    mat a = r21 - xy_int;
                    a = a.t();
                    mat b = r2 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0)
                    {
                        lam12 = min(1.0, lambda12);
                        lam21 = max(0.0, lambda21);
                    } else {
                        double aa = arma::norm(r22 - xy_int);
                        double bb = -aa * cos(dtheta) + sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                        lam12 = min(1.0, lambda1 + bb/L1);
                        aa = arma::norm(r1 - xy_int);
                        bb = -aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                        lam21 = max(0.0, lambda2 - bb/L2);
                    }
                    result.Flag = 1;  
                } 
                
            } else if (lambda1 >= 1.0)
            {   
                // cout << "Process: lambda1 >=1 !!!" <<endl;
                if (arma::norm(r22 - r2) < 2 * rho_D)
                {
                    lam12 = 1;
                    lam22 = 1;
                    mat a = r21 - xy_int;
                    a = a.t();
                    mat b = r1 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0)
                    {
                        lam11 = max(lambda11, 0.0);
                        lam21 = max(lambda21, 0.0);
                    } else {
                        double aa = arma::norm(r22-xy_int);
                        double bb = -aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                        lam11 = max(1.0, lambda1 - bb/L1);
                        aa = arma::norm(r2 - xy_int);
                        bb = -aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                        lam21 = max(0.0, lambda2 - bb/L2);
                    }
                    result.Flag = 1;
                }
                
            }
            
            
            
        } else if (lambda2 > 0.0 && lambda2 < 1.0)
        {
            // cout << "Process: 0 < lambda2 < 1  !!!" <<endl;
            if (lambda1 > 0.0 && lambda1 < 1.0)
            {
                // cout << "Process: 0 < lambda1 < 1 !!!" <<endl;
                if (lambda11 < 0.0)
                {
                    double aa = arma::norm(r1 - xy_int);
                    double bb = aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                    mat a = r1 - xy_int;
                    a = a.t();
                    mat b = r21 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0) {
                        lambda21 = max(lambda21, lambda2 - bb/L2);
                    } else {
                        lambda22 = min(lambda22, lambda2 + bb/L2);
                    }
                }
                if(lambda12 > 1.0) {
                    double aa = arma::norm(r2 - xy_int);
                    double bb = aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                    mat a = r2 - xy_int;
                    a = a.t();
                    mat b = r21 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0.0) {
                        lambda21 = max(lambda21, lambda2 - bb/L2);
                    } else {
                        lambda22 = min(lambda22, lambda2 + bb/L2);
                    }
                }
                if(lambda21 < 0.0) {
                    double aa = arma::norm(r21 - xy_int);
                    double bb = aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                    mat a = r21 - xy_int;
                    a = a.t();
                    mat b = r1 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0) {
                        lambda11 = max(lambda11, lambda1 - bb/L1);
                    } else {
                        lambda12 = min(lambda12, lambda1 + bb/L1);
                    }
                }
                if(lambda22 > 1.0) {
                    double aa = arma::norm(r22 - xy_int);
                    double bb = aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                    mat a = r22 - xy_int;
                    a = a.t();
                    mat b = r1 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0) {
                        lambda11 = max(lambda11, lambda1 - bb/L1);
                    } else {
                        lambda12 = min(lambda12, lambda1 + bb/L1);
                    }
                }
                lam11 = max(lambda11,0.0);
                lam12 = min(1.0, lambda12);
                lam21 = max(lambda21,0.0);
                lam22 = min(1.0, lambda22);
                result.Flag = 1;
            } else if (lambda1 <= 0.0)
            {
                // cout << "Process: lambda1 <= 0 !!!" <<endl;
                if (arma::norm(r1 - xy_int) < 2 * rho_D/fabs(sin(dtheta)))
                {
                    result.Flag = 1;
                    lam11 = 0;
                    rc = r1;
                    vec lambda0 = lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2);
                    mat a = r22 - xy_int;
                    a = a.t();
                    mat b = r2 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0)
                    {
                        lam21 = max(0.0, lambda0.min());
                        if (lambda22 > 1)
                        {
                            lam22 = 1;
                            rc = r22;
                            vec lambda00 = lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry);
                            if(lambda00.max() >= 0) 
                            {
                                lam12 = min(1.0, lambda00.max());
                            } else {
                                lam12 = 0;
                                result.Flag = 0;
                            }
                        } else {
                            lam22 = lambda22;
                            lam12 = min(1.0, lambda12);
                        }
                        
                    } else 
                    {
                        lam22 = min(1.0, lambda0.max());
                        if (lambda21 < 0)
                        {
                            lam21 = 0;
                            rc = r21;
                            vec lambda00 = lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry);
                            if (lambda00.min() >= 0)
                            {
                                lam12 = min(1.0, lambda00.max());
                            } else {
                                lam12 = 0;
                                result.Flag = 0;
                            }
                            
                        } else {
                            lam21 = lambda21;
                            lam12 = min(1.0, lambda12);
                        }

                    }
                    
                }
                
            } else if (lambda1 >=1.0)
            {   
                // cout << "Process: lambda1 >= 1 !!!" <<endl;
                if (arma::norm(r2 - xy_int) < 2 * rho_D/fabs(sin(dtheta)))
                {
                    result.Flag = 1;
                    lam12 = 1;
                    rc = r2;
                    vec lambda0 = lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2);
                    mat a = r22 - xy_int;
                    a = a.t();
                    mat b = r1 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0)
                    {
                        lam21 = max(0.0, lambda0.min());
                        if (lambda22 > 1)
                        {
                            lam22 = 1;
                            rc = r22;
                            vec lambda00 = lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry);
                            if (lambda00.min() <= 1)
                            {
                                lam11 = max(0.0, lambda00.min());
                            } else {
                                lam11 = 1;
                                result.Flag = 0;
                            }
                            
                        } else {
                            lam22 = lambda22;
                            lam11 = max(0.0, lambda11);
                        }
                        
                    } else {
                        lam22 = min(1.0, lambda0.max());
                        if (lambda21 < 0) 
                        {
                            lam21 = 0;
                            rc = r21;
                            vec lambda00 = lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry);
                            if (lambda00.min() <= 1) {
                                lam11 = max(0.0, lambda00.min());
                            } else {
                                lam11 = 1;
                                result.Flag = 0;
                            }
                        } else {
                            lam21 = lambda21;
                            lam11 = max(0.0, lambda11);
                        }
                    }
                    
                }
                
            }
            
            
            
        } else { //lambda2 < 0
            // cout << "Process: lambda2 <= 0 !!!" <<endl;
            if (lambda1 > 0.0 && lambda1 < 1.0)
            {
                // cout << "Process: 0 < lambda1 < 1 !!!" <<endl;
                if (arma::norm(r21 - xy_int) < 2*rho_D/fabs(sin(dtheta))) {
                    result.Flag = 1;
                    lam21 = 0;
                    rc = r21;
                    vec lambda0 = lambdaInterSecCircLine(rc,2*rho_D,r1,mL1,cL1,drx,dry);
                    mat a = r22 - xy_int;
                    a = a.t();
                    mat b = r2 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0)
                    {
                        lam11 = max(0.0, min(lambda0));
                        if (lambda12 > 1)
                        {
                            lam12 = 1;
                            rc = r2;
                            vec lambda00 = lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2);
                            if (lambda00.max() >=0 )
                            {
                                lam22 = min(1.0, lambda00.max());
                            } else {
                                lam22 = 0;
                                result.Flag = 0;
                            }
                        } else {
                            lam12 = lambda12;
                            lam22 = min(1.0, lambda22);
                        }
                        
                    } else {
                        lam12 = min(1.0, lambda0.max());
                        if (lambda11 < 0) {
                            lam11 = 0;
                            rc = r1;
                            vec lambda00 = lambdaInterSecCircLine(rc,2*rho_D,r21,mL2,cL2,drx2,dry2);
                            if (lambda00.max() >= 0) {
                                lam22 = min(1.0, lambda00.max());
                            } else {
                                lam22 = 0;
                                result.Flag = 0;
                            }
                        } else {
                            lam11 = lambda11;
                            lam22 = min(1.0,lambda22);
                        }
                    }
                    
                }
            } else if (lambda1 <= 0.0)
            {
                // cout << "Process: lambda1 <= 0 !!!" <<endl;
                if (arma::norm(r21 - r1) < 2 * rho_D)
                {
                    lam21 = 0;
                    lam11 = 0;
                    mat a = r22 - xy_int;
                    a = a.t();
                    mat b = r1 - xy_int;
                    mat c = a * b;
                    if (c(0,0) > 0)
                    {
                        lam12 = min(1.0, lambda12);
                        lam22 = min(1.0, lambda22);
                    } else {
                        double aa = arma::norm(r1 - xy_int);
                        double bb = -aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                        lam22 = min(1.0,lambda2 + bb/L2);
                        aa = arma::norm(r21 - xy_int);
                        bb = -aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                        lam12 = min(1.0, lambda1 + bb/L1);
                    }
                    result.Flag = 1;
                } else if (lambda1 >= 1.0)
                {
                    // cout << "Process: lambda1 >= 1 !!!" <<endl;
                    if (arma::norm(r2 - r21) < 2 * rho_D)
                    {
                        lam21 = 0;
                        lam12 = 1;
                        mat a = r22 - xy_int;
                        a = a.t();
                        mat b = r1 - xy_int;
                        mat c = a * b;
                        if (c(0,0) > 0)
                        {
                            lam11 = max(0.0, lambda11);
                            lam22 = min(1.0, lambda22);   
                        } else {
                            double aa = arma::norm(r2 - xy_int);
                            double bb = -aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                            lam22 = min(1.0, lambda2 + bb/L2);
                            aa = arma::norm(r21 - xy_int);
                            bb = -aa*cos(dtheta)+sqrt(4*pow((rho_D),2)-pow(aa,2)*pow(sin(dtheta),2));
                            lam11 = max(0.0, lambda1 - bb/L1);
                        }
                        result.Flag = 1;
                    }
                    
                }
                
                
            }
            
            
        }

        if (result.Flag == 1)
        {
            result.Pos1.col(0) = {r1(0)+lam11*drx, r1(1)+lam11*dry};
            result.Pos1.col(1) = {r1(0)+lam12*drx, r1(1)+lam12*dry};
            result.Pbar1 = (lam12-lam11)*L1;
            result.S10(0) = lam11*L1;
            result.S10(1) = lam12*L1;

            result.Pos2.col(0) = {r21(0)+lam21*drx2, r21(1)+lam21*dry2};
            result.Pos2.col(1) = {r21(0)+lam22*drx2, r21(1)+lam22*dry2};
            result.Pbar2 = (lam22-lam21)*L2;
            result.S20(0) = lam21*L2;
            result.S20(1) = lam22*L2;

        }
        
        
        
    }
    // cout << "in interSecLineLine function Pbar1: " << result.Pbar1 << endl;
    // cout << "in interSecLineLine function Pbar2: " << result.Pbar2 << endl;
    return result;
}


// int main() {
//     // vec r1 = {11.9907,-10.7580};
//     // vec r2 = {14.3142742593192,
//     //         -18.2808998346075};
//     // double mL1 = -3.23764123502204;
//     // double mL2 = -2.25685625521036;
//     // double cL1 = 28.0635847567788;
//     // double cL2 = 14.0243595663336;
//     // double theta1 = 5.01196050370651;
//     // double theta2 = 5.12948525219316;
//     // double drx = 2.32357425931916;
//     // double dry = -7.52289983460749;
//     // double drx2 = 4.94187425931916;
//     // double dry2 = -11.1530998346075;
//     // double dtheta = 0.117524748486650;
//     // double L1 = 7.87356458410725;
//     // double L2 = 12.1989244245402;
//     // vec r21 = {9.37240000000000,
//     //             -7.12780000000000};
//     // vec r22 = {14.3142742593192,
//     //             -18.2808998346075};

//     vec r1, r2, r21, r22;
//     r1.load("../../../../../Downloads/swarm_matlab/LineV/r1.txt");
//     r2.load("../../../../../Downloads/swarm_matlab/LineV/r2.txt");
//     r21.load("../../../../../Downloads/swarm_matlab/LineV/r21.txt");
//     r22.load("../../../../../Downloads/swarm_matlab/LineV/r22.txt");

//     std::ifstream fin("../../../../../Downloads/swarm_matlab/LineV/cL1.txt");
//     double cL1;
//     fin >> cL1;
//     std::cout << "cL1: "<<cL1 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/cL2.txt");
//     double cL2;
//     fin >> cL2;
//     std::cout << "cL2: "<<cL2 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/drx.txt");
//     double drx;
//     fin >> drx;
//     std::cout << "drx: "<< drx << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/drx2.txt");
//     double drx2;
//     fin >> drx2;
//     std::cout << "drx2: "<<drx2 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/dry.txt");
//     double dry;
//     fin >> dry;
//     std::cout << "dry: "<<dry << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/dry2.txt");
//     double dry2;
//     fin >> dry2;
//     std::cout << "dry2: "<<dry2 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/dtheta.txt");
//     double dtheta;
//     fin >> dtheta;
//     std::cout << "dtheata: "<<dtheta << std::endl;
//     fin.close();
        
//     fin.open("../../../../../Downloads/swarm_matlab/LineV/L1.txt");
//     double L1;
//     fin >> L1;
//     std::cout << "L1: "<<L1 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/L2.txt");
//     double L2;
//     fin >> L2;
//     std::cout << "L2: "<<L2 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/mL1.txt");
//     double mL1;
//     fin >> mL1;
//     std::cout << "mL1: "<<mL1 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/mL2.txt");
//     double mL2;
//     fin >> mL2;
//     std::cout << "mL2: "<<mL2 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/theta1.txt");
//     double theta1;
//     fin >> theta1;
//     std::cout << "theta1: "<<theta1 << std::endl;
//     fin.close();

//     fin.open("../../../../../Downloads/swarm_matlab/LineV/theta2.txt");
//     double theta2;
//     fin >> theta2;
//     std::cout << "theta2: "<<theta2 << std::endl;
//     fin.close();    



//     lineIntersec a = interSecLineLine(r1,r2, mL1,cL1, theta1, drx, dry, L1,
//                                       r21, r22, mL2, cL2, theta2, drx2, dry2, L2, dtheta);
//     cout << "Flag: " << a.Flag << endl;
//     cout << "Pbar1: " <<a.Pbar1 << endl;
//     cout << "Pbar2: " <<a.Pbar2 << endl;
//     a.Pos1.print("Pos1: ");
//     a.Pos2.print("Pos2: ");
//     a.S10.print("S10: ");
//     a.S20.print("S20: ");
// }