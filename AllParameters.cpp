#include "iostream"
#include "AllParameters.hpp" 
#include <math.h>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
 
// Driver Code
double dt = 0.01;
double Rs = 1000;
double rho_P = 45;

int rP[2] = {0, 200};
int rS[2] = {500, 750};

double rho_A = 0.5;
double rho_c_A = 150;

double C_d = 0.2;
double vma = 6;
double uma;
double rhoA_safe;

void calControlLimits() {
    if (C_d != 0)
    {
        uma = C_d * vma * vma; //temporary safe distances (not based on the quadratic drag term
        //analysis
        rhoA_safe = log((uma+ C_d * pow((2*vma),2))/uma)/ (2*C_d); // additional aafety distance required between the attackers due to bounded acceleration
    } else {
        uma = 0.35 * vma;
        rhoA_safe = pow(vma+vma,2) / (2*uma);
    }
}

double R_bar_AA1;
double R_bar_AA2;
double R_m_AA;
double R_bar_AA;
double R_u_AA;
double dR_AA_cube;
double A_A_A;
double B_A_A;
double C_A_A;
double D_A_A;

void VectorFields_A() {
    R_bar_AA1=18*rho_A;
    R_bar_AA2=25*rho_A;
    R_m_AA=1*(rho_A+rho_A)+rhoA_safe+1;
    R_bar_AA=R_m_AA+R_bar_AA1+rhoA_safe;
    R_u_AA=R_m_AA+R_bar_AA2+rhoA_safe;
    dR_AA_cube=pow((R_u_AA-R_bar_AA),3);
    A_A_A=2/dR_AA_cube;
    B_A_A=-3*(R_u_AA+R_bar_AA)/dR_AA_cube;
    C_A_A=6*R_u_AA*R_bar_AA/dR_AA_cube;
    D_A_A=pow(R_u_AA,2)*(R_u_AA-3*R_bar_AA)/dR_AA_cube;
}

double v_maxA[NA], u_maxA[NA];

void initialAttackersVel(){
    for (size_t i = 0; i < NA; i++)
    {
       v_maxA[i] = vma;
       u_maxA[i] = uma;
    }
}

double NA_sep = 0;

int axS[2] = {1, 1};
// rDmin = rho_P+10
int rDmin = 7;


double RA0;
double potential;
double Rii00[NA], Rik00[NA];
std::vector<double> Rjk00;
void cal(){
    RA0 = R_m_AA + 4;
    potential = RA0 * sqrt(2*(1-cos(2*M_PI/NA))); //for formation potential
    for (size_t i = 0; i < NA; i++)
    {
        Rii00[i] = potential;
        Rik00[i] = 5 * Rii00[i];
    }
}

int dRA0 = 10;

int sceneraio = 3;
double RA01 = RA0;
double rA0[2],rA[2][NA],vA[2][NA],thetaA0;

double XA0[4][NA];
double rA_follow[2][NA];
double RA;
double rho_S;
double rho_Acon;

void InitializeAttackers(int s) {
    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < NA; j++)
        {
            rA[i][j] = 0;
            vA[i][j] = 0;
        }
        
    }
    
    if (s == 1) {
        if (0) {
            double rA0[2] = {-955, 710};
            rA[0][0] =rA0[0];
            rA[1][0] = rA0[1];
            thetaA0=atan2(rA0[2]-rP[2],rA0[1]-rP[1]);
            cout << "1: thetaA0: "<<thetaA0 << endl;
            rA[0][1] = rA0[0] + RA01 * cos(thetaA0-M_PI/10);
            rA[1][1] = rA0[1] + RA01 * sin(thetaA0-M_PI/10);
            rA[0][2] = rA0[0] + RA01 * cos(thetaA0+2*M_PI/10);
            rA[1][2] = rA0[1] + RA01 * sin(thetaA0+2*M_PI/10);
            for (size_t i = 3; i < NA; i++)
            {
                double thetaA = thetaA0;               //////////////////////////////////////////////////////////////////////////////
                RA=sqrt(2)*(RA0)+RA0*(i-3);     //  + 1*randn; did not figure out normal distributed random number in C++   //
                rA[0][i] = rA0[0] + RA * cos(thetaA);  //////////////////////////////////////////////////////////////////////////////
                rA[1][i] = rA0[1] + RA * sin(thetaA);
                std::complex<double> vector (rA[0][i] - rP[0], rA[1][i] - rP[1]);
                double calculatedNorm = sqrt(std::norm(vector));
                vA[0][i] = (rA[0][i] - rP[0]) * 0.01 * v_maxA[i] / calculatedNorm;
                vA[1][i] = (rA[1][i] - rP[1]) * 0.01 * v_maxA[i] / calculatedNorm;
            }
        } else {
            double rA0[2]= {-845,80};
            //rA0=[-875,200]';
            thetaA0=atan2(rA0[1]-rP[1],rA0[0]-rP[0])+M_PI/3;
            for (size_t i = 0; i < NA; i++)
            {
                double thetaA = thetaA0;
                double RA = RA0 * (i);
                rA[0][i] = rA0[0] + RA*cos(thetaA);
                rA[1][i] = rA0[1] + RA*sin(thetaA);
                std::complex<double> vector (rA[0][i] - rP[0], rA[1][i] - rP[1]);
                double calculatedNorm = sqrt(std::norm(vector));
                vA[0][i] = -(rA[0][i] - rP[0]) * 0.01 * v_maxA[i] / calculatedNorm;
                vA[1][i] = -(rA[1][i] - rP[1]) * 0.01 * v_maxA[i] / calculatedNorm;
            }
            
        }
    } else if (s == 2) {
        double rA0[2] = {-645,820};
        thetaA0 = atan2(rA0[1]-rP[1],rA0[0]-rP[0])-M_PI/3;
        for (size_t i = 0; i < NA; i++)
        {
            double thetaA=thetaA0-M_PI/5*(i+1);
            RA=RA0/1.5*(i);// +1*randn
            rA[0][i] = rA0[0] + RA*cos(thetaA);
            rA[1][i] = rA0[1] + RA*sin(thetaA);
            std::complex<double> vector (rA[0][i] - rP[0], rA[1][i] - rP[1]);
            double calculatedNorm = sqrt(std::norm(vector));
            vA[0][i] = -(rA[0][i] - rP[0]) * 0.01 * v_maxA[i] / calculatedNorm;
            vA[1][i] = -(rA[1][i] - rP[1]) * 0.01 * v_maxA[i] / calculatedNorm;
        }
    } else if (s == 3) {
        double rA0[2] = {-740,-620};
        // rA0=[-875,200]';
        double thetaA0=atan2(rA0[1]-rP[1],rA0[0]-rP[0])-M_PI/3;
        for (size_t i = 0; i < NA; i++)
        {
            double thetaA=thetaA0+M_PI/5*(i+1);
            RA= RA0 / 1.5* (i) ;
            rA[0][i] = rA0[0] + RA * cos(thetaA);
            rA[1][i] = rA0[1] + RA * sin(thetaA);
            std::complex<double> vector (rA[0][i] - rP[0], rA[1][i] - rP[1]);
            double calculatedNorm = sqrt(std::norm(vector));
            vA[0][i] = -(rA[0][i] - rP[0]) * 0.01 * v_maxA[i] / calculatedNorm;
            vA[1][i] = -(rA[1][i] - rP[1]) * 0.01 * v_maxA[i] / calculatedNorm;
        }
    }

    for (size_t i = 0; i < NA; i++)
    {
        XA0[0][i] = rA[0][i];
        XA0[1][i] = rA[1][i];
        XA0[2][i] = vA[0][i];
        XA0[3][i] = vA[1][i];
    }
    
    rA_follow[0][0] = 0;
    rA_follow[1][0] = 0;
    for (size_t i = 1; i < NA; i++)
    {
        rA_follow[0][i] = rA[0][i] - rA[0][0];
        rA_follow[1][i] = rA[1][i] - rA[1][0];
    }
    
    //XA0(:,NA)=[-700,0,0,0]';
    double rAcm0[2] = {0, 0};
    double vAcm0[2] = {0, 0};
    
    for (size_t i = 0; i < NA; i++)
    {
        rAcm0[0] += XA0[0][i];
        rAcm0[1] += XA0[1][i];
        vAcm0[0] += XA0[2][i];
        vAcm0[1] += XA0[3][i];
    }
    rAcm0[0] = rAcm0[0] / NA;
    rAcm0[1] = rAcm0[1] / NA;
    vAcm0[0] = vAcm0[0] / NA;
    vAcm0[1] = vAcm0[1] / NA;
    
    double thetaAcm0 = atan2(rAcm0[1]-rP[1],rAcm0[0]-rP[0]);

    rho_Acon = 3*RA0+rhoA_safe+dRA0+10;  //radius of the connectivity region
    
    rho_S= 2*rho_Acon;   //radius of the safe area

}
// Defenders

double rho_D = 0.5;
double rho_sn_max = 65;
double ND;
double N;
double rho_sn;

std::vector<double> v_maxD,v_maxDC,u_maxD,u_maxD1,u_maxD2,u_maxDr1,u_maxDr2;
double rho_safe,dthetai,alphaD_v;
double rhoAD_safe,rhoD_safe;
double umd1,umd2,umdf_s1,umdf_s2,umd,umdf_h1,vmd,vmdc,vmdf_s,umd_e1,umd_e2;
obstacle obs;
void defenders() {
    ND = ceil(M_PI/acos(rho_Acon/rho_sn_max));
    N = NA + ND;
    alphaD_v=0.5;
    umdf_h1=uma;
    double kDv=.1;
    double alphaD_r=alphaD_v/(2-alphaD_v);
    umdf_s2=0.6*uma;
    umdf_s1=umdf_s2;
    double eta_max=vma+2;
    umd2=5;
    umd1=1.5*umd2;
    umd=umd1+umd2+umdf_s1+umdf_s2+1;
    umd_e1=0.5*umd;
    umd_e2=0.5*umd;
    if (C_d != 0) {
        vmd=sqrt(umd/C_d);
        rhoD_safe=log((umd + C_d*pow((2*vmd),2))/umd)/(2*C_d);   //additional aafety distance required between the defenders due to bounded acceleration
        rhoAD_safe=log((uma +C_d*pow((vma+vmd),2))/uma)/(2*C_d);  //additional safty for the attackers to stay safe from the defenders due to bounded acceleration
        rho_safe = rhoD_safe;
        vmdc=sqrt((sqrt(pow(umd,2)/(1/pow(rho_safe,2)+pow(C_d,2)))));
    } else {
        u_maxD.push_back(0.4*vmd);    // Ask here: what is the purpose of doing this since it will be assign value below
        rhoD_safe=pow((vmd+vmd),2)/(2*umd);
        rhoAD_safe=pow((vma+vmd),2)/(2*uma);
        rho_safe = rhoD_safe;
        vmdc=sqrt(umd*rhoD_safe);
    }
    vmdf_s=sqrt((umdf_s1+umdf_s2)/C_d);
    for (size_t i = 0; i < ND; i++)
    {
        v_maxD.push_back(vmd);
        v_maxDC.push_back(vmdc);
        u_maxD.push_back(umd);
        u_maxD1.push_back(umd1);
        u_maxD2.push_back(umd2);
        u_maxDr1.push_back(umdf_s1);
        u_maxDr2.push_back(umdf_s2);
    }
    obs.rho_safe = rho_safe;
    dthetai=acos(1-pow((2*rho_D),2)/(2*pow(rho_safe,2)));   //angular shift for two agents colliding on a ciruclar arc segment
    rho_sn=rho_Acon+rhoD_safe+10;   //radius of the stringNet
    for (size_t i = 0; i < ND; i++)
    {
        Rjk00.push_back(Rik00[1]);
    }
}

int largeP = 100000;

////////////////////////////////////
//  RD_con fsolve problem!        //
////////////////////////////////////

double kr1=0.5;
double kv1=0.1;
double kv2=0.5;

//control gains for attackers
double kAFr=10;
double kAOr=1.5;
double kAOr2=0.5;
double kADr=.5;  //kADr=1;
double kADv=.5;  //kADv=2;
double kAFv=.1;
double alphaAFv=1;
double kAOv=2;
double kAOv2=.2;
double alphaAOv=0.5;
double alphaADv=0.5;
double kAPr=0.003;
double kAPv=0.1;

//finite time gains for the defenders
double kDOr=1;
double kDDr=.5;
double kDOv1=2;
double kDOv2=0.08;
double kDFr=5.5;
double kDFv=.45;
double alphaDFv=.9;
double alphaDFr;
double kDFr2;
double alphaDOv=1;
double kDDv=0.3;
double alphaDDv=0.5;
double kDRr=.0001;
double kDRv=0.1;
double kDFphi=0.003;
double kDFphid=0.007;
double kDFphir=0.02;
double kDFphiv=0.1;
double kDDesr=1.5;
double kDDesv=1;

//////////////////////////////////////////////////////////////////////
////      %bound on the convergence error during tracking         ////
////      A_tilde=[zeros(2),eye(2);-kDFr*eye(2),-kDFv*eye(2)];    ////
////      Qd=1*eye(4);                                            ////
////      Pd=lyap(A_tilde',Qd);                                   ////
////      c1=min(eig(Pd));                                        ////
////      c2=max(eig(Pd));                                        ////
////      c3=min(eig(Qd));                                        ////
////      c4=2*max(eig(Pd));                                      ////
////      p=c4/c3*sqrt(c2/c1);                                    ////
////      bd=15;                                                  ////
////   Problem: lyap conuld not use coder to covert,              ////
////            and did not see p use anywhere else, it is not    ////
////            define as global                                  ////
//////////////////////////////////////////////////////////////////////


///////////////////////////////////
//Convex polygonal Obstalces///////
///////////////////////////////////
std::vector<std::vector<std::vector<double>>> rVO = {
        {{100,200,200,100},
        {410,410,550,550}},

        {{-350,-240,-240,-350},
        {450,450,660,660}},

        {{-230,-100,-100,-230},
        {-160,-160,-60,-60}},

        {{-500,-300,-400,-500},
        {-400,-400,-550,-550}},

        {{-750,-520,-520,-620,-680},
        {0,0,90,150,90}},

        {{200,450,400,250},
        {-50,-50,160,160}}

      };


double NO = rVO.size();
double NVOk;
double rOc[2] = {0,0};
std::vector<double> dVO;

std::vector<double> R_bar_AcOc, R_u_AcOc, A_Ac_Oc, B_Ac_Oc, C_Ac_Oc, D_Ac_Oc;
void FormationOri() {
    alphaDFr=alphaDFv/(2-alphaDFv);
    kDFr2=(umd2+C_d*pow(vmd,2))/pow((vmd+vma),(alphaDFv));
    for (size_t k = 0; k < NO; k++)
    {
        NVOk = rVO[k][0].size();
        for (size_t i = 0; i < rVO[k][0].size();i++) {
            rOc[0] += rVO[k][0][i];
            rOc[1] += rVO[k][1][i];
        }
        rOc[0] = rOc[0] / NVOk;
        rOc[1] = rOc[1] / NVOk;
        for (size_t j = 0; j < NVOk; j++) {
            double a = rOc[0] - rVO[k][0][j];
            double b = rOc[1] - rVO[k][1][j];
            std::complex<double> comp = {a,b};
            dVO.push_back(sqrt(std::norm(comp)));
        }
        
        std::vector<double>::iterator result = std::max_element(dVO.begin(), dVO.end());
        double max_dVO = result[0];
        R_bar_AcOc.push_back(rho_Acon + max_dVO + 5);
        R_u_AcOc.push_back(rho_Acon + max_dVO +155);
        double dR_AcOc_cube=pow((R_u_AcOc[k]-R_bar_AcOc[k]),3);
        A_Ac_Oc.push_back(2/dR_AcOc_cube);
        B_Ac_Oc.push_back(-3*(R_u_AcOc[k]+R_bar_AcOc[k])/dR_AcOc_cube);
        C_Ac_Oc.push_back(6*R_u_AcOc[k]*R_bar_AcOc[k]/dR_AcOc_cube);
        D_Ac_Oc.push_back(pow(R_u_AcOc[k],2)*(R_u_AcOc[k]-3*R_bar_AcOc[k])/dR_AcOc_cube);

        rOc[0] = 0;
        rOc[1] = 0;
        if (k != NO - 1) {
        dVO.clear();
        }
    }


}
double R_bar_AD1;
double R_bar_AD2;  //75
double R_m_AD;
//R_m_AD=0;   For testing purposes
double R_bar_AD;
double R_u_AD;
double dR_AD_cube;
double A_A_D;
double B_A_D;
double C_A_D;
double D_A_D;

double R_m_AD2; //1*(rho_A+rho_D)+rhoAD_safe;
//R_m_AD=0;  % For testing purposes
double R_u_AD2;
double dR_AD_cube2;
double A_A_D2;
double B_A_D2;
double C_A_D2;
double D_A_D2;

void calVfieldA() {
    R_bar_AD1=95*rho_A;  //65
    R_bar_AD2=105*rho_A;  //75
    R_m_AD=1*(rho_A+rho_D)+rhoAD_safe+2;
    //R_m_AD=0;   For testing purposes
    R_bar_AD=R_m_AD+R_bar_AD1+rhoAD_safe;
    R_u_AD=R_m_AD+R_bar_AD2+rhoAD_safe;
    dR_AD_cube =pow((R_u_AD-R_bar_AD),3);
    A_A_D=2/dR_AD_cube;
    B_A_D=-3*(R_u_AD+R_bar_AD)/dR_AD_cube;
    C_A_D=6*R_u_AD*R_bar_AD/dR_AD_cube;
    D_A_D=pow(R_u_AD,2)*(R_u_AD-3*R_bar_AD)/dR_AD_cube;

    R_m_AD2=6*R_m_AD; //1*(rho_A+rho_D)+rhoAD_safe;
    //R_m_AD=0;  % For testing purposes
    R_bar_AD2=R_m_AD2+R_bar_AD1+rhoAD_safe;
    R_u_AD2=R_m_AD2+R_bar_AD2+rhoAD_safe;
    dR_AD_cube2=pow((R_u_AD2-R_bar_AD2),3);
    A_A_D2=2/dR_AD_cube2;
    B_A_D2=-3*(R_u_AD2+R_bar_AD2)/dR_AD_cube2;
    C_A_D2=6*R_u_AD2*R_bar_AD2/dR_AD_cube2;
    D_A_D2=pow(R_u_AD2,2)*(R_u_AD2-3*R_bar_AD2)/dR_AD_cube2;
}


std::vector<double> Rij0,Rjj0;
double rho_Fmax;
void fillR() {
    for (size_t i = 0; i < ND; i++)
    {
        Rij0.push_back(R_u_AD);
        Rjj0.push_back(20);
    }
    double RAD_max=0.2*R_bar_AD+.8*R_m_AD;
    rho_Fmax=RAD_max+rho_D;
}

double alpha=1;

double GO=500;
double R_safe=5*rho_D;
double E_OF1=2;  //repulsive
double E_OF2=16;  //blending of repulsive and attractive

void calVfield_formation() {
    for (size_t k = 0; k < NO; k++)
    {
        w.push_back(rVO[k][0][1] - rVO[k][0][0]);
        h.push_back(rVO[k][1][2] - rVO[k][1][1]);
        w_bar.push_back(w[k] + 2*(rho_Fmax+R_safe));
        h_bar.push_back(2*(rho_Fmax+R_safe));
        //////////////////////////////////////
        //                                  //
        //         fsolve problem           //
        //         368- 381                 //
        //////////////////////////////////////
    }
    
}

std::vector<std::vector<double>> R_m_AO, R_bar_AO, R_u_AO, R_v_AO;
std::vector<std::vector<double>> A_A_O, B_A_O, C_A_O, D_A_O;
std::vector<std::vector<double>> A_bar_A_O, B_bar_A_O, C_bar_A_O, D_bar_A_O;
void calVfield_attackers() {
    double R_bar_O1=7*rho_A;   //repulsive
    double R_bar_O2=17*rho_A;  //blending of repulsive and attractive
    double R_bar_O3=12;

    int default_value = 1;
 
    // resize the vector to `NA` elements of type std::vector<int>,
    // each having size `NO` and default value
    // NA*NO matrix
    R_m_AO.resize(NA, std::vector<double>(NO, default_value));
    R_bar_AO.resize(NA, std::vector<double>(NO, default_value));
    R_u_AO.resize(NA, std::vector<double>(NO, default_value));
    R_v_AO.resize(NA, std::vector<double>(NO, default_value));
    A_A_O.resize(NA, std::vector<double>(NO, default_value));
    B_A_O.resize(NA, std::vector<double>(NO, default_value));
    C_A_O.resize(NA, std::vector<double>(NO, default_value));
    D_A_O.resize(NA, std::vector<double>(NO, default_value));
    A_bar_A_O.resize(NA, std::vector<double>(NO, default_value));
    B_bar_A_O.resize(NA, std::vector<double>(NO, default_value));
    C_bar_A_O.resize(NA, std::vector<double>(NO, default_value));
    D_bar_A_O.resize(NA, std::vector<double>(NO, default_value));
    


    for (size_t j = 0; j < NA; j++)
    {
        for (size_t k = 0; k < NO; k++)
        {
            double rho_bar=2*10*rho_A;
            R_m_AO[j][k]=0;
            R_bar_AO[j][k]=2*30*rho_A;
            R_u_AO[j][k]=2*40*rho_A;
            R_v_AO[j][k]=2*40*rho_A;
            double dR_AO_cube=pow((R_u_AO[j][k]-R_bar_AO[j][k]),3);
            A_A_O[j][k]=2/dR_AO_cube;
            B_A_O[j][k]=-3*(R_u_AO[j][k]+R_bar_AO[j][k])/dR_AO_cube;
            C_A_O[j][k]=6*R_u_AO[j][k]*R_bar_AO[j][k]/dR_AO_cube;
            D_A_O[j][k]=pow(R_u_AO[j][k],2)*(R_u_AO[j][k]-3*R_bar_AO[j][k])/dR_AO_cube;
            
            double dR_bar_AO_cube=pow((R_v_AO[j][k]-R_u_AO[j][k]),3);
            A_bar_A_O[j][k]=2/dR_bar_AO_cube;
            B_bar_A_O[j][k]=-3*(R_v_AO[j][k]+R_u_AO[j][k])/dR_bar_AO_cube;
            C_bar_A_O[j][k]=6*R_v_AO[j][k]*R_u_AO[j][k]/dR_bar_AO_cube;
            D_bar_A_O[j][k]=pow(R_v_AO[j][k],2)*(R_v_AO[j][k]-3*R_u_AO[j][k])/dR_bar_AO_cube;
        }
        
    }
    
}


std::vector<std::vector<double>> A_D_O, B_D_O, C_D_O, D_D_O;
std::vector<std::vector<double>> R_m_DO, R_bar_DO, R_u_DO;
double R_m_DD, R_m_DDO, R_m2_DD, R_u_DD, R_bar_DD;
double A_D_D, B_D_D, C_D_D, D_D_D;

void calVfield_defenders() {
    double E_bar_O1=10*rho_D;  //repulsive
    double E_bar_O2=20*rho_D;  //blending of repulsive and attractive
    double E_bar_O3=5;

    int default_value = 1;
 
    // resize the vector to `NA` elements of type std::vector<int>,
    // each having size `NO` and default value
    // NA*NO matrix
    R_m_DO.resize(ND, std::vector<double>(NO, default_value));
    R_bar_DO.resize(ND, std::vector<double>(NO, default_value));
    R_u_DO.resize(ND, std::vector<double>(NO, default_value));
    A_D_O.resize(ND, std::vector<double>(NO, default_value));
    B_D_O.resize(ND, std::vector<double>(NO, default_value));
    C_D_O.resize(ND, std::vector<double>(NO, default_value));
    D_D_O.resize(ND, std::vector<double>(NO, default_value));
    
    for (size_t j = 0; j < ND; j++)
    {
        for (size_t k = 0; k < NO; k++)
        {
            double rho_bar=2*10*rho_D;
            R_m_DO[j][k]=0;
            R_bar_DO[j][k]=2*30*rho_D;
            R_u_DO[j][k]=2*40*rho_D;
            double dR_DO_cube=pow((R_u_DO[j][k]-R_bar_DO[j][k]),3);
            A_D_O[j][k]=2/dR_DO_cube;
            B_D_O[j][k]=-3*(R_u_DO[j][k]+R_bar_DO[j][k])/dR_DO_cube;
            C_D_O[j][k]=6*R_u_DO[j][k]*R_bar_DO[j][k]/dR_DO_cube;
            D_D_O[j][k]=pow(R_u_DO[j][k],2)*(R_u_DO[j][k]-3*R_bar_DO[j][k])/dR_DO_cube;
            
        }
        
    }

    //for other defenders
    double R_bar_DD1=6*rho_D;
    double R_bar_DD2=6*rho_D;
    R_m_DD = 1*(rho_D+rho_D)+rhoD_safe;
    R_m_DDO=3*R_m_DD;
    R_m2_DD=100;
    R_bar_DD=R_m_DD+R_bar_DD1+rhoD_safe;
    R_u_DD=R_bar_DD+R_bar_DD2+rhoD_safe;
    double dR_DD_cube = pow((R_u_DD-R_bar_DD),3);
    A_D_D=2/dR_DD_cube;
    B_D_D=-3*(R_u_DD+R_bar_DD)/dR_DD_cube;
    C_D_D=6*R_u_DD*R_bar_DD/dR_DD_cube;
    D_D_D=pow(R_u_DD,2)*(R_u_DD-3*R_bar_DD)/dR_DD_cube;

    
}

//Initialize the defenders
double rD0[2] = {-150,-300};
std::vector<std::vector<double>> rD;
std::vector<std::vector<double>> XD0;
std::vector<std::vector<double>> rSD_goal;
void InitializeDefenders() {
    rD.resize(2, std::vector<double>(ND, 1));
    XD0.resize(4, std::vector<double>(ND, 0));
    rSD_goal.resize(2, std::vector<double>(ND, 0));
    for (size_t j = 0; j < ND; j++)
    {
        double RD0=20+(j)*8;//3*rDmin;
        double thetaD=(j)*(2*M_PI-M_PI/3)/ND;
        rD[0][j]=rD0[0] + RD0*cos(thetaD);
        rD[1][j]=rD0[1] + RD0*sin(thetaD);
        XD0[0][j]=rD[0][j];
        XD0[1][j]=rD[1][j];
        rSD_goal[0][j]=rS[0]+RD0*cos(thetaD);
        rSD_goal[1][j]=rS[1]+RD0*sin(thetaD);
    }
    
}

void rD_value(int s) {
    if (s == 1)
    {
        rD = {{-204,-60,-500,300,-150},
             {-250,-365,500,-100,500}};
    } else if (s == 2)
    {
        rD = {{-740,-180,-800,270,-150},
            {-200,-258,100,450,500}};
    } else if (s == 3)
    {
        rD = {{-740,-180,-800,370,-150,-500},
             {-200,-258,100,-390,500,400}};
    }
    XD0.resize(4, std::vector<double>(rD[0].size(), 0));
    for (size_t i = 0; i < rD[0].size(); i++)
    {
        XD0[0][i]=rD[0][i];
        XD0[1][i]=rD[1][i];
    }
    
}

int main()
{
 
    // Function declared in header
    // file to find the sum
    calControlLimits();
    VectorFields_A();
    initialAttackersVel();
    cal();
    InitializeAttackers(sceneraio);
    defenders();
    FormationOri();
    calVfieldA();
    fillR();
    calVfield_attackers();
    calVfield_defenders();
    InitializeDefenders();
    rD_value(sceneraio);
    for (size_t i = 0; i < XD0.size(); i++)
    {
        for (size_t j = 0; j < XD0[0].size(); j++)
        {
            cout << XD0[i][j] << "   ";
        }
        cout << "\n";
    }
    
    
}