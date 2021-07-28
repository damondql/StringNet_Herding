#pragma once

#include "iostream"
#include "AllParametersExperiment.hpp" 
#include <math.h>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <armadillo>
using namespace std;
using namespace arma;
 
// Driver Code
double dt = 0.25;
double Rs = 1000;
double rho_P = 1;

vec rP(2);
vec rS(2);

int wW=20; // World width
int lW=30; //World length

double rho_A = 0.4;
double rho_c_A = 10;

double C_d = 0.4;
double vma = 0.22;
double uma;
double rhoA_safe;

void calControlLimits() {
    rP = {14,-1};
    rS = {22,-36};
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
    R_m_AA=1*(rho_A+rho_A)+rhoA_safe;
    R_bar_AA=R_m_AA+R_bar_AA1+rhoA_safe;
    R_u_AA=R_m_AA+R_bar_AA2+rhoA_safe;
    dR_AA_cube=pow((R_u_AA-R_bar_AA),3);
    A_A_A=2/dR_AA_cube;
    B_A_A=-3*(R_u_AA+R_bar_AA)/dR_AA_cube;
    C_A_A=6*R_u_AA*R_bar_AA/dR_AA_cube;
    D_A_A=pow(R_u_AA,2)*(R_u_AA-3*R_bar_AA)/dR_AA_cube;
}

vec v_maxA(NA), u_maxA(NA);

void initialAttackersVel(){
    for (size_t i = 0; i < NA; i++)
    {
       v_maxA(i) = vma;
       u_maxA(i) = uma;
    }
}

double NA_sep = 0;

int axS[2] = {1, 1};
// rDmin = rho_P+10
int rDmin = 7;


double RA0;
double potential;
vec Rii00(NA), Rik00(NA);

void cal(){
    RA0 = R_m_AA + 0.1;
    potential = RA0 * sqrt(2*(1-cos(2*M_PI/NA))); //for formation potential
    for (size_t i = 0; i < NA; i++)
    {
        Rii00(i) = potential;
        Rik00(i) = 1 * Rii00(i);
    }
}

int dRA0 = 1;

int sceneraio = 1;
double RA01 = RA0;
double thetaA0;
vec rA0(2);
mat rA(2,NA, fill::zeros);
mat vA(2,NA, fill::zeros);
mat XA0(4,NA);
mat rA_follow(2,NA);
double RA;
double rho_S;
double rho_Acon;

void InitializeAttackers() {

    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 0; j < NA; j++)
        {
            rA(i,j) = 0;
            vA(i,j) = 0;
        }
        
    }
    rA0(0) = 6;
    rA0(1) = -lW+4;
    rA(0,0) = 6.1548;
    rA(1,0) = -33.0553;
    vA(0,0) = 0;
    vA(1,0) = 0;

    for (size_t i = 0; i < NA; i++)
    {
        XA0(0,i) = rA(0,i);
        XA0(1,i) = rA(1,i);
        XA0(2,i) = vA(0,i);
        XA0(3,i) = vA(1,i);
    }
    
    for (size_t i = 0; i < NA; i++)
    {
        rA_follow(0,i) = rA(0,i) - rA(0,0);
        rA_follow(1,i) = rA(1,i) - rA(1,0);
    }
    
    //XA0(:,NA)=[-700,0,0,0]';
    double rAcm0[2] = {0, 0};
    double vAcm0[2] = {0, 0};
    
    for (size_t i = 0; i < NA; i++)
    {
        rAcm0[0] += XA0(0,i);
        rAcm0[1] += XA0(1,i);
        vAcm0[0] += XA0(2,i);
        vAcm0[1] += XA0(3,i);
    }
    rAcm0[0] = rAcm0[0] / NA;
    rAcm0[1] = rAcm0[1] / NA;
    vAcm0[0] = vAcm0[0] / NA;
    vAcm0[1] = vAcm0[1] / NA;

    double thetaAcm0 = atan2(rAcm0[1]-rP(1),rAcm0[0]-rP(0));
    rho_Acon = 1.3*RA0;  //radius of the connectivity region
    
    rho_S= 13;   //radius of the safe area

}
// Defenders

double rho_D = 0.4;
double rho_sn_max = 65;
double ND = 3;
double N;
double rho_sn;
vec Rjk00(ND);
arma::vec v_maxD,v_maxDC,u_maxD,u_maxD1,u_maxD2,u_maxDr1,u_maxDr2;
double rho_safe,dthetai,alphaD_v;
double rhoAD_safe,rhoD_safe;
double umd1,umd2,umdf_s1,umdf_s2,umd,umdf_h1,vmd,vmdc,vmdf_s,umd_e1,umd_e2;
obstacle obs;
void defenders() {
    v_maxD.resize(ND);
    v_maxD.resize(ND);
    v_maxDC.resize(ND);
    u_maxD.resize(ND);
    u_maxD1.resize(ND);
    u_maxD2.resize(ND);
    u_maxDr1.resize(ND);
    u_maxDr2.resize(ND);
    N = NA + ND;
    alphaD_v=0.5;
    umdf_h1=uma;
    double kDv=.1;
    double alphaD_r=alphaD_v/(2-alphaD_v);
    umdf_s2=0.6*uma;
    umdf_s1=umdf_s2;
    double eta_max=vma+2;
    umd2=0.03;
    umd1=1.5*umd2;
    umd=umd1+umd2+umdf_s1+umdf_s2;
    umd_e1=0.75*umd;
    umd_e2=0.25*umd;
    if (C_d != 0) {
        vmd=sqrt(umd/C_d);
        rhoD_safe=log((umd + C_d*pow((2*vmd),2))/umd)/(2*C_d);   //additional aafety distance required between the defenders due to bounded acceleration
        rhoAD_safe=log((uma +C_d*pow((vma+vmd),2))/uma)/(2*C_d);  //additional safty for the attackers to stay safe from the defenders due to bounded acceleration
        rho_safe = rhoD_safe;
        vmdc=sqrt((sqrt(pow(umd,2)/(1/pow(rho_safe,2)+pow(C_d,2)))));
    } else {
        u_maxD(0)=(0.4*vmd);    // Ask here: what is the purpose of doing this since it will be assign value below
        rhoD_safe=pow((vmd+vmd),2)/(2*umd);
        rhoAD_safe=pow((vma+vmd),2)/(2*uma);
        rho_safe = rhoD_safe;
        vmdc=sqrt(umd*rhoD_safe);
    }
    vmdf_s=sqrt((umdf_s1+umdf_s2)/C_d);
    for (size_t i = 0; i < ND; i++)
    {
        v_maxD(i) = vmd;
        v_maxDC(i) = vmdc;
        u_maxD(i) = umd;
        u_maxD1(i) = umd1;
        u_maxD2(i) = umd2;
        u_maxDr1(i) = umdf_s1;
        u_maxDr2(i) = umdf_s2;
    }
    obs.rho_safe = rho_safe;
    dthetai=acos(1-pow((2*rho_D),2)/(2*pow(rho_safe,2)));   //angular shift for two agents colliding on a ciruclar arc segment
    rho_sn=rho_Acon+rhoD_safe;   //radius of the stringNet
    for (size_t i = 0; i < ND; i++)
    {
        Rjk00(i)=(Rik00(0));
    }
}

int largeP = 100000;

////////////////////////////////////
//  RD_con fsolve problem!        //
////////////////////////////////////

double kr1=0.5;
double kr2=0.5;
double kv1=0.1;
double kv2=0.5;

//control gains for attackers
double kAFr=1;
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
double kAPr=0.009;
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
double kDDesr=2;
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
cube rVO(2,4,0);

double NO = rVO.n_slices;
double NVOk;
vec rOc(2,fill::zeros);
vec dVO;

vec R_bar_AcOc, R_u_AcOc, A_Ac_Oc, B_Ac_Oc, C_Ac_Oc, D_Ac_Oc;
void FormationOri() {
    // rVO.slice(0) = {{60,160,160,60},
    //                 {310,310,390,390}};
    // rVO.slice(1) = {{-300,-190,-190,-300},
    //                 {350,350,560,560}};
    // cout << "set rV0 vaule success" << endl;
    alphaDFr=alphaDFv/(2-alphaDFv);
    kDFr2=(umd2+C_d*pow(vmd,2))/pow((vmd+vma),(alphaDFv));
    for (size_t k = 0; k < NO; k++)
    {
        NVOk = rVO.n_cols;
        for (size_t i = 0; i < rVO.n_cols;i++) {
            rOc(0) += rVO(0,i,k);
            rOc(1) += rVO(1,i,k);
        }
        rOc(0) = rOc(0) / NVOk;
        rOc(1) = rOc(1) / NVOk;
        for (size_t j = 0; j < NVOk; j++) {
            double a = rOc(0) - rVO(0,j,k);
            double b = rOc(1) - rVO(1,j,k);
            std::complex<double> comp = {a,b};
            dVO(j)=sqrt(std::norm(comp));
        }
        
        double max_dVO = dVO.max();
        R_bar_AcOc(k)= (rho_Acon + max_dVO + 5);
        R_u_AcOc(k)=(rho_Acon + max_dVO +155);
        double dR_AcOc_cube=pow((R_u_AcOc[k]-R_bar_AcOc[k]),3);
        A_Ac_Oc(k) = (2/dR_AcOc_cube);
        B_Ac_Oc(k) = (-3*(R_u_AcOc[k]+R_bar_AcOc[k])/dR_AcOc_cube);
        C_Ac_Oc(k) = (6*R_u_AcOc[k]*R_bar_AcOc[k]/dR_AcOc_cube);
        D_Ac_Oc(k) = (pow(R_u_AcOc[k],2)*(R_u_AcOc[k]-3*R_bar_AcOc[k])/dR_AcOc_cube);

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
    R_bar_AD1=0*rho_A;  //65
    R_bar_AD2=0.5*rho_A;  //75
    R_m_AD= 2*(rho_A+rho_D)+rhoA_safe;
    //R_m_AD=0;   For testing purposes
    R_bar_AD=R_m_AD+R_bar_AD1+rhoA_safe;
    R_u_AD=R_m_AD+R_bar_AD2+rhoA_safe;
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


vec Rij0(ND),Rjj0(ND);
double rho_Fmax;
double RAD_max;
void fillR() {
    for (size_t i = 0; i < ND; i++)
    {
        Rij0(i) = (R_u_AD);
        Rjj0(i) = 20;
    }
    RAD_max=0.2*R_bar_AD+.8*R_m_AD;
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
        w.push_back(rVO(0,1,k) - rVO(0,0,k));
        h.push_back(rVO(1,2,k) - rVO(1,1,k));
        w_bar.push_back(w[k] + 2*(rho_Fmax+R_safe));
        h_bar.push_back(2*(rho_Fmax+R_safe));
        //////////////////////////////////////
        //                                  //
        //         fsolve problem           //
        //         368- 381                 //
        //////////////////////////////////////
    }
    
}

// std::vector<std::vector<double>> R_m_AO, R_bar_AO, R_u_AO, R_v_AO;
// std::vector<std::vector<double>> A_A_O, B_A_O, C_A_O, D_A_O;
// std::vector<std::vector<double>> A_bar_A_O, B_bar_A_O, C_bar_A_O, D_bar_A_O;

mat R_m_AO, R_bar_AO, R_u_AO, R_v_AO;
mat A_A_O, B_A_O, C_A_O, D_A_O;
mat A_bar_A_O, B_bar_A_O, C_bar_A_O, D_bar_A_O;
void calVfield_attackers() {
    double R_bar_O1=7*rho_A;   //repulsive
    double R_bar_O2=17*rho_A;  //blending of repulsive and attractive
    double R_bar_O3=12;

    int default_value = 1;
 
    // resize the vector to `NA` elements of type std::vector<int>,
    // each having size `NO` and default value
    // NA*NO matrix
    R_m_AO.resize(NA,NO);
    R_bar_AO.resize(NA,NO);
    R_u_AO.resize(NA,NO);
    R_v_AO.resize(NA,NO);
    A_A_O.resize(NA,NO);
    B_A_O.resize(NA,NO);
    C_A_O.resize(NA,NO);
    D_A_O.resize(NA,NO);
    A_bar_A_O.resize(NA,NO);
    B_bar_A_O.resize(NA,NO);
    C_bar_A_O.resize(NA,NO);
    D_bar_A_O.resize(NA,NO);
    


    for (size_t j = 0; j < NA; j++)
    {
        for (size_t k = 0; k < NO; k++)
        {
            double rho_bar=2*10*rho_A;
            R_m_AO(j,k)=0;
            R_bar_AO(j,k)=2*30*rho_A;
            R_u_AO(j,k)=2*40*rho_A;
            R_v_AO(j,k)=2*40*rho_A;
            double dR_AO_cube=pow((R_u_AO(j,k)-R_bar_AO(j,k)),3);
            A_A_O(j,k)=2/dR_AO_cube;
            B_A_O(j,k)=-3*(R_u_AO(j,k)+R_bar_AO(j,k))/dR_AO_cube;
            C_A_O(j,k)=6*R_u_AO(j,k)*R_bar_AO(j,k)/dR_AO_cube;
            D_A_O(j,k)=pow(R_u_AO(j,k),2)*(R_u_AO(j,k)-3*R_bar_AO(j,k))/dR_AO_cube;
            
            double dR_bar_AO_cube=pow((R_v_AO(j,k)-R_u_AO(j,k)),3);
            A_bar_A_O(j,k)=2/dR_bar_AO_cube;
            B_bar_A_O(j,k)=-3*(R_v_AO(j,k)+R_u_AO(j,k))/dR_bar_AO_cube;
            C_bar_A_O(j,k)=6*R_v_AO(j,k)*R_u_AO(j,k)/dR_bar_AO_cube;
            D_bar_A_O(j,k)=pow(R_v_AO(j,k),2)*(R_v_AO(j,k)-3*R_u_AO(j,k))/dR_bar_AO_cube;
        }
        
    }
    
}


mat A_D_O, B_D_O, C_D_O, D_D_O;
mat R_m_DO, R_bar_DO, R_u_DO;
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
    R_m_DO.resize(ND,NO);
    R_bar_DO.resize(ND,NO);
    R_u_DO.resize(ND,NO);
    A_D_O.resize(ND,NO);
    B_D_O.resize(ND,NO);
    C_D_O.resize(ND,NO);
    D_D_O.resize(ND,NO);
    
    for (size_t j = 0; j < ND; j++)
    {
        for (size_t k = 0; k < NO; k++)
        {
            double rho_bar=2*10*rho_D;
            R_m_DO(j,k)=0;
            R_bar_DO(j,k)=2*30*rho_D;
            R_u_DO(j,k)=2*40*rho_D;
            double dR_DO_cube=pow((R_u_DO(j,k)-R_bar_DO(j,k)),3);
            A_D_O(j,k)=2/dR_DO_cube;
            B_D_O(j,k)=-3*(R_u_DO(j,k)+R_bar_DO(j,k))/dR_DO_cube;
            C_D_O(j,k)=6*R_u_DO(j,k)*R_bar_DO(j,k)/dR_DO_cube;
            D_D_O(j,k)=pow(R_u_DO(j,k),2)*(R_u_DO(j,k)-3*R_bar_DO(j,k))/dR_DO_cube;
            
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
mat rD;
arma::mat XD0;
std::vector<std::vector<double>> rSD_goal;
mat XD;
void rD_value() {

    rD = {{11.9907, 9.3724,	16.1838,},
          {-10.7580, -7.1278, -5.8411}};
   
    XD0.resize(4, rD.n_cols);
    XD0.zeros();
    XD.set_size(4, rD.n_cols);
    XD.zeros();
    for (size_t i = 0; i < rD.n_cols; i++)
    {
        XD0(0,i)=rD(0,i);
        XD0(1,i)=rD(1,i);
        XD(0,i)=rD(0,i);
        XD(1,i)=rD(1,i);
    }
    
}

void calAllparametersExperiment() {
 
    // Function declared in header
    // file to find the sum
    calControlLimits();
    VectorFields_A();
    initialAttackersVel();
    cal();
    InitializeAttackers();
    defenders();
    FormationOri();
    calVfieldA();
    fillR();
    calVfield_attackers();
    calVfield_defenders();
    rD_value();
    
    
}