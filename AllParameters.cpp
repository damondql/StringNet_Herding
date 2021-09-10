#pragma once
#include "AllParameters.hpp" 

using namespace std;
using namespace arma;
 
// Driver Code
double dt = 0.01;
double Rs = 1000;
double rho_P = 45;

vec rP(2);
mat rS(2,4);

double wW=20.0; // World width
double lW=30.0; //World length

double rho_A = 0.5;
double rho_c_A = 150;

double C_d = 0.2;
double vma = 6;
double uma;
double rhoA_safe;

void calControlLimits() {
    rP = {0,200};
    rS = {{1000, -1000, -1000,  1000},
          {1000,  1000, -1000, -1000}};
    rS = 1.5*rS;
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
    R_bar_AA1=4*rho_A;
    R_bar_AA2=8*rho_A;
    R_m_AA=1*(rho_A+rho_A)+rhoA_safe+1;
    R_bar_AA=R_m_AA+R_bar_AA1;
    R_u_AA=R_m_AA+R_bar_AA2;
    dR_AA_cube=pow((R_u_AA-R_bar_AA),3);
    A_A_A=2/dR_AA_cube;
    B_A_A=-3*(R_u_AA+R_bar_AA)/dR_AA_cube;
    C_A_A=6*R_u_AA*R_bar_AA/dR_AA_cube;
    D_A_A=pow(R_u_AA,2)*(R_u_AA-3*R_bar_AA)/dR_AA_cube;
}

vec v_maxA, u_maxA;

void initialAttackersVel(int NA){
    v_maxA.resize(NA);
    u_maxA.resize(NA);
    for (size_t i = 0; i < NA; i++)
    {
       v_maxA(i) = vma;
       u_maxA(i) = uma;
    }
}

double NA_sep = 0;

int axS[2] = {1, 1};
// // rDmin = rho_P+10
double rDmin = 7.0;


double RA0;
double potential;
vec Rii00, Rik00;

void cal(int NA){
    Rii00.resize(NA);
    Rik00.resize(NA);
    RA0 = R_m_AA + 8;
    potential = RA0 * sqrt(2*(1-cos(2*M_PI/NA))); //for formation potential
    for (size_t i = 0; i < NA; i++)
    {
        Rii00(i) = potential;
        Rik00(i) = 5 * Rii00(i);
    }
}

int dRA0 = 10;

int sceneraio = 1;
double RA01;
double thetaA0;
vec rA0(2);
mat rA;
mat vA;
mat XA0;
mat rA_follow;
double RA;
double rho_S;
double rho_Acon;
arma::mat clusteridA;
arma::field<vec> indAinClusterA(1);
arma::mat indALeader;
// arma::mat rA_follow;
arma::mat leaderIDA;
vec NAinClusterA;

void InitializeAttackers(int NA) {
    // cout << "start initialized attackers" << endl;
    RA01 = RA0;
    rA.resize(2,NA);
    vA.resize(2,NA);
    XA0.resize(4,NA);
    rA_follow.resize(2,NA);

    rA0 = {1345, -210};
    thetaA0 = atan2(rA0(1) - rP(1), rA0(0)-rP(0)) + 5 * M_PI/6;
    // cout << "enter for loop" << endl;
    for (int i = 0; i < NA; i++)
    {
        double thetaA = thetaA0 + 2*M_PI*(i+1)/NA;
        RA = 2 * RA0;
        // cout << "RA: " << RA << endl;
        // cout << "thetaA" << thetaA << endl;
        rA.col(i) = {RA * thetaA, RA * sin(thetaA)};
        rA.col(i) += rA0;
        vA.col(i) = -(rA.col(i) - rP)*0.01 * v_maxA(i) / arma::norm(rA.col(i) - rP);
    }
    // rA.print("rA: ");
    // vA.print("vA: ");
    XA0 = join_cols(rA,vA);
    int NClusterA = 1;

    
    clusteridA = arma::ones(NA,1);
    indALeader = arma::zeros(NA,1);
    rA_follow = arma::zeros(2,NA);
    indAinClusterA(0) = {9,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18};
    indAinClusterA(0) = indAinClusterA(0) - 1;
    NAinClusterA.resize(NClusterA);
    arma::mat ind;
    // cout << "enter for loop" << endl;
    for (int i = 0; i < NClusterA; i++)
    {
        NAinClusterA(i) = indAinClusterA(i).n_elem;
        indALeader(indAinClusterA(i)(0)) = 1; //ind of leader is 9, so ind in matrix is 8
        leaderIDA.resize(NAinClusterA(i),1);
        for (int j = 0; j < NAinClusterA(i); j++)
        {
            leaderIDA(j) = indAinClusterA(i)(0);
            // cout << "LeadID: " << indAinClusterA(i)(0) <<endl;
        }
        ind = indAinClusterA(i).subvec(1,NAinClusterA(i)-1);
        // ind.print("ind: ");
        for (int j = 0; j < ind.n_elem; j++)
        {
            rA_follow.col(ind(j)) = rA.col(ind(j)) - rA.col(indAinClusterA(i)(0));
        }
        
        
    }
    // leaderIDA.print();
    // rA_follow.print();
    // indAinClusterA.print();
    // indALeader.print("indAleader:");

    vec rAcm0 = arma::sum(XA0.submat(0,0,1,NA-1),1)/NA;
    vec vAcm0 = arma::sum(XA0.submat(2,0,3,NA-1),1)/NA;

    // rAcm0.print("rAcm0: ");
    // vAcm0.print("vAcm0:");
    
    double thetaAcm0 = atan2(rAcm0(1)-rP(1), rAcm0(0)-rP(0));
    rho_Acon = 2 * RA0 + rhoA_safe + dRA0 + 10;
    double RAcon = RA0;
    rho_S = 4*rho_Acon;
    // cout << "rho_S: "<<rho_S<< endl;
}









// Defenders

double rho_D = 0.5;
double rho_sn_max = 60;

double N;
double rho_sn;
vec Rjk00;
arma::vec v_maxD,v_maxDC,u_maxD,u_maxD1,u_maxD2,u_maxDr1,u_maxDr2;
double rho_safe,dthetai,alphaD_v;
double rhoAD_safe,rhoD_safe;
// double umd1,umd2,umdf_s1,umdf_s2,umd,umdf_h1,vmd,vmdc,vmdf_s,umd_e1,umd_e2;
// obstacle obs;
void defenders(int NA, int ND) {
    double vmd, vmdc;
    Rjk00.resize(ND);
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
    // umdf_h1=uma;
    double kDv=.1;
    double alphaD_r=alphaD_v/(2-alphaD_v);
    double umdr2=uma+kDv*vma;
    double umdr1=0.5*umdr2;
    double eta_max=vma+2;
    double umd2 = max(umdr1+umdr2,C_d*pow(eta_max,2)+pow(eta_max,(alphaD_v)));
    double umd1 = umd2;
    double umd = umd1+umd2;
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
    // vmdf_s=sqrt((umdf_s1+umdf_s2)/C_d);
    for (size_t i = 0; i < ND; i++)
    {
        v_maxD(i) = vmd;
        v_maxDC(i) = vmdc;
        u_maxD(i) = umd;
        u_maxD1(i) = umd1;
        u_maxD2(i) = umd2;
        u_maxDr1(i) = umdr1;
        u_maxDr2(i) = umdr2;
    }
    // obs.rho_safe = rho_safe;
    dthetai=acos(1-pow((2*rho_D),2)/(2*pow(rho_safe,2)));   //angular shift for two agents colliding on a ciruclar arc segment
    rho_sn=rho_Acon+rhoD_safe +10;   //radius of the stringNet
    for (size_t i = 0; i < ND; i++)
    {
        Rjk00(i)=(Rik00(0));
    }
    // cout << "rhoD_safe: " <<rhoD_safe << endl;
    // cout << "rhoAD_safe: " <<rhoAD_safe << endl;
    // cout << "rho_safe: " <<rho_safe << endl;
    // cout << "vmdc: " <<vmdc << endl;

    // v_maxD.print("v_maxD");
    // v_maxDC.print("v_maxDC");
    // u_maxD.print("u_maxD");
    // u_maxD1.print("u_maxD1");
    // u_maxD2.print("u_maxD2");
    // u_maxDr1.print("u_maxDr1");
    // u_maxDr2.print("u_maxDr2");
    // cout << "dthetai: "<<dthetai << endl;
    // cout << "rho_sn: " << rho_sn << endl;
    // Rjk00.print("Rjk00");
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
double kAPr=0.006;
double kAPv=0.1;
double kAPr2 = 0.1;
double kAPv2 = 1;

//finite time gains for the defenders
double kDOr=1;
double kDDr=.5;
double kDOv1=2;
double kDOv2=0.1;
double kDFr=5.5;
double kDFr2 = 0.3;
double kDFv=.45;
double alphaDFv=.9;
double alphaDFr;
// double kDFr2;
double alphaDOv=1;
double kDDv=0.3;
double alphaDDv=0.5;
double kDRr=.0001;
double kDRv=0.1;
double kDFphi=0.003;
double kDFphid=0.007;
double kDFphir=0.3;
double kDFphiv=0.1;
double kDDesr=2;
double kDDesv=1;

// //////////////////////////////////////////////////////////////////////
// ////      %bound on the convergence error during tracking         ////
// ////      A_tilde=[zeros(2),eye(2);-kDFr*eye(2),-kDFv*eye(2)];    ////
// ////      Qd=1*eye(4);                                            ////
// ////      Pd=lyap(A_tilde',Qd);                                   ////
// ////      c1=min(eig(Pd));                                        ////
// ////      c2=max(eig(Pd));                                        ////
// ////      c3=min(eig(Qd));                                        ////
// ////      c4=2*max(eig(Pd));                                      ////
// ////      p=c4/c3*sqrt(c2/c1);                                    ////
// ////      bd=15;                                                  ////
// ////   Problem: lyap conuld not use coder to covert,              ////
// ////            and did not see p use anywhere else, it is not    ////
// ////            define as global                                  ////
// //////////////////////////////////////////////////////////////////////

double bd=15;   
///////////////////////////////////
//Convex polygonal Obstalces///////
///////////////////////////////////
cube rVO(2,4,0);

double NO = rVO.n_slices;
// double NVOk;
// vec rOc(2,fill::zeros);
// vec dVO;

// vec R_bar_AcOc, R_u_AcOc, A_Ac_Oc, B_Ac_Oc, C_Ac_Oc, D_Ac_Oc;
// void FormationOri() {
//     // rVO.slice(0) = {{60,160,160,60},
//     //                 {310,310,390,390}};
//     // rVO.slice(1) = {{-300,-190,-190,-300},
//     //                 {350,350,560,560}};
//     // cout << "set rV0 vaule success" << endl;
//     alphaDFr=alphaDFv/(2-alphaDFv);
//     kDFr2=(umd2+C_d*pow(vmd,2))/pow((vmd+vma),(alphaDFv));
//     for (size_t k = 0; k < NO; k++)
//     {
//         NVOk = rVO.n_cols;
//         for (size_t i = 0; i < rVO.n_cols;i++) {
//             rOc(0) += rVO(0,i,k);
//             rOc(1) += rVO(1,i,k);
//         }
//         rOc(0) = rOc(0) / NVOk;
//         rOc(1) = rOc(1) / NVOk;
//         for (size_t j = 0; j < NVOk; j++) {
//             double a = rOc(0) - rVO(0,j,k);
//             double b = rOc(1) - rVO(1,j,k);
//             std::complex<double> comp = {a,b};
//             dVO(j)=sqrt(std::norm(comp));
//         }
        
//         double max_dVO = dVO.max();
//         R_bar_AcOc(k)= (rho_Acon + max_dVO + 5);
//         R_u_AcOc(k)=(rho_Acon + max_dVO +155);
//         double dR_AcOc_cube=pow((R_u_AcOc[k]-R_bar_AcOc[k]),3);
//         A_Ac_Oc(k) = (2/dR_AcOc_cube);
//         B_Ac_Oc(k) = (-3*(R_u_AcOc[k]+R_bar_AcOc[k])/dR_AcOc_cube);
//         C_Ac_Oc(k) = (6*R_u_AcOc[k]*R_bar_AcOc[k]/dR_AcOc_cube);
//         D_Ac_Oc(k) = (pow(R_u_AcOc[k],2)*(R_u_AcOc[k]-3*R_bar_AcOc[k])/dR_AcOc_cube);

//         rOc[0] = 0;
//         rOc[1] = 0;
//         if (k != NO - 1) {
//         dVO.clear();
//         }
//     }


// }
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
    R_bar_AD1= 95*rho_A;  //65
    R_bar_AD2= 105*rho_A;  //75
    R_m_AD= 1*(rho_A+rho_D)+rhoAD_safe+2;
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
    // cout << "rho_A: " <<rho_A << endl;
    // cout << "R_bar_AD1: " <<R_bar_AD1 << endl;
    // cout << "R_bar_AD2: " <<R_bar_AD2 << endl;
    // cout << "R_u_AD: " <<R_u_AD << endl;
    // cout << "R_bar_AD: " <<R_bar_AD << endl;
    // cout << "R_m_AD: " <<R_m_AD << endl;
    // cout << "A_A_D: " <<A_A_D << endl;
    // cout << "B_A_D: " <<B_A_D << endl;
    // cout << "C_A_D: " <<C_A_D << endl;
    // cout << "D_A_D: " <<D_A_D << endl;
    // cout << "A_A_D2: " <<A_A_D2 << endl;
    // cout << "B_A_D2: " <<B_A_D2 << endl;
    // cout << "C_A_D2: " <<C_A_D2 << endl;
    // cout << "D_A_D2: " <<D_A_D2 << endl;

}


vec Rij0,Rjj0;
double rho_Fmax;
double RAD_max;
void fillR(int ND) {
    Rij0.resize(ND);
    Rjj0.resize(ND);
    for (size_t i = 0; i < ND; i++)
    {
        Rij0(i) = (R_u_AD);
        Rjj0(i) = 20;
    }
    RAD_max=0.2*R_bar_AD+.8*R_m_AD;
    rho_Fmax=RAD_max+rho_D;
    // Rij0.print("Rij0");
    // Rjj0.print("Rjj0");
    // cout << "RAD_max: " << RAD_max << endl;
    // cout << "rho_Fmax: " << rho_Fmax << endl;  
}

double alpha=1;

double GO=500;
double R_safe;
double E_OF1=2;  //repulsive
double E_OF2=16;  //blending of repulsive and attractive

// void calVfield_formation() {
//     R_safe = 5*rho_D;
//     for (size_t k = 0; k < NO; k++)
//     {
//         w.push_back(rVO(0,1,k) - rVO(0,0,k));
//         h.push_back(rVO(1,2,k) - rVO(1,1,k));
//         w_bar.push_back(w[k] + 2*(rho_Fmax+R_safe));
//         h_bar.push_back(2*(rho_Fmax+R_safe));
//         //////////////////////////////////////
//         //                                  //
//         //         fsolve problem           //
//         //         368- 381                 //
//         //////////////////////////////////////
//     }
    
// }

// // std::vector<std::vector<double>> R_m_AO, R_bar_AO, R_u_AO, R_v_AO;
// // std::vector<std::vector<double>> A_A_O, B_A_O, C_A_O, D_A_O;
// // std::vector<std::vector<double>> A_bar_A_O, B_bar_A_O, C_bar_A_O, D_bar_A_O;

// mat R_m_AO, R_bar_AO, R_u_AO, R_v_AO;
// mat A_A_O, B_A_O, C_A_O, D_A_O;
// mat A_bar_A_O, B_bar_A_O, C_bar_A_O, D_bar_A_O;
// void calVfield_attackers(int NA) {
//     double R_bar_O1=7*rho_A;   //repulsive
//     double R_bar_O2=17*rho_A;  //blending of repulsive and attractive
//     double R_bar_O3=12;

//     int default_value = 1;
 
//     // resize the vector to `NA` elements of type std::vector<int>,
//     // each having size `NO` and default value
//     // NA*NO matrix
//     R_m_AO.resize(NA,NO);
//     R_bar_AO.resize(NA,NO);
//     R_u_AO.resize(NA,NO);
//     R_v_AO.resize(NA,NO);
//     A_A_O.resize(NA,NO);
//     B_A_O.resize(NA,NO);
//     C_A_O.resize(NA,NO);
//     D_A_O.resize(NA,NO);
//     A_bar_A_O.resize(NA,NO);
//     B_bar_A_O.resize(NA,NO);
//     C_bar_A_O.resize(NA,NO);
//     D_bar_A_O.resize(NA,NO);
    


//     for (size_t j = 0; j < NA; j++)
//     {
//         for (size_t k = 0; k < NO; k++)
//         {
//             double rho_bar=2*10*rho_A;
//             R_m_AO(j,k)=0;
//             R_bar_AO(j,k)=2*30*rho_A;
//             R_u_AO(j,k)=2*40*rho_A;
//             R_v_AO(j,k)=2*40*rho_A;
//             double dR_AO_cube=pow((R_u_AO(j,k)-R_bar_AO(j,k)),3);
//             A_A_O(j,k)=2/dR_AO_cube;
//             B_A_O(j,k)=-3*(R_u_AO(j,k)+R_bar_AO(j,k))/dR_AO_cube;
//             C_A_O(j,k)=6*R_u_AO(j,k)*R_bar_AO(j,k)/dR_AO_cube;
//             D_A_O(j,k)=pow(R_u_AO(j,k),2)*(R_u_AO(j,k)-3*R_bar_AO(j,k))/dR_AO_cube;
            
//             double dR_bar_AO_cube=pow((R_v_AO(j,k)-R_u_AO(j,k)),3);
//             A_bar_A_O(j,k)=2/dR_bar_AO_cube;
//             B_bar_A_O(j,k)=-3*(R_v_AO(j,k)+R_u_AO(j,k))/dR_bar_AO_cube;
//             C_bar_A_O(j,k)=6*R_v_AO(j,k)*R_u_AO(j,k)/dR_bar_AO_cube;
//             D_bar_A_O(j,k)=pow(R_v_AO(j,k),2)*(R_v_AO(j,k)-3*R_u_AO(j,k))/dR_bar_AO_cube;
//         }
        
//     }
    
// }


// mat A_D_O, B_D_O, C_D_O, D_D_O;
// mat R_m_DO, R_bar_DO, R_u_DO;
double R_m_DD, R_m_DDO, R_m2_DD, R_u_DD, R_bar_DD;
double A_D_D, B_D_D, C_D_D, D_D_D;

void calVfield_defenders(int ND) {
    double E_bar_O1=10*rho_D;  //repulsive
    double E_bar_O2=20*rho_D;  //blending of repulsive and attractive
    double E_bar_O3=5;

    int default_value = 1;
 
    // resize the vector to `NA` elements of type std::vector<int>,
    // each having size `NO` and default value
    // NA*NO matrix
    // R_m_DO.resize(ND,NO);
    // R_bar_DO.resize(ND,NO);
    // R_u_DO.resize(ND,NO);
    // A_D_O.resize(ND,NO);
    // B_D_O.resize(ND,NO);
    // C_D_O.resize(ND,NO);
    // D_D_O.resize(ND,NO);
    
    // for (size_t j = 0; j < ND; j++)
    // {
    //     for (size_t k = 0; k < NO; k++)
    //     {
    //         double rho_bar=2*10*rho_D;
    //         R_m_DO(j,k)=0;
    //         R_bar_DO(j,k)=2*30*rho_D;
    //         R_u_DO(j,k)=2*40*rho_D;
    //         double dR_DO_cube=pow((R_u_DO(j,k)-R_bar_DO(j,k)),3);
    //         A_D_O(j,k)=2/dR_DO_cube;
    //         B_D_O(j,k)=-3*(R_u_DO(j,k)+R_bar_DO(j,k))/dR_DO_cube;
    //         C_D_O(j,k)=6*R_u_DO(j,k)*R_bar_DO(j,k)/dR_DO_cube;
    //         D_D_O(j,k)=pow(R_u_DO(j,k),2)*(R_u_DO(j,k)-3*R_bar_DO(j,k))/dR_DO_cube;
            
    //     }
        
    // }

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

    // cout << "A_D_D: " <<A_D_D << endl;
    // cout << "B_D_D: " <<B_D_D << endl;
    // cout << "C_D_D: " <<C_D_D << endl;
    // cout << "D_D_D: " <<D_D_D << endl;
}

//Initialize the defenders
// double rD0[2] = {-150,-300};
vec rD0(2);
mat rD;
arma::mat XD0;
mat rSD_goal;
mat XD;
void InitializeDefenders(int NA, int ND) {
    rD0 = {-150,-300};
    // cout << "start initialize defenders" << endl;
    rD.resize(2,ND);
    rD = {{-204, -60, -800, 600, -150, -330, 340, 11, -100, -200, -500, 200, 300,  -120, -240, 100, -300, -400},
          {-250,-365, 70  , 250, 500 , 23,   -10, 11, -100, 300,  100,  150, -130, 120,  70,   120, -30,  130}};

    XD0 = join_cols(rD, arma::zeros(2,ND));
    // XD0.print("XD0: ");
    // cout << "finish initialized defenders" << endl;
    
}

void calAllparametersExperiment(int NA, int ND) {
 
    // Function declared in header
    // file to find the sum
    calControlLimits();
    VectorFields_A();
    initialAttackersVel(NA);
    cal(NA);
    InitializeAttackers(NA);
    defenders(NA,ND);
    // FormationOri();
    calVfieldA();
    fillR(ND);
    // calVfield_attackers(NA);
    calVfield_defenders(ND);
    InitializeDefenders(NA, ND);
}

// int main(){
//     calControlLimits();
//     VectorFields_A();
//     int NA = 18;
//     initialAttackersVel(NA);
//     cal(NA);
//     InitializeAttackers(NA);
//     int ND = NA;
//     defenders(NA, ND);
//     calVfieldA();
//     fillR(ND);
//     calVfield_defenders(ND);
//     InitializeDefenders(NA, ND);

//     // rA.print("rA: ");
//     // vA.print("vA: ");
// }