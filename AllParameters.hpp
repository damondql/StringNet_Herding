#pragma once

#include <math.h>
#include "iostream"
#include <math.h>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <armadillo>

#include <string>
#include <vector>
#include <list>
#include <deque>
#include <iterator>
#include "boost/bind.hpp"
#include "boost/function.hpp"


// Plotting
/* SetPlotDefaults;
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

global options

options = optimset('Display','off','MaxIter',1000,'TolFun',1e-3);
*/

// General
extern double dt, Rs,alpha, rho_B, rho_c_A ,rho_P, rho_D, rho_A;
extern double rho_P;    //Radius of the protected area


//Initial Positions of the protected and safe area
extern arma::vec rP;
extern arma::mat rS;
extern int NS;
extern double rho_S;

// Attackers

// Parameters for the attackers
extern double rho_A;
extern double rho_c_A;

// Control Limits
extern double C_d;
extern double vma;

extern double uma, rhoA_safe;

//Parameters for the vector fields for the attacker corresponding to the
//other attackers
extern double R_m_AA, R_bar_AA, R_u_AA;
extern double A_A_A, B_A_A, C_A_A, D_A_A;

//Initialize the attackers
extern int NA, ND;    //number of attackers and defenders
extern arma::vec v_maxA, u_maxA;
extern double NA_sep;

extern int axS[2];
extern double rDmin;

extern double RA0, rho_Acon;
extern arma::vec Rii00, Rik00;
extern arma::vec Rjk00;

extern int dRA0;
extern double rho_sn;

// Defenders
extern arma::vec v_maxD,v_maxDC,u_maxD,u_maxD1,u_maxD2,u_maxDr1,u_maxDr2;
extern double rho_safe,dthetai,alphaD_v;
extern double rhoAD_safe,rhoD_safe;
extern double umd1,umd2,umdf_s1,umdf_s2,umd,umdf_h1,vmd,vmdc,vmdf_s,umd_e1,umd_e2;
extern double rho_sn;
struct obstacle
{
    double rho_safe;
    int NO;
    std::vector<std::vector<std::vector<double>>> rVO;
    std::vector<double> rVO1;
    std::vector<std::vector<double>> rCO2;
};

extern obstacle obs;
extern int largeP;
extern double kr0, kr1, kr2, kv1, kv2, RD_con;

//control gains for attackers
extern double  kADr,kAFr, kAOr, kAOr2, kAPr, kAPv, kAFv, alphaAFv, kAOv, kAOv2;
extern double  alphaAOv, kADv, alphaADv;

//finite time gains for the defenders
extern double kDOr, kDDr, kDOv1, kDOv2, alphaDOv, kDDv, alphaDDv;
extern double kDFr,  kDFr2, kextern, doubleDFv, alphaDFr, alphaDFv, kDRr, kDRv;
extern double kDFphi, kDFphiv, kDFphir, kDFphid;
extern double kDDesr, kDDesv;

//bound on the convergence error during tracking


///////////////////////////////////
//Convex polygonal Obstalces///////
///////////////////////////////////

extern arma::cube rVO;

//For formation orientation
extern arma::vec R_bar_AcOc, R_u_AcOc, A_Ac_Oc, B_Ac_Oc, C_Ac_Oc, D_Ac_Oc;

//Parameters for the vector fields for the attacker corresponding to the
//defenders
extern double R_m_AD,R_bar_AD,R_u_AD;
extern double global,A_A_D,B_A_D,C_A_D,D_A_D;

extern double R_m_AD2,R_bar_AD2,R_u_AD2;
extern double A_A_D2,B_A_D2,C_A_D2,D_A_D2;

extern arma::vec Rij0,Rjj0;
extern double RAD_max;
extern  double nO, A, B, C, D, aO, bO, E_bar_O, E_m_O, E_u_O, GO;

std::vector<double> w, h, w_bar, h_bar;

//Parameters for the vector fields for the attacker corresponding to the
//obstacles
extern arma::mat R_m_AO, R_bar_AO, R_u_AO, R_v_AO;
extern arma::mat A_A_O, B_A_O, C_A_O, D_A_O;
extern arma::mat A_bar_A_O, B_bar_A_O, C_bar_A_O, D_bar_A_O;

extern arma::mat A_D_O, B_D_O, C_D_O, D_D_O;
extern arma::mat R_m_DO, R_bar_DO, R_u_DO;
extern double R_m_DD, R_m_DDO, R_m2_DD, R_u_DD, R_bar_DD;
extern double A_D_D, B_D_D, C_D_D, D_D_D;

extern arma::mat rSD_goal;



void calControlLimits();
void VectorFields_A();
void initialAttackersVel(int NA);
void cal();
void InitializeAttackers(int NA);
void defenders(int NA, int ND);
void FormationOri();
void calVfieldA();
void fillR(int NA);
void calVfield_formation();
void calVfield_attackers(int NA);
void calVfield_defenders(int ND);
void InitializeDefenders();
void calAllparametersExperiment(int NA, int ND);