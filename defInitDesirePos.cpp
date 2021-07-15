#include <iostream>
#include <armadillo>
#include "AllParametersExperiment.cpp"
#include "findShortestPath.cpp"

using namespace std;
using namespace arma;
struct defDesForm
{
    mat rDFc0;
    mat XD_des0;
    mat XD_des_dot0;
    double phi;
    double phi_dot;
};

struct path_elem
{
    mat nodes;
    mat rV;
    double P;
    mat S;
    mat rVC;
    int segType;
    int NS;
};

struct tanG_prime_elem
{
    mat G;
    mat G_pathType;
    cube rVO;
    cube rVO2;
    mat gTO_all;
    mat rTO_all;
    mat obsld_all;
    mat obsVertId_all;
    mat obsVertPos_all;
};

struct interSec_elem
{
    int Flag;
    int Pbar1;
    int Pbar2;
    int flag0;
};

struct pathVel_elem
{
    mat v;
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

struct motionPlan
{
    mat assign;
    double assignCost;
    double maxOptT;
    std::vector<path_elem> path;
    std::vector<tanG_prime_elem> tanG_prime;
    std::vector<::vector<interSec_elem>> interSec;
    std::vector<pathVel_elem> pathVel;
    mat optT;
    mat leadTime;
    mat startTime;
    double tCompMIQP;
    double tCompMILP;
    double tComp;
    mat XD;
    mat XD_des;
};

struct DesiredPos {
    defDesForm dDf;
    motionPlan mP;
};

DesiredPos defInitDesiredPos(mat XA0, mat XD, int Na, int ND, double RD0, double v_maxA, double DeltaT_s) {
    colvec rAcm = arma::sum(XA0.submat(0,0,1,XA0.n_cols-1),1) / Na;
    colvec vAcm = arma::sum(XA0.submat(2,0,3,XA0.n_cols-1),1) / Na;
    double thetaAcm0=atan2(rAcm(2)-rP(2),rAcm(1)-rP(1))+20*M_PI/180;
    if (thetaAcm0 < 0) {
        thetaAcm0 = thetaAcm0 + 2 * M_PI;
    }
    






    
}