#include "AllParametersExperiment.cpp"
#include "symDerivative.cpp"
#include "findCommGraphAndFormDist.cpp"
#include <armadillo>

using namespace std;
using namespace arma;

int flagExp = 1;
int flagGazeboExp = 0;


auto t = 0;
auto Time = 0;
auto t_max = 2000;
auto endGame = 0;
double Niter;
mat vA_dot(2,NA, fill::zeros);
mat X;
mat XA(4, NA, fill::zeros);
//Allocate memory
void AllocateMemory() {
    Niter = floor(t_max/dt);
    X.set_size(4*(ND+NA+1), Niter+1);
    X.fill(0);
    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < NA; j++)
        {
            XA(i,j) = XA0[i][j];
        }
        
    }
    for (size_t i = 0; i < XD0.size(); i++)
    {
        for (size_t j = 0; j < XD0[0].size(); j++)
        {
            XD(i,j) = XD0[i][j];
        }
        
    }
    
}

mat Cov_YA(4,4,fill::zeros);
mat YA;
mat SD;
mat SD_arr;
mat rAcm;
mat rDcm;
mat vDcm;
double flagNDeven;
void measurements(int s) {
    if(s == 1 || s == 3) {
        Cov_YA = {{0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0}};
    } else if (s == 2){ 
        Cov_YA = {{pow(0.9/3,2),0,0,0},
                  {0,pow(0.9/3,2),0,0},
                  {0,0,pow(0.12,2),0},
                  {0,0,0,pow(0.12,2)}};
    }
    YA = XA + arma::mvnrnd(arma::zeros(4,1), Cov_YA, NA);
    SD = arma::zeros(ND,1);
    SD_arr = SD;
    rAcm = arma::sum(YA,1)/NA;
    rAcm = rAcm(2,0, arma::size(2,1));
    rDcm = arma::sum(XD,1)/ND;
    vDcm = rDcm(2,0,arma::size(2,1));
    rDcm = rDcm(0,0, arma::size(2,1));
    flagNDeven = int(ND) % 2;
}


//Initialization of vectors and matrices
//Initial control
mat uA;
mat RPA_arr,RAD_arr;
mat Psi, Psi_prime, Psi_des, Psi_des_dot, Psi_des_ddot, Psi_prime_dot, Psi_prime_ddot;
mat vA_dot_arr, FA_dot_arr, XD_Des, U, times, minEAO, minEDO, minRDD, minRAD, minRAA, arr_norm_F_AO;
mat AOut;
double NAout = 0;
double AIcount = 0; //counter for Attacker interceptions;
double safe_flag=0;
double RPA, epsilon, RAD, RA0_des, Rii0, Rii1;
double RDF_open = 0;
double R_DD_string, RDF_closed;
// Initial string chain
int flagHerd = 0;
int flagDForm = 0;
mat WDString;
void initial_contorl() {
    uA = vA;
    RPA_arr = arma::zeros(1, Niter+1);
    RAD_arr=arma::zeros<mat>(1,Niter+1);
    Psi=arma::zeros<mat>(1,Niter+1);
    Psi_prime=arma::zeros<mat>(1,Niter+1);
    Psi_des=arma::zeros<mat>(1,Niter+1);
    Psi_des_dot=arma::zeros<mat>(1,Niter+1);
    Psi_des_ddot=arma::zeros<mat>(1,Niter+1);
    Psi_prime=arma::zeros<mat>(1,Niter+1);
    Psi_prime_dot=arma::zeros<mat>(1,Niter+1);
    Psi_prime_ddot=arma::zeros<mat>(1,Niter+1);
    vA_dot_arr=arma::zeros<mat>(2*NA,Niter+1);
    FA_dot_arr=arma::zeros<mat>(2*NA,Niter+1);
    XD_Des=arma::zeros<mat>(4*ND,Niter+1);
    U=arma::zeros<mat>(2*(ND+NA+1),Niter+1);
    times=arma::zeros<mat>(1,Niter+1);
    minEAO=arma::zeros<mat>(1,Niter+1);
    minEDO=minEAO;
    minRDD=minEAO;
    minRAD=minEAO;
    minRAA=minRAD;
    arr_norm_F_AO=minEAO;
    std::complex<double> comp = {rA(0,0)-rP[0],rA(1,0) - rP[1]};
    RPA=sqrt(std::norm(comp));
    cout << "RPA: " << RPA << endl;
    epsilon=M_PI/100;
    RAD=RAD_max;
    RA0_des=0;
    Rii0 = RA0*sqrt(2*(1-cos(2*M_PI/NA)));
    Rii1=Rii0*sqrt(2*(1-cos(180-2*M_PI/NA)));
    output attacker_graph = findCommGraphAndFormDist(NA,1,RA);
    output defender_graph_close = findCommGraphAndFormDist(ND+1, 2, rho_sn);
    double RDF_open = 0;
    if (flagExp == 1) {
        RDF_open = 4.4;
    } else if (flagGazeboExp == 1)
    {
        RDF_open = 1.5 * rho_sn;
    } else {
        RDF_open = rho_sn + 20;
    }
    output defender_graph_open = findCommGraphAndFormDist(ND+1, 4, RDF_open);
    R_DD_string = 1.5 * defender_graph_close.Rij_tilde(0,1);
    RDF_closed = 1.1 * rho_sn;
    WDString = arma::zeros<mat>(ND, ND);
}



int main() {
    calAllparametersExperiment();
    AllocateMemory();
    measurements(2);
    initial_contorl();
}