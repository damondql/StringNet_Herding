// #include "AllParametersExperiment.cpp"
#include "symDerivative.cpp"
#include "findCommGraphAndFormDist.cpp"
#include "defInitDesirePos.cpp"
// #include "findCoordOnPath.cpp"
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
            XA(i,j) = XA0(i,j);
        }
        
    }
    for (size_t i = 0; i < XD0.n_rows; i++)
    {
        for (size_t j = 0; j < XD0.n_cols; j++)
        {
            XD(i,j) = XD0(i,j);
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
    epsilon=M_PI/100;
    RAD=RAD_max;
    RA0_des=0;
    Rii0 = RA0*sqrt(2*(1-cos(2*M_PI/NA)));
    Rii1=Rii0*sqrt(2*(1-cos(180-2*M_PI/NA)));
    CommGraph attacker_graph = findCommGraphAndFormDist(NA,1,RA);
    CommGraph defender_graph_close = findCommGraphAndFormDist(ND+1, 2, rho_sn);
    // double RDF_open = 0;
    if (flagExp == 1) {
        RDF_open = 4.4;
    } else if (flagGazeboExp == 1)
    {
        RDF_open = 1.5 * rho_sn;
    } else {
        RDF_open = rho_sn + 20;
    }
    CommGraph defender_graph_open = findCommGraphAndFormDist(ND+1, 4, RDF_open);
    R_DD_string = 1.5 * defender_graph_close.Rij_tilde(0,1);
    RDF_closed = 1.1 * rho_sn;
    WDString = arma::zeros<mat>(ND, ND);
    // cout  << "RDF_closed: " << RDF_closed << endl;
}

DesiredPos motionP_result;

void getMoitonPlan(){
    motionP_result =  defInitDesiredPos(YA, XD0, NA, ND, RDF_open, v_maxA[0], 105);
    motionP_result.dDf.phi += M_PI;
    motionP_result.mP.startTime += ones<mat>(ND,1) * 25;
    // XD.print("XD:");
    arma::mat tempM;
    tempM.reshape(XD.n_rows, XD.n_cols+1);
    tempM.submat(0,0,XD.n_rows-1,XD.n_cols-1) = XD;
    tempM.submat(0,tempM.n_cols-1,motionP_result.dDf.rDFc0.n_rows-1, tempM.n_cols-1) = motionP_result.dDf.rDFc0;
    tempM.submat(motionP_result.dDf.rDFc0.n_rows,tempM.n_cols-1, tempM.n_rows-1,tempM.n_cols-1) = zeros<mat>(2,1);
    XD = tempM;
    // XD.print("XD:");
    tempM.resize(X.n_rows,1);
    tempM.submat(0,0,XA.n_rows-1,0) = XA.as_col();
    tempM.submat(XA.n_rows,0,tempM.n_rows-1,0) = XD.as_col();
    // tempM.print("tempM:");
    X.col(0) = tempM;
    // X.submat(0,0,19,9).print("X sub:");

}

void calDistance(){
    vec indDef = regspace(1,ND+1);
    int na = motionP_result.mP.assign.n_elem;
    int nid = indDef.n_elem;
    vec tempV(indDef.n_elem);
    int indx_num = na + regspace(na+1,nid).n_elem;
    vec indx(indx_num);
    indx.subvec(0,na-1) = motionP_result.mP.assign;
    indx.subvec(na,indx_num-1) = regspace(na+1,nid);
    // indx.print("indx:");
    for (int i = 0; i < tempV.n_elem; i++)
    {
        tempV(indx(i)-1) = indDef(i);
    }
    indDef = tempV;
    // cout << "na: " << na <<endl;
    // cout<< "nid: " <<nid <<endl;
    vec RDjDl(ND);
    if (ND >1)
    {
        for (int j = 0; j < ND-1; j++)
        {
            for (int i = j+1 ; i < ND; i++)
            {
                
            }
            
        }
        
    }
    


}


int main() {
    calAllparametersExperiment();
    AllocateMemory();
    measurements(2);
    initial_contorl();
    // arma::vec rAcm = {5.9807,-33.2064};
    // arma::vec rP = {14,-1};
    // path_elem p = findShortestPath(rAcm, rP);
    // p.rV.print("rV: ");
    // p.S.print("S: ");
    // double q = 16.0949;
    // CoorOnPath A = findCoordOnPath(q, p);
    // A.rp.print("rDFc0: ");
    // A.thetap.print("thetaAcom0: ");
    // v_maxD.print("v_maxD: ");
    // v_maxDC.print("v_maxDC: ");
    // u_maxD.print("u_maxD: ");
    // cout << "rho_P:" << rho_P << endl;
    // cout << "rho_safe:" << rho_safe << endl;
    // cout << "rho_sn:" << rho_sn << endl;
    // cout << "rho_Acon:" << rho_Acon << endl;
    // YA.print("YA:");
    YA(0,0) = 6.15469799011143;
    YA(1,0) = -33.3348469481075;
    YA(2,0) = -0.116564575844461;
    YA(3,0) = 0.207257244260275;
    YA.print("YA: ");
    XD0 = {{11.9907000000000,9.37240000000000,16.1838000000000},
{-10.7580000000000,-7.12780000000000,-5.84110000000000},
{0,0,0},
{0,0,0}};
//     XD0.print("XD0: ");
    // cout << "RDF_open:" << RDF_open << endl;
    getMoitonPlan();
    calDistance();
}