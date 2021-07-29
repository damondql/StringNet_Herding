// #include "AllParametersExperiment.cpp"
#include "symDerivative.cpp"
#include "findCommGraphAndFormDist.cpp"
#include "defInitDesirePos.cpp"
#include "controlAttacker4.cpp"
#include "controlDefender5.cpp"
// #include "findCoordOnPath.cpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <istream>


using namespace std;
using namespace arma;

int flagExp = 1;
int flagGazeboExp = 0;


double t = 0;
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
mat vAcm;
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
    YA.load("../../../../../Downloads/swarm_matlab/controlD/YA.txt");
    // XD0.load("../../../../../Downloads/swarm_matlab/controlD/XD0.txt");
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
CommGraph attacker_graph;
CommGraph defender_graph_close;
CommGraph defender_graph_open;
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
    attacker_graph = findCommGraphAndFormDist(NA,1,RA);
    defender_graph_close = findCommGraphAndFormDist(ND+1, 2, rho_sn);
    // double RDF_open = 0;
    if (flagExp == 1) {
        RDF_open = 4.4;
    } else if (flagGazeboExp == 1)
    {
        RDF_open = 1.5 * rho_sn;
    } else {
        RDF_open = rho_sn + 20;
    }
    defender_graph_open = findCommGraphAndFormDist(ND+1, 4, RDF_open);
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

vec indDef;
void calDistance(){
    //Defenders indices for the formation
    indDef = regspace(1,ND+1);
    int na = motionP_result.mP.assign.n_elem;
    int nid = indDef.n_elem;
    vec tempV(indDef.n_elem);
    int indx_num = na + regspace(na+1,nid).n_elem;
    vec indx(indx_num);
    indx.subvec(0,na-1) = motionP_result.mP.assign;
    indx.subvec(na,indx_num-1) = regspace(na+1,nid);
    
    for (int i = 0; i < tempV.n_elem; i++)
    {
        tempV(indx(i)-1) = indDef(i);
    }
    indDef = tempV;
    

    //minimum distance between the defenders
    vec RDjDl(ND);
    vec arr_minRDD(ND-1);
    vec arr_minRAD(ND-1);
    if (ND >1)
    {
        for (int j = 0; j < ND-1; j++)
        {
            for (int i = j+1 ; i < ND; i++)
            {
                RDjDl(i) = norm(rD.col(j) - rD.col(i));
            }
            arr_minRDD(j) = RDjDl.subvec(j+1,RDjDl.n_rows-1).min();
            arr_minRAD(j) = norm(rD.col(j) - rA);
        }
        minRDD(0,0) = arr_minRDD.min();
        minRAD(0,0) = arr_minRAD.min();
        
    }
    

}

double dpsiAi0=0.02;
double vA_ref_max=4;
extern double dpsi_var;
int flagDefReachOpen=0;
int flagDefReachClosed=0;
int flagDefConnect=0;
int flagAttInSight=0;
mat defReachCount(ND,1, fill::zeros);
mat flagAttInObs(NA,NO, fill::zeros);
mat FlagDefInObs(ND,NO, fill::zeros);
mat betaAv0(NA,NO,fill::zeros);
mat betaDv0(ND,NO,fill::zeros);
mat mDv0=betaDv0;
mat cDv0=mDv0;
mat mAv0=betaAv0;
mat cAv0=mAv0;
mat speedA0=betaAv0;
mat assignment; // needs to be initialized later
int distTol=3;
int countDDes=0;

int flagGather=1;
int flagSeek=0;
int flagEnclose=0;
int flagAttackerStayTogether=1;

void checkFormation(){
    assignment = regspace(1,1,ND);
    mat RefTraj(4*ND, Niter);
    mat SigmaProdD_arr;
    mat rAcm_arr;
    mat vAcm_arr;
    for (int ti = 0; ti < 1; ti++)
    {
        mat Psi(NA,1, fill::zeros);
        mat Psi_dot;
        mat Psi_ddot;
        Psi_dot = Psi;
        Psi_ddot = Psi;
        
        int countAinS = 0;
        for (int ii = 0; ii < NA; ii++)
        {
            if (norm(rA.col(ii) - rS) < rho_S)
            {
                countAinS++;
            }
        }

        int countDinS = 0;
        for (int jj = 0; jj < ND; jj++)
        {
            if (norm(rD.col(jj) - rS) < rho_S)
            {
                countDinS++;
            }
        }
        if (countAinS > NA-1 && countDinS > ND-1)
        {
            break;
        }
        //////////////////////////////////////////
        /////////   safe flag protion     ////////
        //////////////////////////////////////////

        mat xd0(XD.n_rows, ND);
        for (int i = 0; i < ND; i++)
        {
            xd0.col(i) = XD.col(i);
        }
        XD_Des.col(ti+1) = xd0.as_col();
        
        vec pairDA(ND);
        pairDA = regspace(1,1,ND);
        double kref = 0.1;

        mat refTraj(4,ND, fill::ones);
        RefTraj.col(ti) = refTraj.as_col();


        vec XA_lead_goal(4,fill::zeros);
        XA_lead_goal.subvec(0,1) = rP;
        mat XA_goal(XA_lead_goal.n_rows, NA);
        mat XA_goal_dot(XA_lead_goal.n_rows, NA);
        XA_goal.col(0) = XA_lead_goal;
        XA_goal_dot.col(0) = XA_lead_goal;

        control_attacker_t control_A_result = controlAttacker4(XA,XA_goal,XA_goal_dot,flagEnclose, flagHerd, XD,attacker_graph.W,WDString,NA,ND);
        SigmaProdD_arr.insert_cols(ti,control_A_result.SigmaProdD);

        rAcm = arma::sum(XA.submat(0,0, 1,XA.n_cols-1),1)/NA;
        vAcm = arma::sum(XA.submat(2,0, 3,XA.n_cols-1),1)/NA;
        
        WDString = zeros<mat>(ND,ND);
        
        for (int j = 1; j <= ND-1; j++)
        {
            uvec j11 = find(motionP_result.mP.assign == j);
            uvec j22 = find(motionP_result.mP.assign == j+1);
            int j1 = j11(0);
            int j2 = j22(0);
            if (arma::norm(XD.submat(0,j1,1,j1) - XD.submat(0,j2,1,j2)) <= R_DD_string)
            {
                WDString(j1,j2) = 1;
                WDString(j2,j1) = 1;
            } else {
                WDString(j1,j2) = 0;
                WDString(j2,j1) = 0;
            } 
        }
        mat uD;
        if (flagGather == 1 && flagSeek != 1 && flagEnclose != 1 && flagHerd !=1)
        {
            mat XD_des = motionP_result.dDf.XD_des0;
            mat XD_des_dot = motionP_result.dDf.XD_des_dot0;
            uD = controlDefender5(XD, SD, regspace(1,1,ND+1), motionP_result.mP.assign, XD_des, XD_des_dot, motionP_result.mP, times(ti), ND);
            for (int j = 0; j < count; j++)
            {
                if (arma::norm(XD.submat(0,indDef(j), 1,indDef(j))- XD_des(0,j,1,j)) < 1e-3) {
                    defReachCount(j) = 1;
                    if (accu(defReachCount) >= ND)
                    {
                        flagDefReachOpen = 1;
                        flagSeek = 1;
                        flagGather = 0;
                    }
                    
                }
            }
            extern double ti_g = ti;
        } else if (flagGather != 1 && flagSeek != 1 && flagEnclose == 1 && flagHerd != 1) {
            if (flagDForm != 1)
            {
                flagDForm = 1;
            }
            for (int j = 0; i < ND; j++)
            {
                double RD0 = RDF_closed;
                double thetaD;
                thetaD = motionP_result.dDf.phi + 2*M_PI *(j+1)/ND - M_PI / ND;
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
    
//     XD0.print("XD0: ");
    // cout << "RDF_open:" << RDF_open << endl;
    // YA.load("../../../../../Downloads/swarm_matlab/controlD/YA.txt");
    // XD0.load("../../../../../Downloads/swarm_matlab/controlD/XD0.txt");
    getMoitonPlan();
    calDistance();
    checkFormation();
    // mat XD_des = motionP_result.dDf.XD_des0;
    // mat XD_des_dot = motionP_result.dDf.XD_des_dot0;
    // std::ifstream fin("../../../../../Downloads/swarm_matlab/controlD/t.txt");
    // double time;
    // fin >> time;
    // // time = 10;
    // std::cout << "t: "<<time << std::endl;

    // fin.close();
    // XD.load("../../../../../Downloads/swarm_matlab/controlD/XD.txt");
    // XD.print("XD: ");
    // SD.load("../../../../../Downloads/swarm_matlab/controlD/SD.txt");
    // SD.print("SD: ");
    // // indDef.print("indDef: ");
    // motionP_result.mP.assign.print("assign: ");
    // XD_des.print("XD_des: ");
    // XD_des_dot.print("XD_des_dot: ");
    // mat uD = controlDefender5(XD,SD, regspace(1,1,4), motionP_result.mP.assign, XD_des, XD_des_dot, motionP_result.mP, time,3);
    // uD.print("uD: ");

}