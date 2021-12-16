#include "AllParameters.cpp"
// #include "symDerivative.cpp"
#include "findCommGraphAndFormDist.cpp"
#include "defInitDesirePos.cpp"
#include "controlAttacker4.cpp"
#include "controlDefender5.cpp"
#include "defDesiredOpenForm.cpp"
// #include "controlFiniteTimeTrajTracking.cpp"
#include "inhull.cpp"
#include "controlDefenderFormation4.cpp"
#include "defDesiredClosedForm.cpp"
#include "modifiedDIDynamics.cpp"
#include "findCoordOnPath.cpp"
#include "assignDefenders2ClustersMIQCP.cpp"
#include "helperFunction.cpp"
#include "DBSCAN.cpp"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <istream>


using namespace std;
using namespace arma;

int flagExp = 0;
int flagGazeboExp = 0;


double t = 0;
auto Time = 0;
auto t_max = 2000;
auto endGame = 0;
double Niter;
mat vA_dot;
mat X;
mat XA;
//Allocate memory
void AllocateMemory(int NA, int ND) {
    vA_dot.resize(2,NA);
    vA_dot.fill(0);
    XA.resize(4,NA);
    Niter = floor(t_max/dt);
    X.set_size(4*(ND+NA), Niter+1);
    X.fill(0);
    XA = XA0;
    XD = XD0;
    // XA.print("XA: ");
    // XD.print("XD: ");
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
void measurements(int NA, int ND) {
    Cov_YA = {{pow(0.9/3,2),0,0,0},
                  {0,pow(0.9/3,2),0,0},
                  {0,0,pow(0.12,2),0},
                  {0,0,0,pow(0.12,2)}};
    // YA = XA + arma::mvnrnd(arma::zeros(4,1), Cov_YA, NA);
    // YA.load("../../../../../Downloads/swarm_matlab/controlD/YA.txt");
    YA = XA;
    // XD0.load("../../../../../Downloads/swarm_matlab/controlD/XD0.txt");
    SD = arma::zeros(ND,1);
    SD_arr = arma::zeros(ND,Niter+1);
    rAcm = arma::sum(YA,1)/NA;
    rAcm = rAcm(2,0, arma::size(2,1));
    rDcm = arma::sum(XD,1)/ND;
    vDcm = rDcm(2,0,arma::size(2,1));
    rDcm = rDcm(0,0, arma::size(2,1));
    // vDcm.print("vDcm: ");
    // rDcm.print("rDcm: ");
    // rAcm.print("rAcm: ");
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

int flagDForm = 0;
double epsilon_clust;
mat WDString;
CommGraph attacker_graph;
CommGraph defender_graph_close;
CommGraph defender_graph_open;
mat rhoA_con, rhoA_con_A, rho_SN;
int MinPts;
double sf_RDF_open;
double sf_RDF_closed;
void initial_contorl(int NA, int ND) {
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
    U=arma::zeros<mat>(2*(ND+NA),Niter+1);
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
    attacker_graph = findCommGraphAndFormDist(NA,1,RA0);
    defender_graph_close = findCommGraphAndFormDist(ND+1, 2, rho_sn);
    // attacker_graph.Rij_tilde.print("attacker R: ");
    // attacker_graph.W.print("attacker W: ");
    // defender_graph_close.Rij_tilde.print("defender close R: ");
    // defender_graph_close.W.print("defender close W: ");
    // double RDF_open = 0;
    if (flagExp == 1) {
        RDF_open = 2.2;
    } else if (flagGazeboExp == 1)
    {
        RDF_open = 1.5 * rho_sn;
    } else {
        RDF_open = 1.5 * rho_sn;
    }
    defender_graph_open = findCommGraphAndFormDist(ND+1, 4, RDF_open);
    // defender_graph_open.Rij_tilde.print("d open R: ");
    // defender_graph_open.W.print("d open W: ");
    R_DD_string = 80;
    // RDF_closed = 1.1 * rho_sn;
    WDString = arma::zeros<mat>(ND, ND);
    vec tempV(1);
    tempV(0) = rhoA_con_fun(NA,NA,ND,R_DD_string, rhoD_safe, R_m_DD);
    rhoA_con.insert_cols(rhoA_con.n_cols, tempV);
    rhoA_con_A.insert_cols(rhoA_con_A.n_cols,tempV);
    tempV(0) = (rho_sn_fun(ND,R_DD_string));
    rho_SN.insert_cols(rho_SN.n_cols,tempV);
    sf_RDF_open=1.05;
    sf_RDF_closed=1.1;
    RDF_open = sf_RDF_open * M_PI / 2 *rho_SN(0);
    RDF_closed=sf_RDF_closed*rho_SN(0);
    double coeff=0.9;
    MinPts=4;
    epsilon_clust=coeff*R_DD_string/2*(1/tan(M_PI/ND))*floor(MinPts/2)/(NA-1);

    epsilon_clust *= 1.3;
    // rhoA_con.print("rhoA_con: ");
    // rhoA_con_A.print("rhoA_con_A: ");
    // rho_SN.print("rho_SN: ");
    // cout  << "RDF_closed: " << RDF_closed << endl;
    // cout << "RDF_open: " <<RDF_open << endl;
    // cout << "epsilon_clust: " <<epsilon_clust << endl;
}

DesiredPos motionP_result;
std::vector<mat> XDFc;
vec assignment;
void getMotionPlan(int NA, int ND){
    XDFc.resize(ND);
    motionP_result =  defInitDesiredPos(YA, XD0, NA, ND, RDF_open, rhoA_con(0),v_maxA[0], 50);
    // cout << "got motion plan result " << endl;
    motionP_result.dDf.phi += M_PI;
    // motionP_result.mP.startTime += ones<mat>(ND,1) * 25;
    // XD.print("XD:");
    // arma::mat tempM;
    // tempM.reshape(XD.n_rows, XD.n_cols+1);
    // tempM.submat(0,0,XD.n_rows-1,XD.n_cols-1) = XD;
    // tempM.submat(0,tempM.n_cols-1,motionP_result.dDf.rDFc0.n_rows-1, tempM.n_cols-1) = motionP_result.dDf.rDFc0;
    // tempM.submat(motionP_result.dDf.rDFc0.n_rows,tempM.n_cols-1, tempM.n_rows-1,tempM.n_cols-1) = zeros<mat>(2,1);
    // XD = tempM;
    // // XD.print("XD:");
    // tempM.resize(X.n_rows,1);
    // tempM.submat(0,0,NA*XA.n_rows-1,0) = XA.as_col();
    // tempM.submat(NA*XA.n_rows,0,tempM.n_rows-1,0) = XD.as_col();
    // tempM.print("tempM:");
    X.col(0) = join_cols(XA.as_col(), XD.as_col());
    XDFc[0] = join_cols(motionP_result.dDf.rDFc0, zeros<mat>(2,1));
    vec tempV = regspace(1,1,ND);
    assignment.copy_size(tempV);
    for (int i = 0; i < tempV.n_elem; i++)
    {
        assignment(motionP_result.mP.assign(i)-1) = tempV(i);
    }
    assignment = assignment - 1;
    assignment.print("assignment:");
    // X.submat(0,0,19,9).print("X sub:");
}

vec indDef;
vec RDjDl;
void calDistance(int NA, int ND){
    RDjDl.resize(ND);
    //Defenders indices for the formation
    // indDef = regspace(1,ND+1);
    // int na = motionP_result.mP.assign.n_elem;
    // int nid = indDef.n_elem;
    // vec tempV(indDef.n_elem);
    // int indx_num = na + regspace(na+1,nid).n_elem;
    // vec indx(indx_num);
    // indx.subvec(0,na-1) = motionP_result.mP.assign;
    // indx.subvec(na,indx_num-1) = regspace(na+1,nid);
    
    // for (int i = 0; i < tempV.n_elem; i++)
    // {
    //     tempV(indx(i)-1) = indDef(i);
    // }
    // indDef = tempV;
    

    //minimum distance between the defenders
    vec arr_minRDD(ND-1);
    if (ND >1)
    {
        for (int j = 0; j < ND-1; j++)
        {
            for (int i = j+1 ; i < ND; i++)
            {
                RDjDl(i) = norm(rD.col(j) - rD.col(i));
            }
            arr_minRDD(j) = RDjDl.subvec(j+1,RDjDl.n_rows-1).min();
        }
        minRDD(0,0) = arr_minRDD.min();
        // minRAD(0,0) = arr_minRAD.min();
        // arr_minRDD.print("arr_minRDD: ");
    }

    mat arr_minRAD(NA, ND);
    for (int i = 0; i < NA; i++)
    {
        for (int j = i+1; j < ND; j++)
        {
            // cout <<"i: " << i << " j: " << j << endl;
            arr_minRAD(i,j) = norm(rD.col(j) - rA.col(i));
        }
        
    }
    minRAD(0,0) = arr_minRAD.min();
    
    // minRAD.print("minRAD: ");
    vec RAiAl(NA);
    vec arr_minRAA(NA-1);
    for (int i = 0; i < NA-1; i++)
    {
        for (int l = i+1; l<NA;l++)
        {
            RAiAl(l) = arma::norm(rA.col(i)-rA.col(l));
        }
        arr_minRAA(i) = RAiAl.subvec(i+1, RAiAl.n_rows-1).min();
    }
    minRAA(0,0) = arr_minRAA.min();
    // arr_minRAA.print("arr_minRAA: ");
    


    

}

double dpsiAi0=0.02;
double vA_ref_max=4;
extern double dpsi_var;
int flagDefReachOpen=0;
vec flagDefReachClosed(1,fill::zeros);
int flagDefConnect=0;
vec flagAttInSight(1,fill::zeros);
mat defReachCount;
mat flagAttInObs;
mat FlagDefInObs;
mat betaAv0;
mat betaDv0;
mat mDv0;
mat cDv0;
mat mAv0;
mat cAv0;
mat speedA0;
// mat assignment; // needs to be initialized later
int distTol=3;
int countDDes=0;

mat NClusterAD;
mat ClusterIdAD;
mat clusterIdAD;
vec flagClusterAEnclosed;
vec flagAEnclosed;
vec flagToSplitDefender;

vec NDinCluster(1);
vec NClusterD;
vec assignSafeArea;

vec flagSemiCircFormNotAchieved;
vec flagGather;
vec flagSeek;
vec flagEnclose;
vec flagHerd;

int flagAttackerStayTogether=1;
int splitCout = 0;

field<vec> ClusterRad;


void control_loop(int NA, int ND){
    cout<< "enter control loop" << endl;
    flagDefReachClosed.resize(1);
    flagDefReachClosed(0)=0;
    defReachCount = zeros<mat>(ND,1);
    NClusterAD = zeros(Niter+1,1);
    ClusterIdAD = zeros(NA, Niter+1);
    NClusterAD(0) = 1;
    clusterIdAD = zeros<mat>(NA,1);
    flagClusterAEnclosed = zeros<vec>(2);
    flagAEnclosed = zeros<vec>(NA);
    flagToSplitDefender.resize(1);
    flagToSplitDefender(0) = 0;
    NDinCluster(0) = ND;
    NClusterD.resize(1);
    NClusterD(0) = 1;
    // cout<< "111111111111111111111111" << endl;
    assignSafeArea = zeros<vec>(1);
    flagSemiCircFormNotAchieved.resize(1);
    flagSemiCircFormNotAchieved(0) = 0;
    flagGather.resize(1);
    flagGather(0) = 1;
    flagSeek.resize(1);
    flagSeek(0) = 0;
    flagEnclose = zeros<vec>(1);
    flagHerd = zeros<vec>(1);
    flagAttackerStayTogether = 1;
    splitCout = 0;
    ClusterRad.set_size(Niter,1);
    // cout<< "222222222222222222222222" << endl;
    mat RefTraj(4*ND, Niter);
    mat SigmaProdD_arr;
    mat rAcm_arr;
    mat vAcm_arr;
    mat sigmaProd;
    mat minRAAProjS(Niter,1);
    int out_ti;
    mat uD;
    mat XD_des;
    mat XD_des_dot;
    mat uDFc_trans;
    mat XDF_des;
    cube WDString_mat(ND,ND, Niter);
    vec RDF_OPEN;
    vec RDF_CLOSED;
    int bound = Niter;
    int countAinS = 0;
    int countAinP = 0;
    std::vector<mat> rAcm_new;
    std::vector<vec> indexOfNewClusterAD;
    // cout<< "33333333333333333333333" << endl;
    vec phi;
    vec phi_dot;
    vec ti_e;
    for (int ti = 0; ti < bound; ti++)
    {
        cout << "ti ===================== " << ti << endl;
        mat Psi = zeros<mat>(NA,1);
        mat Psi_dot = Psi;
        mat Psi_ddot = Psi;

        countAinS = 0;
        countAinP = 0;

        for (int ii = 0; ii < NA; ii++)
        {
            // cout << "clusterIdAD(ii): " << clusterIdAD(ii) << endl;
            // cout << "assignSafeArea(clusterIdAD(ii)-1): "<< assignSafeArea(clusterIdAD(ii)-1) << endl;
            if( arma::norm(rA.col(ii) - rS.col(assignSafeArea(clusterIdAD(ii)))) < rho_S )
            {
                countAinS++;
            }
            if(arma::norm(rA.col(ii) - rP ) < rho_P)
            {
                countAinP ++;
            }
            int countDinS = 0;
            for (int jj = 0; jj < ND; jj++)
            {
                if( arma::norm(rD.col(jj) - rS.col(assignSafeArea(clusterIdAD(ii))) )  < rho_S)
                {
                    countDinS++;
                }
            }
            // cout << "countAinS: " << countAinS << endl;
            // cout << "countAinP: " << countAinP << endl;
            // cout << "countDinS: " << countDinS << endl;
            if( countAinS > NA-1 || countDinS > ND-1)
            {
                break;
            }
            if(countAinP >= 1)
            {
                break;
            }

        }

        mat xd0 = XD;
        XD_Des.col(ti+1) = xd0.as_col();
        vec pairDA = regspace(1,1,ND);
        double kref = 0.1;
        mat refTraj = ones<mat>(4,ND);
        RefTraj.col(ti) = refTraj.as_col();

        mat XA_lead_goal = join_cols(rP, zeros<vec>(2));
        mat XA_goal = zeros<mat>(4,NA);
        uvec indAL = find(indALeader == 1);
        // indAL.print("indAL:");
        for (int i = 0; i < indAL.n_elem; i++)
        {
            XA_goal.col(indAL(i)) = XA_lead_goal;
        }
        // XA_goal.print("XA_goal:");
        mat XA_goal_dot = zeros<mat>(4,NA);
        // leaderIDA.print("leadIDA:");
        uA.print("uA: ");
        XA_goal.print("XA_goal");
        for (int i = 0; i < NA; i++)
        {
            if(indALeader(i) != 1)
            {
                XA_goal.submat(0,i,1,i) = XA.submat(0,leaderIDA(i),1,leaderIDA(i))+rA_follow.submat(0,i,1,i);
                XA_goal.submat(2,i,3,i) = XA.submat(2,leaderIDA(i),3,leaderIDA(i));
                XA_goal_dot.submat(0,i,1,i) = XA_goal.submat(2,i,3,i);
                XA_goal_dot.submat(2,i,3,i) = uA.submat(0,leaderIDA(i),1,leaderIDA(i))-C_d*arma::norm(XA.submat(2,leaderIDA(i),3,leaderIDA(i)))*XA.submat(2,leaderIDA(i),3,leaderIDA(i));
            }
        }
        // XA_goal.print("XA_goal: ");
        // XA_goal_dot.print("XA_goal_dot: ");
        vec clusterRad(NClusterAD(ti), fill::zeros);
        std::vector<uvec> indAinClusterAD(NClusterAD(ti));
        for (int c = 0; c < NClusterAD(ti); c++)
        {
            ti_e.resize(NClusterAD.n_elem);
            indAinClusterAD[c] = find(clusterIdAD == c);
            // indAinClusterAD[c].print("indAinClusterAD:");
            rAcm.resize(2,c+1);
            rAcm.col(c) = arma::sum(XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)) ,1) / indAinClusterAD[c].n_elem;
            for (int ii = 0; ii < indAinClusterAD[c].n_elem; ii++)
            {
                int i = indAinClusterAD[c](ii);
                double RAcA = arma::norm(rAcm.col(c) - XA.submat(0,i,1,i));
                if (RAcA > clusterRad(c))
                {
                    clusterRad(c) = RAcA;
                }
            }
            // clusterRad.print("clusterRad: ");
            
        }
        // rAcm.print("rAcm:");
        ClusterRad(ti) = clusterRad;

        int flagSpliAttackers = 0;
        int NClusterAD0 = NClusterAD(ti);

        for (int c = 0; c < NClusterAD(ti); c++)
        {
            cout << "enter for loop" << endl;
            if(c == 0){
                rAcm_new.clear(); //reset std::vector
            }
            if(!flagClusterAEnclosed(c))
            {
                if(clusterRad(c) > rhoA_con(c))
                {
                    flagSpliAttackers = 1;
                    if(flagSpliAttackers){
                        mat tempM1 = XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1));
                        tempM1.print("XA.submat:");
                        tempM1.save("../dbs/XA_submat.csv", csv_ascii);
                        vec clusterIdAD0 = DBSCAN(XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)), epsilon_clust, MinPts);
                        clusterIdAD0.print("clusterIdAD0:");
                        int NNewClusterAD = clusterIdAD0.max()+1;
                        if(NNewClusterAD > 0)
                        {
                            flagToSplitDefender(c) = 1;
                            splitCout++;
                        }
                        cout << "1.1" << endl;
                        // indAinClusterAD.clear();
                        std::vector<uvec> indAinClusterAD0;
                        mat tempM(2, NNewClusterAD);
                        for (int cc = 0; cc < NNewClusterAD; cc++)
                        {
                            indAinClusterAD0.push_back( indAinClusterAD[c](arma::find(clusterIdAD0 == cc)));
                            tempM.col(cc) = arma::sum(XA.submat(0,indAinClusterAD0[c](0),1,indAinClusterAD0[c](indAinClusterAD0[c].n_elem-1)),1)/indAinClusterAD0[cc].n_elem;
                        }
                        rAcm_new.push_back(tempM);
                        cout << "1.2" << endl;
                        uvec noiseId = arma::find(clusterIdAD0 == -1);
                        for (int no1 = 0; no1 < noiseId.n_elem; no1++)
                        {
                            int no = noiseId(no1);
                            double minDistNoiseClust = INFINITY;
                            int noc;
                            for (int cc = 0; cc < NNewClusterAD; cc++)
                            {
                                double disNoiseClust = arma::norm(XA.submat(0,no, 1,no) - rAcm_new[c].col(cc));
                                if(disNoiseClust < minDistNoiseClust)
                                {
                                    minDistNoiseClust = disNoiseClust;
                                    noc = cc;
                                }
                            }
                            clusterIdAD0(no) = noc;
                            indAinClusterAD[noc] = join_cols(indAinClusterAD0[noc], indAinClusterAD[c].subvec(no,no));
                            
                        }
                        //create new clusters of the attackers
                        cout << "1.3" << endl;
                        vec tempV(NNewClusterAD,fill::zeros);
                        assignSafeArea.resize(NNewClusterAD);
                        indAinClusterAD = indAinClusterAD0;
                        for (int cc = 0; cc < NNewClusterAD-1; cc++)
                        {
                            int c1 = NClusterAD0 + cc;
                            // cout << "cc: " << cc << endl;
                            // cout << "c1: " << c1 << endl;
                            // cout << "indAinClusterAD size: " << indAinClusterAD.size() << endl;
                            // cout << "indAinClusterAD0 size: " << indAinClusterAD0.size() << endl;
                            // // indAinClusterAD[c1].print("indAinClusterAD[c1]:");
                            // // indAinClusterAD0[cc+1].print("indAinClusterAD0[cc+1]: ");
                            // for(int a = 0; a < indAinClusterAD0.size(); a++) {
                            //     indAinClusterAD0[a].print("indAinClusterAD0[a]: ");
                            // }
                            // indAinClusterAD[c1] = indAinClusterAD0[cc+1];
                            // cout<<"get cc+1 " << endl;
                            for (int i = 0; i < indAinClusterAD[c1].n_elem; i++)
                            {
                                clusterIdAD(indAinClusterAD[c1](i)) = c1;
                            }
                            tempV(cc+1) = c1;
                            assignSafeArea(c1) = 0; //Initilization of the assigned safe area
                        }
                        indexOfNewClusterAD.push_back(tempV);
                        cout << "1.4" << endl;
                        cout << "indAinClusterAD size: " << indAinClusterAD.size() << endl;
                        cout << "indAinClusterAD0 size: " << indAinClusterAD0.size() << endl;
                        // cout << "ti: " << ti << endl;
                        indAinClusterAD[c] = indAinClusterAD0[0];
                        cout << "1.5" << endl;
                        for (int i = 0; i < indAinClusterAD[c].n_elem; i++)
                        {
                            clusterIdAD(indAinClusterAD[c](i)) = c;
                            cout << "1.6" << endl;
                        }
                        indexOfNewClusterAD[c](0) = c;
                        cout << "1.7" << endl;
                        NClusterAD0 = NClusterAD0+NNewClusterAD-1;
                        cout << "1.8" << endl;
                    }
                }    
            }
        }
        cout << "4444444444444444" << endl;
        NClusterAD(ti) = NClusterAD0; //Number of clusters at time ti after clustering
        arma::vec NAinClusterAD(NClusterAD(ti));
        if(flagSpliAttackers == 1)
        {
            RDF_OPEN.resize(NClusterAD(ti));
            RDF_CLOSED.resize(NClusterAD(ti));
            cout << "4.1" << endl;
            flagClusterAEnclosed = zeros<vec>(NClusterAD(ti));
            clusterRad = zeros<vec>(NClusterAD(ti));
            for (int c = 0; c < NClusterAD(ti); c++)
            {
                // clusterRad(c) = 0;
                indAinClusterAD[c] = arma::find(clusterIdAD == c);
                rAcm.resize(2,c+1);
                vAcm.resize(2,c+1);
                rAcm.col(c) = arma::sum(XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/indAinClusterAD[c].n_elem;
                vAcm.col(c) = arma::sum(XA.submat(2,indAinClusterAD[c](0),3,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/indAinClusterAD[c].n_elem;
                NAinClusterAD(c) = indAinClusterAD[c].n_elem;
                cout << "4.2" << endl;
                for (int ii = 0; ii < NAinClusterAD.n_elem; ii++)
                {
                    int i = indAinClusterAD[c](ii);
                    double RAcA =arma::norm(rAcm.col(c) - XA.submat(0,i,1,i));
                    if(RAcA > clusterRad(c))
                    {
                        clusterRad(c) = RAcA;
                    }
                }
                cout << "4.3" << endl;
                rhoA_con.resize(c+1);
                rho_SN.resize(c+1);
                RDF_OPEN.resize(c+1);
                RDF_CLOSED.resize(c+1);
                rhoA_con(c) = rhoA_con_fun(NAinClusterAD(c), NA, ND, R_DD_string, rhoD_safe, R_m_DD);
                cout << "4.4" << endl;
                rho_SN(c) = rhoA_con(c) + rhoD_safe+5;
                RDF_OPEN(c) = 1.2*M_PI/2*rho_SN(c);
                RDF_CLOSED(c) = 1.1*rho_SN(c);
                cout << "flagClusterAEnclosed size: " << flagClusterAEnclosed.size() << endl;
                cout << "c: " << c << endl;
                // flagClusterAEnclosed(c) = 0;
            }
        }
        cout << "5555555555555555555555" << endl;
        for (int c = 0; c < NClusterAD(ti); c++)
        {
            NAinClusterAD(c) = indAinClusterAD[c].n_elem;
            rAcm.resize(2,c+1);
            vAcm.resize(2,c+1);
            rAcm.col(c) = arma::sum(XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/ NAinClusterAD(c);
            vAcm.col(c) = arma::sum(XA.submat(2,indAinClusterAD[c](0),3,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/ NAinClusterAD(c);
        }
        cout << "666666666666666" << endl;
        ClusterIdAD.col(ti) = clusterIdAD.as_col();
        std::vector<vec> indDes;
        std::vector<vec> assign;
        if(flagGather(0) == 1)
        {
            cout << "flag Gather(0) == 1 if loop" << endl;
            NClusterD.resize(ti+1);
            NClusterD(ti) = 1;
            indDes.push_back(regspace(0,1,ND-1));
            assign.push_back(assignment);
            if(flagSpliAttackers == 1)
            {
                for (int cc = 1; cc < rAcm_new.size(); cc++)
                {
                    if(!rAcm_new[cc].is_empty())
                    {
                        rAcm_new[0].col(cc) = rAcm_new[cc].col(0);
                        int indc1 = rAcm_new[0].row(0).n_elem;
                        int indc2 = rAcm_new[cc].row(0).n_elem;
                        if(indc2 > 1)
                        {
                            vec temp1 = regspace(indc1,1,indc1+indc2-2);
                            vec temp2 = regspace(1,1,indc2-1);
                            for (int i = 0; i < temp1.n_elem; i++)
                            {
                                rAcm_new[0].col(temp1(i)) = rAcm_new[cc].col(temp2(i));
                                indexOfNewClusterAD[0](temp1(i))= indexOfNewClusterAD[cc](temp2(i));
                            }
                        }
                    }
                    
                }
                
            }
        } else 
        {
            cout << "flag Gather(0) == 1 else loop" << endl;
            NClusterD(ti) = NClusterD(ti-1);
            int NClusterD0 = NClusterD(ti);
            for (int c = 0; c < NClusterD(ti); c++)
            {
                if(flagToSplitDefender(c))
                {
                    vec tempV(indexOfNewClusterAD[c].n_elem);
                    for (int i = 0; i < tempV.n_elem; i++)
                    {
                        tempV(i) = NAinClusterAD(indexOfNewClusterAD[c](i));
                    }
                    
                    cout << "before new assign result" << endl;
                    new_assign_elem new_assign_result = assignDefenders2ClustersMIQCP(rAcm_new[c], rD, indDes[c].n_elem, tempV, assign[c], indDes[c], R_DD_string);
                    int NNewClusterD = new_assign_result.assign.size();
                    for (int cc = 1; cc < NNewClusterD; cc++)
                    {
                        int c1 = NClusterAD0 + cc -1;
                        flagGather(c1) = flagGather(c);
                        flagSeek(c1) = flagSeek(c);
                        flagEnclose(c1) = flagEnclose(c);
                        flagHerd(c1) = flagHerd(c);
                        flagAttInSight.resize(c1+1);
                        flagAttInSight(c1) = flagAttInSight(c);
                        NDinCluster.resize(c1+1);
                        NDinCluster(c1) = new_assign_result.assign[cc].n_elem;
                        indDes[c1] = new_assign_result.indDinSwarm[cc];
                        assign[c1] = new_assign_result.assign[cc];
                        double phi0 = atan2(rD(1, new_assign_result.assign[cc](NDinCluster(0))) - rD(1, new_assign_result.assign[cc](0)),rD(0, new_assign_result.assign[cc](NDinCluster(0))) - rD(0, new_assign_result.assign[cc](0))  );
                        phi.resize(c1+1);
                        phi_dot.resize(c1+1);
                        phi(c1) = phi0+M_PI/2;
                        phi_dot(c1) = 0;
                        rho_SN(c1) = rhoA_con(c1) + rhoD_safe + 5;
                        RDF_OPEN(c1) = sf_RDF_open * M_PI/2 * rho_SN(c1);
                        RDF_CLOSED(c1) = sf_RDF_closed * rho_SN(c1);
                        XDFc[c1] = (XD_des.col(new_assign_result.indDinSwarm[cc](0)) + XD_des.col(new_assign_result.indDinSwarm[cc](NDinCluster(c1))))/2;
                        flagSemiCircFormNotAchieved.resize(c+1);
                        flagSemiCircFormNotAchieved(c) = 0;
                        flagDefReachClosed(c) = 0;
                        NClusterD0 = NClusterD0+NNewClusterD-1;
                        flagToSplitDefender(c) = 0;
                    }
                    
                }
            }
            NClusterD(ti) = NClusterD0;
        }
        cout << "777777777777777" << endl;

        
        if(ti == 3000-1)
        {
            int NClusterA = 3;
            clusteridA = {{2,2,2,2,2,1,1,1,1,1,1,3,3,3,3,3,3,3}};
            clusteridA = clusteridA.t();
            clusteridA = clusteridA - 1;
            indAinClusterA.reset();
            indAinClusterA.set_size(3,1);
            mat tempM;
            tempM = {{9,10,11,6,7,8}};
            tempM = tempM.t();
            tempM = tempM -1;
            indAinClusterA(0) = tempM;
            tempM = {{4,5,1,2,3}};
            tempM = tempM.t();
            tempM = tempM -1;
            indAinClusterA(1) = tempM;
            tempM = {{15,12,13,14,16,17,18}};
            tempM = tempM.t();
            tempM = tempM -1;
            indAinClusterA(2) = tempM;
            indALeader = zeros<mat>(NA,1);
            rA_follow = zeros<mat>(2,NA);
            // cout << "attacker splitting initialized" << endl;
            if(NAinClusterA.n_elem < NClusterA) 
            {
                NAinClusterA.resize(NClusterA);
            }
            if(rhoA_con_A.size() < NClusterA) 
            {
                rhoA_con_A.resize(NClusterA);
            }
            for (int i = 0; i < NClusterA; i++)
            {
                // cout << "i: " << i << endl;
                // cout << "NAinClusterA: " <<NAinClusterA.size() << endl;
                // cout <<"rhoA_con_A size: " << rhoA_con_A.size() << endl;
                NAinClusterA(i) = indAinClusterA.n_elem;
                rhoA_con_A(i) = rhoA_con_fun(NAinClusterA(i), NA, ND, R_DD_string, rhoD_safe, R_m_DD);
                // cout << "rhoA_con_A: " << rhoA_con_A <<endl;
                indALeader(indAinClusterA(i,0)(0)) = 1;
                for (int j = 0; j < indAinClusterA(i,0).n_elem; j++)
                {
                    leaderIDA(indAinClusterA(i,0)(j)) = indAinClusterA(i,0)(0);
                }
                // cout << "finished setting leader" << endl;
                vec ind = indAinClusterA(i,0).submat(1,0,indAinClusterA(i,0).n_elem-1,0);
                // cout << "initialized ind" << endl;
                // ind.print("ind:");
                // cout << "indAinClusterA(i,0)(0): " << indAinClusterA(i,0)(0) << endl;
                for (int j = 0; j < ind.n_elem; j++)
                {
                    rA_follow.col(ind(j)) = rA.col(ind(j)) - rA.col(indAinClusterA(i,0)(0));
                }
                // cout << "finished rA_follow" << endl;
            }
        }
        // if (ti == bound-1) 
        // {
        //     cout << "controlAttacker4:" << endl;
        //     XA.print("XA: ");
        //     XA_goal.print("XA_goal: ");
        //     XA_goal_dot.print("XA_goal_dot:");
        //     leaderIDA.print("leaderIDA");
        //     XD.print("XD");
        //     attacker_graph.W.print("attacker_graph.W");
        //     attacker_graph.Rij_tilde.print("attacker_graph.Rij_tilde");
        //     WDString.print("WDString: ");
        //     clusteridA.print("clusteridA: ");
        //     NAinClusterA.print("NAinClusterA: ");
        //     rhoA_con_A.print("rhoA_con_A:");
        //     for(int i = 0; i < assign.size(); i++) {
        //         cout << "assign i: " << i << endl;
        //         assign[i].print("assign:");
        //     }
        // }
        control_attacker_t control_attacker_result = controlAttacker4(XA, XA_goal, XA_goal_dot, leaderIDA, XD, attacker_graph.W, attacker_graph.Rij_tilde, WDString, NA, ND, clusteridA, indAinClusterA, NAinClusterA, rhoA_con_A, flagAEnclosed, assign, ti, bound);
        control_attacker_result.uA.print("control_attacker_result uA:");
        uA = control_attacker_result.uA;
        cout << "8888888888888888888" << endl;
        uD = zeros<mat>(2, ND);
        WDString = zeros<mat>(ND, ND);
        rDcm.resize(2, NClusterD(ti));
        for (int c = 0; c < NClusterD(ti); c++)
        {
            for(int jj = 0; jj < indDes[c].n_elem -1; jj++)
            {
                int j1 = assign[c](jj);
                int j2 = assign[c](jj+1);
                if (arma::norm(XD.submat(0,j1,1,j1) - XD.submat(0,j2,1,j2)) <= R_DD_string)
                {
                    WDString(j1, j2) = 1;
                    WDString(j2, j1) = 1;
                } else
                {
                    WDString(j1, j2) = 0;
                    WDString(j2, j1) = 0;
                }
            }
            vec indD = assign[c];
            uvec indA = indAinClusterAD[c];

            for(int ii = 0; ii < indD.n_elem; ii++)
            {
                rDcm.col(c) = rDcm.col(c) + XD.submat(0,indD(ii),1,indD(ii));
            }
            rDcm.col(c) = rDcm.col(c) / indD.n_elem;
            cout << "9999999999999999999999" << endl;
            // flagGather.print("flagGather: ");
            // flagSeek.print("flagSeek: ");
            // flagEnclose.print("flagEnclose: ");
            // flagHerd.print("flagHerd: ");
            if(flagGather(c) == 1 &&  flagSeek(c) != 1 &&  flagEnclose(c) != 1 && flagHerd(c) != 1)
            {
                cout <<"gathering!" << endl;
                XD_des = motionP_result.dDf.XD_des0;
                XD_des_dot = motionP_result.dDf.XD_des_dot0;
                mat input_XD_des;
                mat input_XD_des_dot;
                for(int k = 0; k < indDes[c].n_elem; k++)
                {
                    input_XD_des.insert_cols(input_XD_des.n_cols,XD_des.col(indDes[c](k)));
                    input_XD_des_dot.insert_cols(input_XD_des_dot.n_cols, XD_des_dot.col(indDes[c](k)));
                }
                uD = controlDefender5(XD, SD, indD, assign[c],input_XD_des, input_XD_des_dot, indDes[c], uD, motionP_result.mP, times(ti), ND);
                uD.print("uD:");
                // if(ti == bound - 1) 
                // {
                //     motionP_result.dDf.XD_des0.print("motionP_result XD_des0:");
                //     XD.print("XD: ");
                //     XD_des.print("XD_des: ");
                // }
                for (int j = 0; j < indD.n_elem; j++)
                {
                    // cout << "defender j: " << j << endl;
                    // cout << "distance: " << arma::norm(XD.submat(0,indD(j),1,indD(j)) - XD_des.submat(0,indDes[c](j),1,indDes[c](j))) << endl;
                    if(arma::norm(XD.submat(0,indD(j),1,indD(j)) - XD_des.submat(0,indDes[c](j),1,indDes[c](j))) < 1e-3)
                    {
                        defReachCount(j) = 1;
                        defReachCount.print("reach count: ");
                        if(accu(defReachCount) >= ND)
                        {
                            flagSeek(c) = 1;
                            flagGather(c) = 0;
                            break;
                        }
                    }
                }
                // cout << "end gathering" << endl;
            } else if (flagGather(c) != 1 && flagSeek(c) == 1 && flagEnclose(c) != 1 && flagHerd(c) != 1)
            {
                cout<< "seeking!" << endl;
                uvec j1 = find(assign[c] == 1);
                uvec jND = find(assign[c] == ND);
                OpenForm OpenForm_result = defDesiredOpenForm(XDFc, c, &XD_des, &XD_des_dot, RDF_OPEN, rho_SN, indDes, NDinCluster(c), NClusterD(ti), YA, rAcm, vAcm, rhoA_con, NClusterAD(ti), c, phi(c), phi_dot(c), NA, ND, 0);
                XDFc[c] = modifiedDIDynamics(XDFc[c], uDFc_trans, dt, C_d);
                phi(c) += OpenForm_result.phi_dot * dt;
                if (phi(c) > 2*M_PI)
                {
                    phi(c) = phi(c) - 2*M_PI;
                }
                phi_dot(c) += dt * OpenForm_result.phi_ddot;
                mat input_XD_des;
                mat input_XD_des_dot;
                for (int i = 0; i < indDes[c].n_elem; i++)
                {
                    input_XD_des.insert_cols(input_XD_des.n_cols, XD_des.col(i));
                    input_XD_des_dot.insert_cols(input_XD_des_dot.n_cols, XD_des_dot.col(i));
                }
                
                uD = controlDefenderFormation4(XD, indD, input_XD_des, input_XD_des_dot, rhoA_con, YA, indAinClusterAD, NClusterAD(ti), 1, NA, ND);
                double ti_2 = ti;
                mat tempM = rAcm.col(c).t() * rDcm.col(c);
                if( ((arma::norm(rAcm.col(c) - rDcm.col(c)) < (2.9*rho_SN(c))) && flagAttInSight(c)) || ((arma::norm(rAcm.col(c) - rDcm.col(c)) < (3.5*rho_SN(c))) && tempM(0,0) > 0)) 
                {
                    flagEnclose(c) = 1;
                    flagSeek(c) = 0;
                    phi(c) = atan2(rD(1,assign[c].tail(1)(0)) - rD(1,assign[c](0)), rD(0,assign[c].tail(1)(0)) - rD(0,assign[c](0))) + M_PI/2;
                }
            } else if (flagGather(c) != 1 && flagSeek(c) != 1 && flagEnclose(c) == 1 && flagHerd(c) != 1)
            {
                cout <<"closing!" << endl;
                if (flagDForm != 1)
                {
                    flagDForm = 1;
                }

                for (int jj = 0; jj < indDes[c].n_elem; jj++)
                {
                    int j = indDes[c](jj);
                    double RD0 = RDF_CLOSED[c];
                    double thetaD;
                    if(flagSemiCircFormNotAchieved(c) == 0)
                    {
                        thetaD = phi(c)+M_PI/2+M_PI*(jj)/(NDinCluster(c)-1);
                    } else
                    {
                        thetaD = phi(c)+2*M_PI*(jj+1)/NDinCluster(c)-M_PI/NDinCluster(c);
                    }
                    mat tempM1(2,1);
                    mat tempM2(2,1);
                    tempM1(0,0) = cos(thetaD);
                    tempM1(1,0) = sin(thetaD);
                    tempM2 = rAcm.col(c) + RD0 * tempM1;
                    XD_des.col(j) = join_cols(tempM2, vAcm.col(c));
                    XD_des_dot.col(j) = zeros<mat>(4,1);

                }
                // Check if the terminal defenders should be connected or not
                int j1=assign[c](1); //Defender corresponding to first desired position
                int jND=assign[c](NDinCluster(c));   //Defender corresponding to last desired position

                int NAEnclosed = 0;
                for (int ii = 0; ii < NAinClusterAD(c); ii++)
                {
                    mat rA_input;
                    mat XD_input;
                    for (int Ai = 0; Ai < indAinClusterAD[c].n_elem; Ai++)
                    {
                        rA_input.insert_cols(rA_input.n_cols, rA.col(indAinClusterAD[c](Ai)));
                    }
                    for (int Di = 0; Di < indD.n_elem; Di++)
                    {
                        XD_input.insert_cols(XD_input.n_cols, XD.submat(0,indD(Di),1,indD(Di)));
                    }
                    if(poly_contain(XD_input, rA_input.col(ii)))
                    {
                        NAEnclosed = NAEnclosed+1;
                        flagAEnclosed(indAinClusterAD[c](ii)) = 1;
                    }
                    
                }

                if(NAEnclosed == NAinClusterAD(c))
                {
                    flagClusterAEnclosed(c) = 1;
                }

                if(((arma::norm(XD.submat(0, j1, 1, j1) - XD_des.submat(0, indDes[c](0), 1, indDes[c](0))) < bd 
                && arma::norm(XD.submat(0, jND, 1, jND) - XD_des.submat(0, indDes[c](NDinCluster(c)-1),1, indDes[c](NDinCluster(c)-1))) < bd ) 
                || arma::norm(XD.submat(0,j1, 1, j1) - XD.submat(0, jND, 1, jND)) <= 1.01*rho_SN(c)*sqrt(2-2*cos(2*M_PI/NDinCluster(c))))
                && NAEnclosed == NAinClusterAD(c)
                )
                {
                    WDString(j1, jND) = 1;
                    WDString(jND, j1) = 1;
                } else
                {
                    WDString(j1, jND) = 0;
                    WDString(jND, j1) = 0;
                }
                
                if((arma::norm(XD.submat(0, j1, 1, j1) - XD_des.submat(0, indDes[c](0), 1, indDes[c](0))) < bd 
                || arma::norm(XD.submat(0, jND, 1, jND) - XD_des.submat(0, indDes[c](NDinCluster(c))-1,1, indDes[c](NDinCluster(c)-1))) < bd )
                && flagSemiCircFormNotAchieved(c) == 0 )
                {
                    flagSemiCircFormNotAchieved(c) = 1;
                }

                int countDefConnect = 0;
                for(int j = 0; j < NDinCluster(c); j++)
                {
                    if(j >= NDinCluster(c))
                    {
                        if(WDString(indD(NDinCluster(c)-1), indD(0)) == 1)
                        {
                            countDefConnect++;
                        }
                    } else
                    {
                        if(WDString(indD(j), indD(j+1)) == 1)
                        {
                            countDefConnect++;
                        }
                    }
                }

                if(countDefConnect == NDinCluster(c))
                {
                    flagHerd(c) = 1;
                    flagEnclose(c) = 0;
                }
                mat XD_des_input;
                mat XD_des_dot_input;
                for(int input_i; input_i < indDes[c].n_elem; input_i++)
                {
                    XD_des_input.insert_cols(XD_des_input.n_cols, XD_des.col(indDes[c](input_i)));
                    XD_des_dot_input.insert_cols(XD_des_dot_input.n_cols, XD_des_dot.col(indDes[c](input_i)));
                }
                mat uD = controlDefenderFormation4(XD, indD, XD_des_input, XD_des_dot_input, rhoA_con, YA, indAinClusterAD, NClusterAD(ti), 1, NA, ND);

                ti_e(c) = ti;
                XDF_des = XD_des;
                XDFc[c] = join_cols(rAcm.col(c), vAcm.col(c));

            } else if (flagGather(c) != 1 && flagSeek(c) != 1 && flagEnclose(c) != 1 && flagHerd(c) == 1)
            {
                cout << "herding!" << endl;
                int j1 = assign[c](0);
                int jND = assign[c](NDinCluster(c)-1);

                if(arma::norm(XD.submat(0,j1,1,j1) - XD.submat(0,jND,1,jND)) <= R_DD_string)
                {
                    WDString(j1, jND) = 1;
                    WDString(jND, j1) = 1;
                } else
                {
                    WDString(j1, jND) = 0;
                    WDString(jND, j1) = 0;
                }

                double deltat = dt * (ti - ti_e(c));
                if(deltat < 0.02)
                {
                    vec distDS(NS);
                    for (int s = 0; s < NS; s++)
                    {
                        distDS(s) = arma::norm(XDFc[c].submat(0,0,1,0) - rS.col(s));
                    }
                    int min_index = distDS.index_min();
                    //938
                    defDesiredClosedForm(XDFc, c, rS.col(min_index), RDF_CLOSED, indDes, NDinCluster(c), NClusterD(ti), 
                                         phi(c), YA, rAcm, vAcm, rhoA_con, NClusterAD(ti), c, NA, ND, flagNDeven, 0, deltat, XD_des, XD_des_dot, uDFc_trans);
                    
                    if(ti - ti_e(c) < 0)
                    {
                        XDFc[c] = join_cols(XDFc[c].submat(0,0,1,0), zeros<vec>(2));
                    } else
                    {
                        XDFc[c] = modifiedDIDynamics(XDFc[c], uDFc_trans, dt, C_d);
                    }
                    mat XD_des_input;
                    mat XD_des_dot_input;
                    for(int input_i; input_i < indDes[c].n_elem; input_i++)
                    {
                        XD_des_input.insert_cols(XD_des_input.n_cols, XD_des.col(indDes[c](input_i)));
                        XD_des_dot_input.insert_cols(XD_des_dot_input.n_cols, XD_des_dot.col(indDes[c](input_i)));
                    }
                    mat uD = controlDefenderFormation4(XD, indD, XD_des_input, XD_des_dot_input, rhoA_con, YA, indAinClusterAD, NClusterAD(ti), 0 ,NA, ND);
                }

            }

        }
        cout << "2.0" << endl;
        XD_Des.col(ti+1) = XD_des.as_col();
        WDString_mat.slice(ti) = WDString;
        vA_dot_arr.col(ti) = vA_dot.as_col();
        cout << "2.1" << endl;
        for(int ii = 0; ii < NA; ii++)
        {
            if(arma::norm(rA.col(ii) - rP ) < rho_P)
            {
                uA.col(ii) = zeros<mat>(2,1);
            }
        }
        cout << "2.2" << endl;
        U.col(ti) = join_cols(uA.as_col(), uD.as_col());
        // double VA1_dot = 0;
        // for(int i = 0; i < NA; i++)
        // {
        //     mat tempM = vA.col(i).t() * uA.col(i);
        //     VA1_dot += tempM(0,0);
        // }
        cout << "2.3" << endl;
        X.col(ti+1) = modifiedDIDynamics(X.col(ti), U.col(ti), dt, C_d);
        XA = reshape(X.submat(0,ti+1, 4*NA-1, ti+1), 4, NA);
        YA = XA;
        rA = XA.submat(0,0,1, XA.n_cols-1);
        vA = XA.submat(2,0,3, XA.n_cols-1);
        cout << "2.4" << endl;
        XD = reshape(X.submat(4*NA, ti+1, 4*N-1, ti+1), 4, ND);
        mat XDp = reshape(X.submat(4*NA, ti, 4*N-1, ti), 4, ND);
        rD = XD.submat(0,0,1, XD.n_cols-1);
        mat vD = XD.submat(2,0,3, XD.n_cols-1);

        NClusterAD(ti+1) = NClusterAD(ti);
        cout << "2.5" << endl;
        // Saturate the velocity if beyond the maximum
        for(int i = 0; i < NA; i++)
        {
            double norm_vA = arma::norm(vA.col(i));
            if(norm_vA > v_maxA(i))
            {
                vA.col(i) = vA.col(i) * v_maxA(i) / norm_vA;
            }
        }
        cout << "2.5.1" << endl;
        for (int j = 0; j < ND; j++)
        {
            double norm_vD = arma::norm(vD.col(j));
            if(norm_vD > v_maxD(j))
            {
                vD.col(j) = vD.col(j) * v_maxD(j) / norm_vD;
            }
        }
        cout << "2.5.2" << endl;
        XA = join_cols(rA, vA);
        XD = join_cols(rD, vD);
        XA.print("XA: ");
        XD.print("XD: ");
        // cout << "2.5.3" << endl;
        // cout << "XD shape: " << XD.n_rows << " ," << XD.n_cols << endl;
        rDcm = arma::sum(XD.submat(0,0,1,ND-1), 1) / ND;
        vDcm = arma::sum(XD.submat(2,0,3,ND-1), 1) / ND;
        // cout << "2.6" << endl;
        for (int j = 0; j < ND; j++)
        {
            SD(j) = SD_arr(j, ti) + arma::norm(XD.submat(0,j,1,j) - XDp.submat(0,j,1,j));
        }
        SD_arr.col(ti+1) = SD;

        mat arr_minRDD(ND, 1);
        mat arr_minRAD(NA, ND);
        mat arr_minRAA(NA, 1);
        for(int j = 0; j < ND; j++)
        {
            if (j < ND-1)
            {
                for (int l = j+1; l < ND; l++)
                {
                    RDjDl(l) = arma::norm(rD.col(j) - rD.col(l));
                }
                arr_minRDD(j) = min(RDjDl.subvec(j+1, RDjDl.n_elem-1));
            }
            for (int i = 0; i < NA; i++)
            {
                arr_minRAD(i, j) = arma::norm(rD.col(j) - rA.col(i));
            }   
        }
        cout << "2.7" << endl;
        if(!arr_minRDD.is_empty())
        {
            minRDD(ti+1) = arr_minRDD.min();
        }

         for (int i = 0; i < NA; i++)
        {
            arr_minRAD(i,ND-1) = arma::norm(rD.col(ND-1)- rA.col(i));
        }
        minRAD(ti+1) = arr_minRAD.min();
        cout << "2.8" << endl;
        mat RAiAl(1,NA);
        for (int i = 0; i < NA; i++)
        {
            if (i < NA-1)
            {
                for (int ii = i+1; ii < NA; ii++)
                {
                    RAiAl.col(ii) = arma::norm(rA.col(i) - rA.col(ii));
                }
                arr_minRAA.resize(1,i+1);
                mat tempM;
                tempM = RAiAl.submat(0,i+1,0,RAiAl.n_cols-1);
                // cout <<"i: " << i << endl;
                arr_minRAA.col(i) = tempM.min();
            }
        }
        cout << "2.9" << endl; 
        if(!arr_minRAA.is_empty())
        {
            minRAA(ti+1) = arr_minRAA.min();
        }

        if (NA > 1)
        {
            minEAO(ti)=(R_m_AA-rhoA_safe)/(control_attacker_result.R_AO_min);
        }

        minRAAProjS(ti)=(R_m_AD-rhoAD_safe)/control_attacker_result.R_AAProjS_min;
        
        t=t+dt;
        times(ti+1)=t;

        out_ti = ti;
        cout << "3.0" << endl;

    }
    if(out_ti == Niter)
    {
        out_ti++;
    }
    cout << "deleting unnecessary elements" << endl;
    vA_dot_arr.col(out_ti) = vA_dot_arr.col(out_ti-1);
    FA_dot_arr.col(out_ti) = FA_dot_arr.col(out_ti-1);
    U.col(out_ti) = U.col(out_ti-1);
    RefTraj.col(out_ti) = RefTraj.col(out_ti-1);
    minEDO(out_ti) = minEDO(out_ti-1);
    minEAO(out_ti) = minEAO(out_ti-1);
    minRAAProjS(out_ti)=minRAAProjS(out_ti-1);
    times = times.submat(0,0,times.n_rows-1,out_ti + 1);
    X = X.submat(0,0,X.n_rows-1, out_ti+1);
    // X.col(out_ti).print("X col out_ti: ");
    // X.col(out_ti+1).print("col out_ti+1: ");
    U = U.submat(0,0,U.n_rows-1,out_ti+1);
    XD_Des = XD_Des.submat(0,0,XD_Des.n_rows-1, out_ti+1);
    vA_dot_arr = vA_dot_arr.submat(0,0,vA_dot_arr.n_rows-1, out_ti+1);
    FA_dot_arr = FA_dot_arr.submat(0,0,FA_dot_arr.n_rows-1, out_ti+1);
    minRDD = minRDD.submat(0,0,minRDD.n_rows-1, out_ti+1);
    minRAD = minRAD.submat(0,0,minRAD.n_rows-1, out_ti+1);
    minRAA = minRAA.submat(0,0,minRAA.n_rows-1, out_ti+1);
    minEDO = minEDO.submat(0,0,minEDO.n_rows-1, out_ti+1);
    minEAO = minEAO.submat(0,0,minEAO.n_rows-1, out_ti+1);
}






int main() {
    int NA = 18;
    int ND = NA;
    calAllparametersExperiment(NA,ND);
    AllocateMemory(NA,ND);
    measurements(NA,ND);
    initial_contorl(NA, ND);
    getMotionPlan(NA,ND);
    // motionP_result.dDf.XD_des0.print("XD_des0");
    calDistance(NA,ND);
    control_loop(NA,ND);

    // XD.load("/home/damon/Downloads/multi_swarm/controlD/XD.txt");
    // SD.load("/home/damon/Downloads/multi_swarm/controlD/SD.txt");
    // XD.print("XD: ");
    // SD.print("XD: ");
    // mat indD;
    // mat XD_des, XD_des_dot;
    // indD.load("/home/damon/Downloads/multi_swarm/controlD/indD.txt");
    // indD = indD-1;
    // XD_des.load("/home/damon/Downloads/multi_swarm/controlD/XD_des.txt");
    // XD_des_dot.load("/home/damon/Downloads/multi_swarm/controlD/XD_des_dot.txt");
    // mat indDes;
    // indDes.load("/home/damon/Downloads/multi_swarm/controlD/indDes.txt");
    // indDes = indDes -1;
    // mat uD;
    // uD.load("/home/damon/Downloads/multi_swarm/controlD/uD.txt");
    // std::vector<vec> assign;
    // assign.push_back(assignment);

    // uD = controlDefender5(XD, SD, indD, assignment, XD_des, XD_des_dot, indDes, uD, motionP_result.mP, 0, ND);
    // uD.print("uD after control Defender");
}
//     // motionP_result.mP.assign.print("assign");
//     // control_loop(NA,ND);
//     // XD.load("/home/damon/Downloads/multi_swarm/controlD/XD.txt");
//     // // XD.print("XD: ");
//     // SD.load("/home/damon/Downloads/multi_swarm/controlD/SD.txt");
//     // mat indD;
//     // mat XD_des, XD_des_dot;
//     // indD.load("/home/damon/Downloads/multi_swarm/controlD/indD.txt");
//     // indD = indD-1;
//     // XD_des.load("/home/damon/Downloads/multi_swarm/controlD/XD_des.txt");
//     // XD_des_dot.load("/home/damon/Downloads/multi_swarm/controlD/XD_des_dot.txt");
//     // mat indDes;
//     // indDes.load("/home/damon/Downloads/multi_swarm/controlD/indDes.txt");
//     // indDes = indDes -1;
//     // mat uD;
//     // uD.load("/home/damon/Downloads/multi_swarm/controlD/uD.txt");
//     // // std::vector<vec> assign;
//     // // assign.push_back(assignment);
    
//     // uD = controlDefender5(XD, SD, indD, assignment, XD_des, XD_des_dot, indDes, uD, motionP_result.mP, 0, ND);
//     // uD.print("uD after control Defender");

//     // mat tempM(4,1,fill::zeros);
//     // tempM(0) = 385.153071102075;
//     // tempM(1) = -176.086389266759;
//     // XDFc.reset();
//     // XDFc.set_size(3,1);
//     // XDFc(0) = tempM;
//     // tempM(0) = 443.375146160435;
//     // tempM(1) = 37.3666215895702;
//     // XDFc(1) = tempM;
//     // tempM(0) = 512.183053047588;
//     // tempM(1) = 289.629270783414;
//     // XDFc(2) = tempM;
//     // XDFc.print("XDFc: ");
//     // int clusterDNum = 1;
//     // mat XD_des, XD_des_dot, RDF_OPEN;
//     // XD_des.load();
//     // XD_des_dot.load();
//     // RDF_OPEN.load();
//     // rho_SN.load();
//     // std::vector<vec> indDes;
//     // vec tempV = {14,15,16,17,18};
//     // tempV = tempV - 1;
//     // indDes.push_back(tempV);
//     // tempV = {8,9,10,11,12,13};
//     // tempV = tempV - 1;
//     // indDes.push_back(tempV);
//     // tempV = {1,2,3,4,5,6,7};
//     // tempV = tempV - 1;
//     // indDes.push_back(tempV);
//     // int NDinCluster = 5;
//     // int NClusterD = 3;
//     // XA.load();
//     // rAcm.load();
//     // vAcm.load();
//     // rhoA_con.load();
//     // int NClusterAD = 3;

    





//     // motionP_result.mP.assign.print("assign: ");
// //     AllocateMemory();
// //     measurements(2);
// //     initial_contorl();
// //     // arma::vec rAcm = {5.9807,-33.2064};
// //     // arma::vec rP = {14,-1};
// //     // path_elem p = findShortestPath(rAcm, rP);
// //     // p.rV.print("rV: ");
// //     // p.S.print("S: ");
// //     // double q = 16.0949;
// //     // CoorOnPath A = findCoordOnPath(q, p);
// //     // A.rp.print("rDFc0: ");
// //     // A.thetap.print("thetaAcom0: ");
// //     // v_maxD.print("v_maxD: ");
// //     // v_maxDC.print("v_maxDC: ");
// //     // u_maxD.print("u_maxD: ");
// //     // cout << "rho_P:" << rho_P << endl;
// //     // cout << "rho_safe:" << rho_safe << endl;
// //     // cout << "rho_sn:" << rho_sn << endl;
// //     // cout << "rho_Acon:" << rho_Acon << endl;
// //     // YA.print("YA:");
    
// // //     XD0.print("XD0: ");
// //     // cout << "RDF_open:" << RDF_open << endl;
// //     // YA.load("../../../../../Downloads/swarm_matlab/controlD/YA.txt");
// //     // XD0.load("../../../../../Downloads/swarm_matlab/controlD/XD0.txt");
// //     getMoitonPlan();
// //     calDistance();
// //     checkFormation();
//     // mat XD_des = motionP_result.dDf.XD_des0;
//     // mat XD_des_dot = motionP_result.dDf.XD_des_dot0;
//     // std::ifstream fin("../../../../../Downloads/swarm_matlab/controlD/t.txt");
//     // double time;
//     // fin >> time;
//     // // time = 10;
//     // std::cout << "t: "<<time << std::endl;

//     // fin.close();
//     // XD.load("../../../../../Downloads/swarm_matlab/controlD/XD.txt");
//     // XD.print("XD: ");
//     // SD.load("../../../../../Downloads/swarm_matlab/controlD/SD.txt");
//     // SD.print("SD: ");
//     // // indDef.print("indDef: ");
//     // motionP_result.mP.assign.print("assign: ");
//     // XD_des.print("XD_des: ");
//     // XD_des_dot.print("XD_des_dot: ");
//     // mat uD = controlDefender5(XD,SD, regspace(1,1,4), motionP_result.mP.assign, XD_des, XD_des_dot, motionP_result.mP, time,3);
//     // uD.print("uD: ");
    
//     // mat XDFc;
//     // XDFc.load("../../../../../Downloads/swarm_matlab/OpenForm/XDFc.txt");
//     // XA.load("../../../../../Downloads/swarm_matlab/OpenForm/XA.txt");
//     // XDFc.print("XDFc: ");
//     // XA.print("XA: ");
//     // std::ifstream fin("../../../../../Downloads/swarm_matlab/OpenForm/RDF0.txt");
//     // double RDF0;
//     // fin >> RDF0;
//     // std::cout << "RDF0: "<< RDF0 << std::endl;
//     // fin.close();

//     // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/phi.txt");
//     // double phi;
//     // fin >> phi;
//     // std::cout << "phi: "<< phi << std::endl;
//     // fin.close();

//     // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/phi_dot.txt");
//     // double phi_dot;
//     // fin >> phi_dot;
//     // std::cout << "phi_dot: "<< phi_dot << std::endl;
//     // fin.close();
//     // mat XD_des, XD_des_dot, uDFf_trans;
//     // double phi_ddot;
//     // defDesiredOpenForm(XDFc, RDF0, XA, phi, phi_dot, NA, ND, &XD_des, &XD_des_dot, &phi_ddot, &uDFf_trans, &flagAttInSight);
//     // XD_des.print("using pointer, XD_des:");
//     // XD_des_dot.print("using pointer, XD_des_dot: ");
//     // cout << "using pointer, phi_ddot: "<<phi_ddot << endl;
    
//     // XD.load("../../../../../Downloads/swarm_matlab/controlDF/XD.txt");
//     // indDef.load("../../../../../Downloads/swarm_matlab/controlDF/indDef.txt");
//     // mat XD_des,XD_des_dot, uDFc_trans;
//     // XD_des.load("../../../../../Downloads/swarm_matlab/controlDF/XD_des.txt");
//     // XD_des_dot.load("../../../../../Downloads/swarm_matlab/controlDF/XD_des_dot.txt");
//     // uDFc_trans.load("../../../../../Downloads/swarm_matlab/controlDF/uDFc_trans.txt");
//     // XA.load("../../../../../Downloads/swarm_matlab/controlDF/XA.txt");
//     // XD.print("XD: ");
//     // indDef.print("indDef: ");
//     // XD_des.print("XD_des: ");
//     // XD_des_dot.print("XD_des_dot: ");
//     // uDFc_trans.print("uDFc_trans: ");
    
//     // mat Abc = controlDefenderFormation4(XD, indDef, motionP_result.mP.assign, XD_des, XD_des_dot, uDFc_trans, XA, NA, 3, 1);
//     // // mat Abc = controlFiniteTimeTrajTracking(XD,indDef, XD_des, XD_des_dot,uDFc_trans, XA, ND, 1);
//     // Abc.print("Abc:");
//     // uD.print("uD: ");

//     // mat XDFc;
//     // XDFc.load("../../../../../Downloads/swarm_matlab/OpenForm/XDFc.txt");
//     // XA.load("../../../../../Downloads/swarm_matlab/OpenForm/XA.txt");

//     // std::ifstream fin("../../../../../Downloads/swarm_matlab/OpenForm/RDF0.txt");
//     // double RDF0;
//     // fin >> RDF0;
//     // fin.close();

//     // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/phi0.txt");
//     // double phi0;
//     // fin >> phi0;
//     // fin.close();
    
//     // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/delta_t.txt");
//     // double delta_t;
//     // fin >> delta_t;
//     // fin.close();

//     // XDFc.print("XDFc:");
//     // cout << "RDF0: " << RDF0 << endl;
//     // cout << "phi0: " << phi0 << endl;
//     // XA.print("XA: ");
//     // cout << "delta_t: " << delta_t << endl;

//     // mat XD_des, XD_des_dot, uDFc_trans;
//     // XD_des.print("XD_des");
//     // XD_des_dot.print("XD_des_dot:");
//     // uDFc_trans.print("uDFc_trans: ");

//     // defDesiredClosedForm(XDFc, RDF0, phi0, XA,NA,ND, 1, delta_t, &XD_des, &XD_des_dot, &uDFc_trans);
//     // XD_des.print("XD_des pointer");
//     // XD_des_dot.print("XD_des_dot pointer:");
//     // uDFc_trans.print("uDFc_trans pointer: ");


//     // mat X0;
//     // X0.load("../../../../../Downloads/swarm_matlab/OpenForm/X0.txt");
//     // U.load("../../../../../Downloads/swarm_matlab/OpenForm/U.txt");
//     // std::ifstream fin("../../../../../Downloads/swarm_matlab/OpenForm/dt.txt");
//     // double dt;
//     // fin >> dt;
//     // fin.close();

//     // X0.print("X0:");
//     // U.print("U");
//     // cout << "dt: " << dt << endl;
//     // cout << "C_d: " << C_d << endl;
//     // mat X1 = modifiedDIDynamics(X0,U, dt, C_d);
//     // X1.print("X1: ");
    
// }