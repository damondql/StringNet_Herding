#include "AllParameters.cpp"
// #include "symDerivative.cpp"
#include "findCommGraphAndFormDist.cpp"
#include "defInitDesirePos.cpp"
// #include "controlAttacker4.cpp"
// #include "controlDefender5.cpp"
// #include "defDesiredOpenForm.cpp"
// #include "controlFiniteTimeTrajTracking.cpp"
// #include "inhull.cpp"
// #include "controlDefenderFormation4.cpp"
// #include "defDesiredClosedForm.cpp"
// #include "modifiedDIDynamics.cpp"
// #include "findCoordOnPath.cpp"
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
    double sf_RDF_open=1.05;
    double sf_RDF_closed=1.1;
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
field<mat> XDFc;
vec assignment;
void getMotionPlan(int NA, int ND){
    XDFc.set_size(ND);
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
    XDFc(0) = join_cols(motionP_result.dDf.rDFc0, zeros<mat>(2,1));
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
void calDistance(int NA, int ND){
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
    vec RDjDl(ND);
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
int flagAttInSight=0;
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

int NDinCluster;
vec NClusterD;
vec assignSafeArea;

int flagSemiCircFormNotAchieved;
vec flagGather;
vec flagSeek;
vec flagEnclose;
vec flagHerd;

int flagAttackerStayTogether=1;
int splitCout = 0;

field<vec> ClusterRad;

void control_loop(int NA, int ND){
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
    NDinCluster = ND;
    NClusterD.resize(1);
    NClusterD(0) = 1;
    assignSafeArea = zeros<vec>(1);
    flagSemiCircFormNotAchieved = 0;
    flagGather.resize(1);
    flagGather(0) = 1;
    flagSeek.resize(1);
    flagSeek(0) = 0;
    flagEnclose = zeros<vec>(1);
    flagHerd = zeros<vec>(1);
    flagAttackerStayTogether = 1;
    splitCout = 0;
    ClusterRad.set_size(Niter,1);

    mat RefTraj(4*ND, Niter);
    mat SigmaProdD_arr;
    mat rAcm_arr;
    mat vAcm_arr;
    mat sigmaProd;
    mat minRAAProjS(Niter,1);
    int out_ti;
    // int bound = Niter;
    mat uD;
    mat XD_des;
    mat XD_des_dot;
    mat uDFc_trans;
    mat XDF_des;
    cube WDString_mat(ND,ND, Niter);
    vec RDF_OPEN;
    vec RDF_CLOSED;
    int bound = 1;
    int countAinS = 0;
    int countAinP = 0;
    std::vector<mat> rAcm_new;
    std::vector<vec> indexOfNewClusterAD;
    cout<< "enter control loop" << endl;
    for (int ti = 0; ti < bound; ti++)
    {
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
            indAinClusterAD[c] = find(clusterIdAD == c);
            indAinClusterAD[c].print("indAinClusterAD:");
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
            clusterRad.print("clusterRad: ");
            
        }
        rAcm.print("rAcm:");
        ClusterRad(ti) = clusterRad;

        int flagSpliAttackers = 0;
        int NClusterAD0 = NClusterAD(ti);

        for (int c = 0; c < NClusterAD(ti); c++)
        {
            if(c == 0){
                rAcm_new.clear(); //reset std::vector
            }
            if(!flagClusterAEnclosed(c))
            {
                flagSpliAttackers = 1;
                if(flagSpliAttackers){
                    vec clusterIdAD0 = DBSCAN(XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)), epsilon, MinPts);
                    int NNewClusterAD = clusterIdAD0.max();
                    if(NNewClusterAD > 0)
                    {
                        flagToSplitDefender(c) = 1;
                        splitCout++;
                    }
                    // indAinClusterAD.clear();
                    std::vector<uvec> indAinClusterAD0;
                    mat tempM(2, NNewClusterAD);
                    for (int cc = 0; cc < NNewClusterAD; cc++)
                    {
                        indAinClusterAD0.push_back( indAinClusterAD[c](arma::find(clusterIdAD0 == cc)));
                        tempM.col(cc) = arma::sum(XA.submat(0,indAinClusterAD0[c](0),1,indAinClusterAD0[c](indAinClusterAD0[c].n_elem-1)),1)/indAinClusterAD0[cc].n_elem;
                    }
                    rAcm_new.push_back(tempM);
                    
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
                    vec tempV(NNewClusterAD,fill::zeros);
                    assignSafeArea.resize(NNewClusterAD);
                    for (int cc = 0; cc < NNewClusterAD-1; cc++)
                    {
                        int c1 = NClusterAD0 + cc;
                        indAinClusterAD[c1] = indAinClusterAD0[cc+1];
                        for (int i = 0; i < indAinClusterAD[c1].n_elem; i++)
                        {
                            clusterIdAD(indAinClusterAD[c1](i)) = c1;
                        }
                        tempV(cc+1) = c1;
                        assignSafeArea(c1) = 0; //Initilization of the assigned safe area
                    }
                    indexOfNewClusterAD.push_back(tempV);

                    indAinClusterAD[c] = indAinClusterAD0[0];
                    for (int i = 0; i < indAinClusterAD[c].n_elem; i++)
                    {
                        clusterIdAD(indAinClusterAD[c](i)) = c;
                    }
                    indexOfNewClusterAD[c](0) = c;
                    NClusterAD0 = NClusterAD0+NNewClusterAD-1;

                }    
            }
        }

        NClusterAD(ti) = NClusterAD0; //Number of clusters at time ti after clustering
        arma::vec NAinClusterAD(NClusterAD(ti));
        if(flagSpliAttackers == 1)
        {
            RDF_OPEN.resize(NClusterAD(ti));
            RDF_CLOSED.resize(NClusterAD(ti));
            for (int c = 0; c < NClusterAD(ti); c++)
            {
                clusterRad(c) = 0;
                indAinClusterAD[c] = arma::find(clusterIdAD == c);
                rAcm.resize(2,c+1);
                vAcm.resize(2,c+1);
                rAcm.col(c) = arma::sum(XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/indAinClusterAD[c].n_elem;
                vAcm.col(c) = arma::sum(XA.submat(2,indAinClusterAD[c](0),3,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/indAinClusterAD[c].n_elem;
                NAinClusterAD(c) = indAinClusterAD[c].n_elem;
                for (int ii = 0; ii < NAinClusterAD.n_elem; ii++)
                {
                    int i = indAinClusterAD[c](ii);
                    double RAcA =arma::norm(rAcm.col(c) - XA.submat(0,i,1,i));
                    if(RAcA > clusterRad(c))
                    {
                        clusterRad(c) = RAcA;
                    }
                }
                rhoA_con.resize(c+1);
                rho_SN.resize(c+1);
                RDF_OPEN.resize(c+1);
                RDF_CLOSED.resize(c+1);
                rhoA_con(c) = rhoA_con_fun(NAinClusterAD(c), NA, ND, R_DD_string, rhoD_safe, R_m_DD);
                rho_SN(c) = rhoA_con(c) + rhoD_safe+5;
                RDF_OPEN(c) = 1.2*M_PI/2*rho_SN(c);
                RDF_CLOSED(c) = 1.1*rho_SN(c);
                flagClusterAEnclosed(c) = 0;
            }
        }

        for (int c = 0; c < NClusterAD(ti); c++)
        {
            NAinClusterAD(c) = indAinClusterAD[c].n_elem;
            rAcm.resize(2,c+1);
            vAcm.resize(2,c+1);
            rAcm.col(c) = arma::sum(XA.submat(0,indAinClusterAD[c](0),1,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/ NAinClusterAD(c);
            vAcm.col(c) = arma::sum(XA.submat(2,indAinClusterAD[c](0),3,indAinClusterAD[c](indAinClusterAD[c].n_elem-1)),1)/ NAinClusterAD(c);
        }

        ClusterIdAD.col(ti) = clusterIdAD.as_col();
        std::vector<vec> indDes;
        std::vector<vec> assign;
        if(flagGather(0) == 1)
        {
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
                                indexOfNewClusterAD[0](temp1(i) = indexOfNewClusterAD[cc](temp2(i));
                            }
                        }
                    }
                    
                }
                
            }
        } else 
        {
            NClusterD(ti) = NClusterD(ti-1);
            int NClusterD0 = NClusterD(ti);
            for (int c = 0; i < NclusterD(ti); i++)
            {
                if(flagToSplitDefender(c))
                {

                }
            }
            
        }
        


    }
    
}






int main() {
    int NA = 18;
    int ND = NA;
    calAllparametersExperiment(NA,ND);
    AllocateMemory(NA,ND);
    measurements(NA,ND);
    initial_contorl(NA, ND);
    getMotionPlan(NA,ND);
    calDistance(NA,ND);
    // motionP_result.mP.assign.print("assign");
    // control_loop(NA,ND);



//     AllocateMemory();
//     measurements(2);
//     initial_contorl();
//     // arma::vec rAcm = {5.9807,-33.2064};
//     // arma::vec rP = {14,-1};
//     // path_elem p = findShortestPath(rAcm, rP);
//     // p.rV.print("rV: ");
//     // p.S.print("S: ");
//     // double q = 16.0949;
//     // CoorOnPath A = findCoordOnPath(q, p);
//     // A.rp.print("rDFc0: ");
//     // A.thetap.print("thetaAcom0: ");
//     // v_maxD.print("v_maxD: ");
//     // v_maxDC.print("v_maxDC: ");
//     // u_maxD.print("u_maxD: ");
//     // cout << "rho_P:" << rho_P << endl;
//     // cout << "rho_safe:" << rho_safe << endl;
//     // cout << "rho_sn:" << rho_sn << endl;
//     // cout << "rho_Acon:" << rho_Acon << endl;
//     // YA.print("YA:");
    
// //     XD0.print("XD0: ");
//     // cout << "RDF_open:" << RDF_open << endl;
//     // YA.load("../../../../../Downloads/swarm_matlab/controlD/YA.txt");
//     // XD0.load("../../../../../Downloads/swarm_matlab/controlD/XD0.txt");
//     getMoitonPlan();
//     calDistance();
//     checkFormation();
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
    
    // mat XDFc;
    // XDFc.load("../../../../../Downloads/swarm_matlab/OpenForm/XDFc.txt");
    // XA.load("../../../../../Downloads/swarm_matlab/OpenForm/XA.txt");
    // XDFc.print("XDFc: ");
    // XA.print("XA: ");
    // std::ifstream fin("../../../../../Downloads/swarm_matlab/OpenForm/RDF0.txt");
    // double RDF0;
    // fin >> RDF0;
    // std::cout << "RDF0: "<< RDF0 << std::endl;
    // fin.close();

    // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/phi.txt");
    // double phi;
    // fin >> phi;
    // std::cout << "phi: "<< phi << std::endl;
    // fin.close();

    // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/phi_dot.txt");
    // double phi_dot;
    // fin >> phi_dot;
    // std::cout << "phi_dot: "<< phi_dot << std::endl;
    // fin.close();
    // mat XD_des, XD_des_dot, uDFf_trans;
    // double phi_ddot;
    // defDesiredOpenForm(XDFc, RDF0, XA, phi, phi_dot, NA, ND, &XD_des, &XD_des_dot, &phi_ddot, &uDFf_trans, &flagAttInSight);
    // XD_des.print("using pointer, XD_des:");
    // XD_des_dot.print("using pointer, XD_des_dot: ");
    // cout << "using pointer, phi_ddot: "<<phi_ddot << endl;
    
    // XD.load("../../../../../Downloads/swarm_matlab/controlDF/XD.txt");
    // indDef.load("../../../../../Downloads/swarm_matlab/controlDF/indDef.txt");
    // mat XD_des,XD_des_dot, uDFc_trans;
    // XD_des.load("../../../../../Downloads/swarm_matlab/controlDF/XD_des.txt");
    // XD_des_dot.load("../../../../../Downloads/swarm_matlab/controlDF/XD_des_dot.txt");
    // uDFc_trans.load("../../../../../Downloads/swarm_matlab/controlDF/uDFc_trans.txt");
    // XA.load("../../../../../Downloads/swarm_matlab/controlDF/XA.txt");
    // XD.print("XD: ");
    // indDef.print("indDef: ");
    // XD_des.print("XD_des: ");
    // XD_des_dot.print("XD_des_dot: ");
    // uDFc_trans.print("uDFc_trans: ");
    
    // mat Abc = controlDefenderFormation4(XD, indDef, motionP_result.mP.assign, XD_des, XD_des_dot, uDFc_trans, XA, NA, 3, 1);
    // // mat Abc = controlFiniteTimeTrajTracking(XD,indDef, XD_des, XD_des_dot,uDFc_trans, XA, ND, 1);
    // Abc.print("Abc:");
    // uD.print("uD: ");

    // mat XDFc;
    // XDFc.load("../../../../../Downloads/swarm_matlab/OpenForm/XDFc.txt");
    // XA.load("../../../../../Downloads/swarm_matlab/OpenForm/XA.txt");

    // std::ifstream fin("../../../../../Downloads/swarm_matlab/OpenForm/RDF0.txt");
    // double RDF0;
    // fin >> RDF0;
    // fin.close();

    // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/phi0.txt");
    // double phi0;
    // fin >> phi0;
    // fin.close();
    
    // fin.open("../../../../../Downloads/swarm_matlab/OpenForm/delta_t.txt");
    // double delta_t;
    // fin >> delta_t;
    // fin.close();

    // XDFc.print("XDFc:");
    // cout << "RDF0: " << RDF0 << endl;
    // cout << "phi0: " << phi0 << endl;
    // XA.print("XA: ");
    // cout << "delta_t: " << delta_t << endl;

    // mat XD_des, XD_des_dot, uDFc_trans;
    // XD_des.print("XD_des");
    // XD_des_dot.print("XD_des_dot:");
    // uDFc_trans.print("uDFc_trans: ");

    // defDesiredClosedForm(XDFc, RDF0, phi0, XA,NA,ND, 1, delta_t, &XD_des, &XD_des_dot, &uDFc_trans);
    // XD_des.print("XD_des pointer");
    // XD_des_dot.print("XD_des_dot pointer:");
    // uDFc_trans.print("uDFc_trans pointer: ");


    // mat X0;
    // X0.load("../../../../../Downloads/swarm_matlab/OpenForm/X0.txt");
    // U.load("../../../../../Downloads/swarm_matlab/OpenForm/U.txt");
    // std::ifstream fin("../../../../../Downloads/swarm_matlab/OpenForm/dt.txt");
    // double dt;
    // fin >> dt;
    // fin.close();

    // X0.print("X0:");
    // U.print("U");
    // cout << "dt: " << dt << endl;
    // cout << "C_d: " << C_d << endl;
    // mat X1 = modifiedDIDynamics(X0,U, dt, C_d);
    // X1.print("X1: ");
    
}