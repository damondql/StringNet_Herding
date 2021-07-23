#pragma once
#include <armadillo>
#include <istream>
#include "findShortestPath.cpp"
#include "findPathSpeeds.cpp"
#include "pathIntersections.cpp"
#include "defGoalAssignMIQP.cpp"
using namespace std;
using namespace arma;

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

// struct interSec_elem
// {
//     int Flag;
//     int Pbar1;
//     int Pbar2;
//     int flag0;
// };


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

void motionPlanForDefOpenForm(mat XD, mat XD_des, int ND, int flagPlotPaths, int figNumber) {
    // cout << "00000000000000000000" << endl;
    // motionPlan motionPlan_resutl;
    double umd = 0.8 * u_maxD(0);
    double vmd = sqrt(umd/C_d);
    std::vector<std::vector<path_elem>> Path(ND, std::vector<path_elem>(ND));
    std::vector<std::vector<pathVel_elem>> pathVel(ND, std::vector<pathVel_elem>(ND));
    mat optT(ND, ND, fill::zeros);
    // cout << "11111111111111111" << endl;
    for (int jj = 0; jj < ND; jj++)
    {
        for (int j = 0; j < ND; j++)
        {
            Path[jj][j] = findShortestPath(XD.submat(0,j,1,j), XD_des.submat(0,jj,1,jj));
            // cout << "1.1" << endl;
            pathVel[jj][j] = findPathSpeeds(Path[jj][j], v_maxDC(j), umd, vmd);
            // cout << "1.2" << endl;
            optT(jj, j) = pathVel[jj][j].T(1); // obs free, T  is only size 2
            // cout << "1.3" << endl;
            // cout << "row: " << jj << " col: " << j <<endl; 
            // Path[jj][j].rV.print("rV:");
            // pathVel[jj][j].T.print("T");
        }
    }

    // cout << "22222222222222222" << endl;
    std::vector<path_elem> Path0;
    for (int jj = 0; jj < ND; jj++) {
        for (int j = 0; j < ND; j++)
        {
            Path0.push_back(Path[j][jj]);
        }
    }

    // for (int i = 0; i < ND*ND; i++)
    // {
    //     cout << i << endl;
    //     Path0[i].rV.print("rV: ");
    // }
    
    mat Pbar;
    Pbar = zeros<mat>(pow(ND,2), pow(ND,2));
    std::vector<std::vector<interSec_elem>> interSec(pow(ND,2), std::vector<interSec_elem>(pow(ND,2)));

    vec optT_vec(ND*ND, fill::zeros);
    for (int i = 0; i < ND; i++)
    {
        for (int j = 0; j < ND; j++)
        {
            optT_vec(i+j) = optT(j,i);
        }
        
    }
    // cout << "3333333333333333333" << endl;
    int count_time=0;
    for (int j = 0; j < pow(ND,2); j++)
    {
        for (int jj = j + 1; jj < pow(ND,2); jj++)
        {   
            count_time++;
            // cout << "3.1" << endl;
            cout << "times: " << count_time << endl;
            interSec[j][jj] = pathIntersections(Path0[j], Path0[jj]);
            // cout << "3.2" << endl;
            Pbar(j,jj) = interSec[j][jj].Pbar1;
            // cout << "3.3" << endl;
            Pbar(jj,j) = interSec[j][jj].Pbar2;
            // cout << "3.4" << endl;
        }
        
    }
    Pbar.print("Pbar: ");
    // cout << "4444444444444444444" << endl;
    goal_assign goal = defGoalAssignMIQP(optT_vec, Pbar, ND);
    // cout << "5555555555555555555" << endl;
    goal.assign.print("assign: ");
    cout<<"cost: " <<goal.cost << endl;
    cout <<"maxOptT: " << goal.maxOptT<<endl;
    std::vector<path_elem> assignedPath(ND);
    std::vector<pathVel_elem> assignedPathVel(ND);
    vec assignedOptT(ND);
    for (int j = 0; j < ND; j++)
    {
        assignedPath[j] = Path[goal.assign(j)-1][j];
        assignedPathVel[j] = pathVel[goal.assign(j)-1][j];
        assignedOptT[j] = pathVel[goal.assign(j)-1][j].T.tail(1)(0);
    }

    
    // return motionPlan_resutl;
}

int main(){
    // vec XD = {11.9907,-10.7580};
    // vec XD_des = {14.808227436999694,-18.0334291835457};
    // path_elem P =  findShortestPath(XD, XD_des);
    // pathVel_elem pV = findPathSpeeds(P,0.3924,0.0786,0.4432);
    // pV.T.print("T:");
    // pV.v.print("v:");
    // cout << "s_bar1: "<<pV.s_bar1<<endl;
    // cout << "s_bar2: "<<pV.s_bar2<<endl;
    // cout << "T_bar1: "<<pV.T_bar1<<endl;
    // cout << "T_bar2: "<<pV.T_bar2<<endl;
    // cout << "v_bar: "<<pV.v_bar<<endl;
    // cout << "v_maxD: "<<pV.v_maxD<<endl;
    // cout << "v_maxDC: "<<pV.v_maxDC<<endl;
    mat XD;
    XD.load("../../../../../Downloads/swarm_matlab/motionPlan/XD.txt");
    mat XD_des;
    XD_des.load("../../../../../Downloads/swarm_matlab/motionPlan/XD_des.txt");
    XD.print("XD: ");
    XD_des.print("XD_des: ");
    u_maxD = {0.0982320000000000, 0.0982320000000000, 0.0982320000000000};
    v_maxD = {0.495560288965934,0.49556028896593,0.495560288965934};
    v_maxDC = {0.392380551560576,0.392380551560576,0.392380551560576};
    motionPlanForDefOpenForm(XD, XD_des, 3, 0, 1);

}