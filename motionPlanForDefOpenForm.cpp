#pragma once
#include <armadillo>
#include <istream>
#include "findShortestPath.cpp"
#include "findPathSpeeds.cpp"
#include "pathIntersections.cpp"
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

motionPlan motionPlanForDefOpenForm(mat XD, mat XD_des, int ND, int flagPlotPaths, int figNumber) {
    double umd = 0.8 * u_maxD(0);
    double vmd = sqrt(umd/C_d);
    std::vector<std::vector<path_elem>> Path(ND, std::vector<path_elem>(ND));
    std::vector<std::vector<pathVel_elem>> pathVel(ND, std::vector<pathVel_elem>(ND));
    mat optT(ND, ND, fill::zeros);
    for (int jj = 0; jj < ND; jj++)
    {
        for (int j = 0; j < ND; j++)
        {
            Path[jj][j] = findShortestPath(XD.submat(0,j,1,j), XD_des.submat(0,jj,1,jj));
            pathVel[jj][j] = findPathSpeeds(Path[jj][j], v_maxDC(j), umd, vmd);
            optT(jj, j) = pathVel[jj][j].T(1); // obs free, T  is only size 2
        }
    }
    std::vector<path_elem> Path0;
    for (int jj = 0; jj < ND; jj++) {
        for (int j = 0; j < ND; j++)
        {
            Path0.push_back(Path[j][jj]);
        }
    }
    mat Pbar;
    Pbar = zeros<mat>(pow(ND,2), pow(ND,2));
    std::vector<std::vector<interSec_elem>> interSec(pow(ND,2), std::vector<interSec_elem>(pow(ND,2)));
    for (int j = 0; j < pow(ND,2); j++)
    {
        for (int jj = j + 1; jj < pow(ND,2); jj++)
        {
            interSec[j][jj] = pathIntersections(Path0[j], Path0[jj]);
            Pbar(j,jj) = interSec[j][jj].Pbar1;
            Pbar(jj,j) = interSec[j][jj].Pbar2;
        }
        
    }
    

}

int main(){
    vec XD = {11.9907,-10.7580};
    vec XD_des = {14.808227436999694,-18.0334291835457};
    path_elem P =  findShortestPath(XD, XD_des);
    pathVel_elem pV = findPathSpeeds(P,0.3924,0.0786,0.4432);
    pV.T.print("T:");
    pV.v.print("v:");
    cout << "s_bar1: "<<pV.s_bar1<<endl;
    cout << "s_bar2: "<<pV.s_bar2<<endl;
    cout << "T_bar1: "<<pV.T_bar1<<endl;
    cout << "T_bar2: "<<pV.T_bar2<<endl;
    cout << "v_bar: "<<pV.v_bar<<endl;
    cout << "v_maxD: "<<pV.v_maxD<<endl;
    cout << "v_maxDC: "<<pV.v_maxDC<<endl;
}