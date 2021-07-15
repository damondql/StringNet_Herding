#include <iostream>
#include <armadillo>
#include <vector>
#include <filesystem>
using namespace std;
using namespace arma;


struct obs
{
    double rho_safe;
    cube rVO;
    int NO;
    mat rCO2;
    cube rVO1;
};

struct tanG
{
    int NO;
    mat nVO;
    // xOk 1x2 cells containing funcion handle
    // yOk 1x2 cells containing funcion handle
    // dxOk_dg 1x2 cells containing funcion handle
    // dyOk_dg 1x2 cells containing funcion handle
    cube GammaOk;
    cube PeriSOk;
    mat PeriO;
    mat G;
    mat G_pathType;
    cube rVO;
    cube rVO2;
    mat gTO_all;
    mat rTO_all;
    mat obsId_all;
    mat vertId_all;
    mat obsVertId_all;
    mat obsVertPos_all;
    mat rCO2;
    double E_m;
};

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

int fileCount(std::filesystem::path p1) {
    int rVO_count { };

    for (auto& p : std::filesystem::directory_iterator(p1))
    {
        ++rVO_count;
    }
    return rVO_count;
}



obs ExpData_obs;
void loadObs()
{
    ExpData_obs.rCO2.load("./tempDataExp/rCO2.txt");
    ExpData_obs.rCO2.print("rCO2:");
    std::ifstream NO("./tempDataExp/NO.txt");
    NO >> ExpData_obs.NO;
    std::filesystem::path p1 { "./tempDataExp/rVO" };
    int rVO_count = fileCount(p1);
    arma::mat tempMat;
    tempMat.load("./tempDataExp/rVO/rVO_1.txt");
    tempMat.print("rVO_1:");
    ExpData_obs.rVO.resize(tempMat.n_rows, tempMat.n_cols,rVO_count);
    tempMat.load("./tempDataExp/rVO/rVO_1.txt");
    ExpData_obs.rVO.slice(0) = tempMat;
    cout << "loaded 1" << endl;
    tempMat.load("./tempDataExp/rVO/rVO_2.txt");
    ExpData_obs.rVO.slice(1) = tempMat;
    cout << "loaded 2" << endl;
    ExpData_obs.rVO.print();
    for (int i = 1; i <= 2; i++)
    {
        string roV_n = "./tempDataExp/rVO/rVO_";
        roV_n.append(".txt");
        cout << roV_n << endl;
    }
    
}



int main()
{
    loadObs();


}