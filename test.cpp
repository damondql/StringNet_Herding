#include <armadillo>
#include <vector>
#include <complex>
using namespace std;
using namespace arma;


struct tanG
{
    int NO;
    arma::mat nVO;
    // xOk 1x2 cells containing funcion handle
    // yOk 1x2 cells containing funcion handle
    // dxOk_dg 1x2 cells containing funcion handle
    // dyOk_dg 1x2 cells containing funcion handle
    cube GammaOk;
    cube PeriSOk;
    arma::mat PeriO;
    arma::mat G;
    arma::mat G_pathType;
    cube rVO;
    cube rVO2;
    arma::mat gTO_all;
    arma::mat rTO_all;
    arma::mat obsId_all;
    arma::mat vertId_all;
    arma::mat obsVertId_all;
    arma::mat obsVertPos_all;
    arma::mat rCO2;
    double E_m;
};

cube rVO(2,4,2);
double rho_safe = 2.0118;
tanG tangentGraph(cube rVO) {
    rVO.slice(0) = {{60,160,160,60},
                    {310,310,390,390}};
    rVO.slice(1) = {{-300,-190,-190,-300},
                    {350,350,560,560}};
    
    tanG result;
    
    int NO = rVO.n_slices;
    
    arma::rowvec PeriO(NO);
    vec nVO(NO);
    arma::mat posVO(rVO.n_rows, rVO.n_cols + 1, fill::zeros);
    arma::mat posVO2(rVO.n_rows, 2 * (rVO.n_cols), fill::zeros);
    int posVO2_cols = 2*rVO.n_cols;
    arma::mat PeriSOk(posVO2_cols, NO, fill::zeros);
    arma::mat AngOk(posVO2_cols, NO, fill::zeros);
    arma::mat GammaOk(posVO2_cols+1, NO, fill::zeros);
    arma::mat GammaCOk(posVO2_cols, NO, fill::zeros);
    rVO.resize(rVO.n_rows, rVO.n_cols + 1, rVO.n_slices);
    arma::cube rVO2(posVO2.n_rows, posVO2.n_cols, NO);
    for (size_t k = 0; k < NO; k++)
    {   
        nVO(k) = rVO.n_cols -1;
        rVO.slice(k).submat(0,rVO.n_cols-1, rVO.n_rows-1, posVO.n_cols-1) = rVO.slice(k).submat(0,0,posVO.n_rows-1,0);
        posVO = rVO.slice(k);
        posVO2.resize(rVO.n_rows, posVO2_cols);
        int countV = 0;
        double xip = 0;
        double yip = 0;
        for (int i = 0; i < nVO(k); i++)
        {
            double mLi = (posVO(1,i+1)-posVO(1,i))/(posVO(0,i+1)-posVO(0,i));
            // cout << "mLi: " << mLi << endl;
            double xi1, xi2, yi1, yi2,dx;
            if (mLi ==0)
            {
                xi1 = posVO(0,i);
                xi2 = posVO(0,i);
                yi1 = posVO(1,i) + rho_safe;
                yi2 = posVO(1,i) - rho_safe;
            } else if (isinf(mLi)) {
                xi1=posVO(0,i)+rho_safe;
                xi2=posVO(0,i)-rho_safe;
                yi1=posVO(1,i);
                yi2=posVO(1,i);
            }else {
                dx = sqrt(pow(rho_safe,2)/(1+1/pow(mLi,2)));
                xi1=posVO(0,i)+dx;
                xi2=posVO(0,i)-dx;
                yi1=posVO(1,i)+dx*(-1/mLi);
                yi2=posVO(1,i)-dx*(-1/mLi);
            }
            vec a = {posVO(0,i+1) - posVO(0,i),posVO(1,i+1)-posVO(1,i), 0};
            vec b = {xi1 - posVO(0,i), yi1 - posVO(1,i), 0};
            vec crossProd = arma::cross(a,b);
            countV++;

            if (crossProd(2) < 0) {
                posVO2.submat(0,countV-1,1,countV-1) = {xi1,yi1};
                countV++;
                if (mLi == 0)
                {
                    xip=posVO(0,i+1);
                    yip=posVO(1,i+1)+rho_safe;
                } else if (isinf(mLi)) {
                    xip=posVO(0,i+1)+rho_safe;
                    yip=posVO(1,i+1);   
                } else {
                    xip=posVO(0,i+1)+dx;
                    yip=posVO(1,i+1)+dx*(-1/mLi);
                }
            } else {
                posVO2.submat(0, countV-1, 1,countV-1) = {xi2, yi2};
                countV++;
                if (mLi == 0)
                {
                    xip=posVO(0,i+1);
                    yip=posVO(1,i+1)-rho_safe;
                } else if (isinf(mLi)) {
                    xip=posVO(0,i+1)-rho_safe;
                    yip=posVO(1,i+1);
                } else {
                    xip=posVO(0,i+1)-dx;
                    yip=posVO(1,i+1)-dx*(-1/mLi);
                }
                
            }
            posVO2.submat(0,countV-1,1,countV-1) = {xip, yip};
            // posVO2.print("posV02");
        }
        arma::mat tempM;
        tempM = posVO2;
        posVO2.submat(0,0,posVO2.n_rows-1,0) = tempM.submat(0,posVO2.n_cols-1, posVO2.n_rows-1, posVO2.n_cols-1);
        posVO2.submat(0,1,posVO2.n_rows-1, posVO2.n_cols-1) = tempM.submat(0,0,posVO2.n_rows-1, posVO2.n_cols-2);
        rVO2.slice(k) = posVO2;
        int nV2 = countV;
        posVO2.resize(posVO2.n_rows, posVO2.n_cols + 1);
        posVO2.submat(0, posVO2.n_cols-1, posVO2.n_rows-1, posVO2.n_cols-1) = posVO2.submat(0,0,posVO2.n_rows-1, 0);
        for (int ii = 1; ii <= nVO(k); ii++)
        {
            double xV = posVO(0,ii-1);
            double yV = posVO(1,ii-1);
            double ang1 = atan2(posVO2(1,2*ii-2)-yV,posVO2(0,2*ii-2)-xV);
            if (ang1 < 0)
            {
                ang1 += 2*M_PI;
            }
            double ang2 = atan2(posVO2(1,2*ii-1)-yV,posVO2(0,2*ii-1)-xV);
            if (ang2 < 0)
            {
                ang2 += 2*M_PI;
            }
            if (ang2 < ang1) {
                ang2 += 2*M_PI;
            }
            PeriSOk(2*ii-2,k) = rho_safe * (ang2 - ang1);
            std::complex<double> comp = {posVO2(0, 2*ii-1) - posVO2(0, 2*ii), posVO2(1,2*ii-1) - posVO2(1,2*ii)};
            PeriSOk(2*ii-1,k)= sqrt(norm(comp));
            AngOk(2*ii-2, k) = ang1;
            AngOk(2*ii-1, k) = ang2; 
        }
        PeriO = arma::sum(PeriSOk);
        GammaOk(0,k) = 0;
        for (int ii = 1; ii <= nVO(k); ii++) {
            GammaOk(2*ii-1,k) = GammaOk(2*ii-2,k) + PeriSOk(2*ii-2,k)/PeriO(k);
            GammaOk(2*ii,k) = GammaOk(2*ii-1,k) + PeriSOk(2*ii-1,k)/PeriO(k);
            GammaCOk(2*ii-2,k) = AngOk(2*ii-2,k)/2/M_PI;
            GammaCOk(2*ii-1,k) = AngOk(2*ii-1,k)/2/M_PI;
        }        
        for (int ii = 1; ii <= nVO(k); ii++)
        {
            int i1 = 2*ii-1;
            int i2 =  2*ii;
            int i3 = 2*ii + 1;
            double xV = posVO(0, ii-1);
            double yV = posVO(1, ii-1);
            double gk1 = GammaOk(i1-1, k);
            double gk2 = GammaOk(i2-1, k);
            double gk3 = GammaOk(i3-1, k);
            double ang1 =  AngOk(2*ii-2,k);
            double ang2 = AngOk(2*ii-1,k); 
        }
    }
    mat A(2*NO, NO, fill::zeros);
    A.submat(0,0, NO -1, NO-1) = arma::eye(2,2);
    A.submat(NO,0, 2*NO -1, NO-1) = -1 * arma::eye(2,2);
    A.print("A:");
    vec gamma = arma::regspace(0,0.01,1);
    int countP = 0;
    int countTotT = 0;
    for (int k = 1; k <= NO -1; k++)
    {
        for (int kk = k+1; kk <= NO; kk++)
        {
            countP++;
        }
        
    }
    



    return result;
}

int main() {
    tanG a = tangentGraph(rVO);
}