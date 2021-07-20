#include "findShortestPath.cpp"
#include "AllParametersExperiment.hpp"

struct interSec_elem
{
    mat flag;
    double Pbar1;
    double Pbar2;
    int flag0;
    mat Pos1;
    mat Pos2;
    mat S10;
    mat S20;
    vec rV11;
    vec rV12;
    vec rV21;
    vec rV22;
    double S11;
    double S12;
    double S21;
    double S22;
    double Si11;
    double Si12;
    double Si21;
    double Si22;
};

interSec_elem pathIntersections(path_elem Path1, path_elem Path2) {
    interSec_elem result;
    int NS1 = 1;
    int NS2 = 1; 
    result.flag = zeros<mat>(NS1, NS2);
    result.Pbar1 = 0;
    result.Pbar2 = 0;
    vec rV1 = Path1.rV;
    vec rV2 = Path2.rV;
    int segType1 = Path1.segType;
    int segType2 = Path2.segType;
    result.flag0 = 0;
    vec r1(2);
    vec r2(2);
    vec dr;

    for (int i = 0; i < NS1; i+=2)
    {
        double Si0 = Path1.S(i);
        r1 = rV1.col(i);
        r2 = rV2.col(i+1);
        dr = r2 - r1;
        double L1 = arma::norm(dr,2);
        double drx = dr(0);
        double dry = dr(1);
        double theta1 = atan2(dry, drx);
        if (theta1 < 0)
        {
            theta1 = theta1 + 2 * M_PI;
        }
        double mL1 = tan(theta1);
        double cL1 = r1(1) - mL1 * r1(0);
        vec dr_hat = dr/L1;

        for (int j = 0; j < NS2; j+=2)
        {
            double Sj0 = Path2.S(j);
            vec r21 = rV2.col(j);
            vec r22 = rV2.col(j+1);
            double drx2 = r22(0) - r21(1);
            double dry2 = r22(1) - r21(1);
            double theta2 = atan2(dry2, drx2);
            if(theta2 < 0) {
                theta2 = theta2 + 2 * M_PI;
            }
            double dtheta = theta2 - theta1;
            if (dtheta < 0) {
                dtheta = dtheta + 2 * M_PI;
            }
            double mL2 = tan(theta2);
            double cL2 = r21(1) - mL2*r21(0); 
            double L2 = norm(r21-r22);


        }
        

    }
    

}
