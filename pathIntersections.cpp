#include "findShortestPath.cpp"
#include "AllParametersExperiment.hpp"
#include "interSecLineLine.cpp"

struct interSec_elem
{
    mat Flag;
    double Pbar1;
    double Pbar2;
    int flag0;
    arma::field<mat> Pos1;
    arma::field<mat> Pos2;
    arma::field<vec> S10;
    arma::field<vec> S20;
    vec rV11;
    vec rV12;
    vec rV21;
    vec rV22;
    double S11;
    double S12;
    double S21;
    double S22;
    vec Si11;
    vec Si12;
    vec Si21;
    vec Si22;
};

interSec_elem pathIntersections(path_elem Path1, path_elem Path2) {
    interSec_elem result;
    int NS1 = 1;
    int NS2 = 1; 
    result.Flag = zeros<mat>(NS1, NS2);
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
        vec dr_hat = dr / L1;

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
            double L2 = arma::norm(r21-r22);

            lineIntersec interSec0 = interSecLineLine(r1,r2,mL1,cL1,theta1,drx,dry,L1,r21,r22,mL2,cL2,theta2,drx2,dry2,L2,dtheta);
            if (interSec0.Flag != 0)
            {
                result.Flag(i,j) = interSec0.Flag;
                result.Pos1(i,j) = interSec0.Pos1;
                result.Pos2(i,j) = interSec0.Pos2;
                result.Pbar1 = result.Pbar1 + interSec0.Pbar1;
                result.Pbar2 = result.Pbar2 + interSec0.Pbar2;
                result.S10(i,j) = Si0 + interSec0.S10;
                result.S20(i,j) = Sj0 + interSec0.S20;
            }
            
        }
        

    }

    for (int j = 0; j < NS2; j+=2)
    {
        double Sj0 = Path2.S(j);
        r1 = rV2.col(j);
        r2 = rV2.col(j+1);
        dr = r2 - r1;
        double L1 = arma::norm(dr);
        double drx = r2(0) - r1(0);
        double dry = r2(1) - r1(1);
        double theta1 = atan2(dry,drx);
        if (theta1 < 0) {
            theta1 = theta1 + 2 * M_PI;
        }
        double mL1 = tan(theta1);
        double cL1 = r1(1) - mL1 * r1(0);
        double dr_hat = dr / (arma::norm(dr));


    }
    int countS1 = -1;
    for (int i = 0; i < NS1; i++)
    {
        for (int j = 0; j < NS2; j++)
        {
            if (result.Flag(i,j) != 0)
            {
                countS1++;
                result.rV11.col(countS1) = result.Pos1(i,j).col(0);
                result.rV12.col(countS1) = result.Pos1(i,j).col(1);
                result.rV21.col(countS1) = result.Pos2(i,j).col(0);
                result.rV22.col(countS1) = result.Pos2(i,j).col(1);
                result.S11 = result.S10(i,j)(0);
                result.S12 = result.S10(i,j)(1);
                result.S21 = result.S20(i,j)(0);
                result.S22 = result.S20(i,j)(1);
            }
            
        }
        
    }
    
    if (countS1 > -1)
    {
        result.flag0 = 1;
        // sort the intervals//
        vec Si11, Si12, Si21, Si22;
        Si11.resize(1);
        Si12.resize(1);
        Si21.resize(1);
        Si22.resize(1);
        int starti = 0;
        int endi = 0;
        int k = 0;
        Si11(k) = result.S11;
        Si12(k) = result.S12;


        //while loop//

        result.Si11 = Si11;
        result.Si12 = Si12;
        
        Si21(k) = result.S21;
        Si22(k) = result.S22;
        result.Si21 = Si21;
        result.Si22 = Si22;

    }
    
    
    return result;
}
