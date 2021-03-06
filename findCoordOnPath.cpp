#pragma once
#include "AllParametersExperiment.cpp"
#include "findShortestPath.cpp"
#include <armadillo>
#include <istream>

using namespace std;
using namespace arma;

struct CoorOnPath
{
    arma::mat rp;
    arma::mat thetap;
};

CoorOnPath findCoordOnPath(double S, path_elem p) {
    CoorOnPath result;
    result.rp.resize(2,1);
    result.thetap.resize(2,1);
    uvec ind0 = arma::find(p.S < S);
    int segld;
    if(!ind0.is_empty()) {
        segld = ind0(ind0.n_elem-1);
    } else {
        segld = 0;
    }
    double S0 = p.S(segld);
    double S1 = p.S(segld+1);
    double dS = S - S0;
    mat r11(p.rV.n_rows,1, fill::zeros);
    mat r12(p.rV.n_rows,1, fill::zeros);
    r11 = p.rV.submat(0,segld,p.rV.n_rows-1 , segld);
    r12 = p.rV.submat(0,segld+1,p.rV.n_rows-1, segld+1);
    if ( (segld+1) % 2 == 1) {
        double lambda = dS / (S1 - S0);
        result.rp.submat(0,0,result.rp.n_rows-1, 0) =  (1-lambda) * r11 + lambda * r12;
        result.thetap = atan2(r12(1)-r11(1), r12(0)-r11(0));
    }
    return result;
}
