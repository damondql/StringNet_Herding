#pragma once
#include <iostream>
#include <armadillo>
#include <complex>


using namespace std;
using namespace arma;

struct path_elem
{
    arma::mat rV;
    double P;
    arma::vec S;
    int NS = 1;
    int segType = 3;
    arma::mat rVc;
};

path_elem findShortestPath(arma::vec rI, arma::vec rF) {
    path_elem resultP;
    resultP.rV.resize(rI.n_rows,2);
    resultP.rV.submat(0,0,rI.n_rows-1,0) = rI;
    resultP.rV.submat(0,1,rI.n_rows-1,1) = rF;
    std::complex<double> comp = {rI(0)-rF(0), rI(1) - rF(1)};
    resultP.P = sqrt(norm(comp));
    resultP.S = {0, resultP.P};
    // resultP.rV.print("rV:");
    // resultP.S.print("S:");
    return resultP;
}
