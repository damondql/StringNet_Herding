#include <iostream>
#include <armadillo>
#include <complex>


using namespace std;
using namespace arma;

struct path
{
    arma::mat rV;
    double P;
    arma::vec S;
    int NS = 1;
    int segType = 3;
    arma::mat rVc;
};

path findPath(arma::vec rI, arma::vec rF) {
    path Path;
    Path.rV.resize(rI.n_rows,2);
    Path.rV.submat(0,0,rI.n_rows-1,0) = rI;
    Path.rV.submat(0,1,rI.n_rows-1,1) = rF;
    std::complex<double> comp = {rI(0)-rF(0), rI(1) - rF(1)};
    Path.P = sqrt(norm(comp));
    Path.S = {0, Path.P};
    return Path;
}

int main() {
    arma:vec a(2);
    a = {0,0};
    arma::vec b(2);
    b = {5,6};
    path p = findPath(a,b);
    p.rV.print("rV:");
    cout << "P: " << p.P << endl;
    p.S.print("S: ");
    cout << "NS: " << p.NS << endl;
    cout << "segType: " << p.segType << endl;
}