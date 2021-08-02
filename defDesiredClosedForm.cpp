#include "AllParametersExperiment.hpp"
#include <armadillo>

using namespace arma;

void defDesiredClosedForm(mat XDFc, double RDF0,double phi0, mat XA, int Na, int ND, int flagNDeven, double delta_t){
    mat rDFc, vDFc;
    rDFc = XDFc.submat(0,0,1,0);
    vDFc = XDFc.submat(2,0,3,0);

    mat uDFc0;
    if (arma::norm(rDFc - rS) < 3*rho_sn)
    {
        uDFc0 = -kDFphir*(rDFc-rS)-kDFphiv*(vDFc);
    } else {
        uDFc0 = -kDFphir*(rDFc-rS);
    }
    
}