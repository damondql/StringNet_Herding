#pragma once
#include <armadillo>

struct output {
    mat W;
    mat Rij_tilde;
};

output findCommGraphAndFormDist(double N, double shapeld, double R0);