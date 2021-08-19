#pragma once

#include "findTimeOnPath.cpp"
#include "pathIntersections.cpp"
#include "findPosOnPath.cpp"
#include "findCoordOnPath.cpp"
vec trajIntersection(path_elem path1, path_elem path2, pathVel_elem pathVel1, pathVel_elem pathVel2, interSec_elem interSec){
    vec Si11 = interSec.Si11;
    vec Si12 = interSec.Si12;
    vec Si21 = interSec.Si21;
    vec Si22 = interSec.Si22;

    vec Ti11 = findTimeOnPath(Si11, path1, pathVel1);
    vec Ti12 = findTimeOnPath(Si12, path1, pathVel1);
    vec Ti21 = findTimeOnPath(Si21, path2, pathVel2);
    vec Ti22 = findTimeOnPath(Si22, path2, pathVel2);
    vec Ti;
    if(! (Ti11(0) <= Ti21(0) || Ti22(0) <= Ti11(0)) )
    {
        Ti = {Ti11(0), Ti12(0), Ti21(0), Ti22(0)};
        vec ascend_Ti = arma::sort(Ti);
        vec Ti11_bar(1);
        Ti11_bar(0) = ascend_Ti(1);
        vec Ti12_bar(1);
        Ti12_bar(0) = ascend_Ti(2);
        mat Si11_bar=findPosOnPath(Ti11_bar,path1,pathVel1);
        mat Si12_bar=findPosOnPath(Ti12_bar,path1,pathVel1);

        vec Si1_bar;
        Si1_bar = arma::regspace(Si11_bar(0), 2.0*rho_D, Si12_bar(0));
        Si1_bar.resize(Si1_bar.n_elem+1);
        Si1_bar(Si1_bar.n_elem-1) = Si12_bar(0);
        CoorOnPath ri1_bar = findCoordOnPath(Si1_bar, path1);
        vec Ti1_bar = findTimeOnPath(Si1_bar, path1, pathVel1);

        int nP = Si1_bar.n_elem;

        Ti.reset();
        Ti = zeros<vec>(2);

        int i = 0;
        while (i < nP)
        {
            i++;
            mat Si2_bar = findPosOnPath(Ti1_bar.subvec(i,i), path2, pathVel2);
            CoorOnPath ri2_bar = findCoordOnPath(Si2_bar, path1);
            double d= arma::norm(ri1_bar.rp.col(i) - ri2_bar.rp);
            if(d <- 2*rho_D)
            {
                Ti(0) = Ti1_bar(i);
            }
            while (d <= 2*rho_D && i < nP)
            {
                i++;
                Si2_bar=findPosOnPath(Ti1_bar.subvec(i,i),path2,pathVel2);
                ri2_bar=findCoordOnPath(Si2_bar,path2);
                d = arma::norm(ri1_bar.rp.col(i) - ri2_bar.rp);
            }
            Ti(2) = Ti1_bar(i);
        }

    } else 
    {
        Ti = zeros<vec>(2);
    }
    return Ti;
}