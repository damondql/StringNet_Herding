#pragma once

#include "findTimeOnPath.cpp"
#include "pathIntersections.cpp"
double leadTimeToAvoidCollision(path_elem path1, path_elem path2, pathVel_elem pathVel1, pathVel_elem pathVel2, interSec_elem interSec, int leadId) {
    vec Si11 = interSec.Si11;
    vec Si12 = interSec.Si12;
    vec Si21 = interSec.Si21;
    vec Si22 = interSec.Si22;

    vec Ti11 = findTimeOnPath(Si11, path1, pathVel1);
    vec Ti12 = findTimeOnPath(Si12, path1, pathVel1);
    vec Ti21 = findTimeOnPath(Si21, path2, pathVel2);
    vec Ti22 = findTimeOnPath(Si22, path2, pathVel2);

    

}