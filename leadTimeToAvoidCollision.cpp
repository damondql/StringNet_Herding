#pragma once

#include "findTimeOnPath.cpp"
#include "pathIntersections.cpp"
#include "trajInterSection.cpp"

double leadTimeToAvoidCollision(path_elem path1, path_elem path2, pathVel_elem pathVel1, pathVel_elem pathVel2, interSec_elem interSec, int leadId) {
    vec Si11 = interSec.Si11;
    vec Si12 = interSec.Si12;
    vec Si21 = interSec.Si21;
    vec Si22 = interSec.Si22;

    vec Ti11 = findTimeOnPath(Si11, path1, pathVel1);
    vec Ti12 = findTimeOnPath(Si12, path1, pathVel1);
    vec Ti21 = findTimeOnPath(Si21, path2, pathVel2);
    vec Ti22 = findTimeOnPath(Si22, path2, pathVel2);
    double leadTime;
    if (leadId == 1)
    {
        if(Ti12(0) > Ti21(0))
        {
            pathVel_elem pathVel20 = pathVel2;
            double t1 = 0;
            vec t2 = Ti12 - Ti21;
            vec Ti0 = trajIntersection(path1, path2, pathVel1, pathVel20, interSec);
            if(!Ti0.is_empty())
            {
                leadTime = t2(0) - t1;
                pathVel20 = pathVel2;
                pathVel20.T = pathVel20.T + leadTime;
                pathVel20.T_bar1 = pathVel20.T_bar1+leadTime;
                pathVel20.T_bar2 = pathVel20.T_bar2 + leadTime;
                vec Ti = trajIntersection(path1,path2, pathVel1, pathVel20, interSec);
                while ((arma::norm(Ti-Ti0) > 1e-5 || arma::norm(Ti) < 1e-5) && leadTime > 1e-5)
                {
                    double dTi0 = Ti0(1) - Ti0(0);
                    double dTi = Ti(1) - Ti(0);
                    if(dTi <= dTi0)
                    {
                        t2 =(t1+t2(0)) / 2;
                        Ti0 = Ti;
                        leadTime= t2(0) - t1;
                        pathVel20=pathVel2;
                        pathVel20.T=pathVel20.T+leadTime;
                        pathVel20.T_bar1=pathVel20.T_bar1+leadTime;
                        pathVel20.T_bar2=pathVel20.T_bar2+leadTime;
                    } else if (dTi > dTi0)
                    {
                        t1 = (t1+t2(0)) / 2;
                        leadTime=t2(0)-t1;
                        pathVel20=pathVel2;
                        pathVel20.T=pathVel20.T+leadTime;
                        pathVel20.T_bar1=pathVel20.T_bar1+leadTime;
                        pathVel20.T_bar2=pathVel20.T_bar2+leadTime;
                        Ti=trajIntersection(path1,path2,pathVel1,pathVel20,interSec);
                    }
                }
                
            } else 
            {
                leadTime = 0;
            }

 
        } else
        {
            leadTime = 0;
        }

    } else
    {
        if(Ti22(0) > Ti11(0))
        {
            pathVel_elem pathVel10 = pathVel1;
            double t1 = 0;
            vec t2 = Ti22 - Ti11;
            vec Ti0 = trajIntersection(path1, path2, pathVel10, pathVel2, interSec);
            if (!Ti0.is_empty())
            {
                leadTime = t2(0) - t1;
                pathVel10=pathVel1;
                pathVel10.T=pathVel10.T+leadTime;
                pathVel10.T_bar1=pathVel10.T_bar1+leadTime;
                pathVel10.T_bar2=pathVel10.T_bar2+leadTime;
                vec Ti=trajIntersection(path1,path2,pathVel10,pathVel2,interSec);
                while((arma::norm(Ti-Ti0)>1e-5 || arma::norm(Ti)<1e-5) && leadTime>1e-5)
                {
                    double dTi0=Ti0(1)-Ti0(0);
                    double dTi=Ti(1)-Ti(0);
                    if (dTi<=dTi0)
                    {
                        t2=(t1+t2(0))/2;
                        Ti0=Ti;
                        
                        leadTime=t2(0)-t1;
                        pathVel10=pathVel1;
                        pathVel10.T=pathVel10.T+leadTime;
                        pathVel10.T_bar1=pathVel10.T_bar1+leadTime;
                        pathVel10.T_bar2=pathVel10.T_bar2+leadTime;
                        Ti=trajIntersection(path1,path2,pathVel10,pathVel2,interSec);
                    } else if (dTi>dTi0)
                    {
                        t1=(t1+t2(0))/2;
                        Ti0=Ti;
                        
                        leadTime=t2(0)-t1;
                        pathVel10=pathVel1;
                        pathVel10.T=pathVel10.T+leadTime;
                        pathVel10.T_bar1=pathVel10.T_bar1+leadTime;
                        pathVel10.T_bar2=pathVel10.T_bar2+leadTime;
                        Ti=trajIntersection(path1,path2,pathVel10,pathVel2,interSec);
                    }
                }

            } else
            {
                leadTime = 0;
            }
        } else {
            leadTime = 0;
        }


    }
    return leadTime;

}