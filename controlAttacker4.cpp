#pragma once
#include "AllParameters.cpp"
#include <armadillo>
#include <math.h>
#include "helperFunction.cpp"
#include "inhull.cpp"
// using namespace std;
using namespace arma;

struct control_attacker_t
{
    arma::mat uA;
    mat uA0;
    double R_AO_min;
    double R_AAProjS_min;
    mat vA_des;
    mat vA_des_dot;
    mat SigmaProdD;
    mat F_A;
    mat F_A_dot;
};


control_attacker_t controlAttacker4(mat XA, mat XA_goal, mat XA_goal_dot,mat leaderIDA,
                      mat XD, mat WA, mat R_tilde_AA,
                      mat WDString,int NA, int ND, mat clusteridA, field<vec>indAinClusterA, vec NAinClusterA, mat rhoA_con,vec flagAEnclosed , std::vector<vec> assign, int ti, int bound){
    // XA.print("XA: ");
    // XA_goal.print("XA_goal: ");
    // XA_goal_dot.print("XA_goal_dot: ");
    // XD.print("XD: ");
    // WA.print("WA: ");
    // WDString.print("WDString: ");
    // cout << "enter function: " << endl;
    double tol = 0.5;
    int NO = 0;

    mat rD, vD;
    if (ND > 0)
    {
        rD.resize(2,ND);
        vD.resize(2,ND);
        rD = XD.submat(0,0,1,ND-1);
        vD = XD.submat(2,0,3,ND-1);
    }
    
    // arma::colvec rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1)/NA;
    // arma::colvec vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1)/NA;
    // cout << "rAcm initializing " << endl;
    mat rAcm(2, NAinClusterA.n_elem, fill::zeros);
    mat vAcm(2, NAinClusterA.n_elem, fill::zeros);
    // cout << "rAcm initialization finished " << endl;
    for (int c = 0; c < NAinClusterA.n_elem; c++)
    {
        int NAinC = NAinClusterA(c);
        for (int ii = 0; ii < NAinC; ii++)
        {
            rAcm.col(c) = rAcm.col(c) + XA.submat(0,indAinClusterA(c)(ii), 1, indAinClusterA(c)(ii));
            vAcm.col(c) = vAcm.col(c) + XA.submat(2,indAinClusterA(c)(ii), 3, indAinClusterA(c)(ii));
        }
        rAcm = rAcm / NAinClusterA(c);
        vAcm = vAcm / NAinClusterA(c);
        
        
    }
    // rAcm.print("rAcm: ");
    // vAcm.print("vAcm: ");

    
    
    mat SigmaProdD(NA,1,fill::ones);
    double R_AO_min = INFINITY;
    double R_AAProjS_min = INFINITY;

    mat F_A(2,NA,fill::zeros);
    mat F_A_dot(2,NA,fill::zeros);
    mat F_AD(2,NA,fill::zeros);
    mat F_AD_dot(2,NA,fill::zeros);
    mat uA;
    mat uA0(2,NA,fill::zeros);
    mat vA_dot;
    mat vA_Des;
    mat vA_Des_dot;
    // cout << "before for loop" << endl;
    for (int i = 0; i < NA; i++)
    {
        mat WDString_temp = WDString;
        mat rA = XA.submat(0,i,1,i);
        mat vA = XA.submat(2,i,3,i);

        if (arma::norm(vA) > v_maxA(i))
        {
            vA = vA * v_maxA(i) / arma::norm(vA);
        }
        
        mat rA_goal = XA_goal.submat(0,i,1,i);
        mat vA_goal = XA_goal.submat(2,i,3,i);

        double digmaProd = 1;
        double sigmaBarProd = 1;
        double sigmaSum = 0;
        double sigmaBarSum = 0;
        mat Sigma = zeros<mat>(NO+ND+NA,1);
        mat Sigma_dot = Sigma;

        vec uAOv = zeros<vec>(2);
        vec uAOr = uAOv;

        vec uAFv = zeros<vec>(2);
        vec uAFr = uAFv;

        uvec find_result = arma::find(WA.row(i) == 1);
        double sigma;
        vec nabla_ri_Vii;
        //check formation controller
        for (int j = 0; j < find_result.n_elem; j++)
        {
            int ii = find_result(j);
            double Rii0 = R_tilde_AA(i,ii) - R_m_AA;
            double Ri_ii = arma::norm(rA - XA.submat(0,ii, 1, ii));
            if(Ri_ii < R_u_AA)
            {
                if(Ri_ii < R_bar_AA)
                {
                    sigma = 1;
                } else if(Ri_ii > R_bar_AA && Ri_ii < R_u_AA)
                {
                    sigma = A_A_A*pow(Ri_ii,3)+B_A_A*pow(Ri_ii,2)+C_A_A*Ri_ii+D_A_A;
                }
            } else 
            {
                sigma = 0;
            }
            if(Ri_ii - R_m_AA > 1e-1)
            {
                nabla_ri_Vii = kAFr*WA(i,ii)*(rA-XA.submat(0,ii,1,ii))/Ri_ii/abs(Ri_ii-R_m_AA)*(pow((Ri_ii-R_m_AA),2)-pow(Rii0,2))/(pow((Ri_ii-R_m_AA),2)+pow(Rii0,2));
            } else
            {
                nabla_ri_Vii=-kAFr*WA(i,ii)*(rA-XA.submat(0,ii,1,ii))/Ri_ii*largeP;
            }
            vec dv;
            if(arma::norm(vA - XA.submat(2,ii,3,ii)) > 1e-16)
            {
                dv = kAFv*(vA-XA.submat(2,ii,3,ii))*pow(arma::norm(vA-XA.submat(2,ii,3,ii)),(alphaAFv-1));
            } else
            {
                dv = zeros<vec>(2);
            }
            uAFv = uAFv - sigma*dv;
            uAFr = uAFr - sigma*nabla_ri_Vii;
        }
        // uAFv.print("uAFv: ");
        // uAFr.print("uAFr: ");



        //check for nearby defenders
        // cout << "111111111111111conA" << endl;
        double sigmaSumD = 0;
        double sigmaProdD = 1;
        double minRAD = INFINITY;
        vec uADv(2,fill::zeros);
        vec uADr(2,fill::zeros);
        
        int countAPS = 0;
        vec uAD_pot(2);
        // cout << "22222222222222222222conA" << endl;
        if (i > NA - NA_sep)
        {
            // cout << "enter if " << endl;
            // double R_m = 15*R_m_AD;
            // double R_underbar = R_m+20;
            // double R_bar = R_m + 25;
            mat sigma_param(4,1);
            sigma_param(0) = A_A_D2;
            sigma_param(1) = B_A_D2;
            sigma_param(2) = C_A_D2;
            sigma_param(3) = D_A_D2;
            vec potentialControl_result = potentialControl(1e-5, XA.col(i),XD,
                                                           rho_c_A,sigma_param,
                                                           R_m_AD2,R_bar_AD2,R_u_AD2, 2*Rij0(0),kADr,kADv, alphaADv);
            uAD_pot = potentialControl_result.subvec(0,1);
            minRAD = potentialControl_result(3);
        } else {
            // cout << "enter else" << endl;
            mat sig_para(4,1);
            sig_para(0,0) = A_A_D;
            sig_para(1,0) = B_A_D;
            sig_para(2,0) = C_A_D;
            sig_para(3,0) = D_A_D;
            // XD.submat(0,0,XD.n_rows-1, ND-1).print("potential input:");
            vec potentialControl_result = potentialControl(1e-5, XA.col(i),XD,
                                                           rho_c_A,sig_para,
                                                           R_m_AD,R_bar_AD, R_u_AD, Rij0(0),kADr,kADv, alphaADv);
            // potentialControl_result.print("finish potential Control:");
            uAD_pot = potentialControl_result.subvec(0,1);
            minRAD = potentialControl_result(2);
        }
        // uAD_pot.print("uAD_pot: ");
        // cout << "minRAD: " << minRAD << endl;
        // cout << "3333333333333333333conA" << endl;
        // uAD_pot.print("uAD_pot: ");
        // cout << "minRAD: " << minRAD << endl;
        mat rAProjS;
        mat vAProjS;
        mat vAProjS0, vAProjS1;
        for(int j = 0; j < ND; j++)
        {
            // cout << "j: " << j << endl;
            uvec find_result2 = arma::find(WDString_temp.row(j) == 1);
            // find_result2.print("find: ");
            if (!find_result2.is_empty())
            {
                
                for (int ii = 0; ii < find_result2.n_elem; ii++)
                {
                    int jj = find_result2(ii);
                    countAPS++;
                    vec projection_resutl = projectionOnLine(rA, XD.submat(0,jj,1,jj), XD.submat(0,j,1,j));
                    rAProjS.insert_cols(countAPS-1, projection_resutl.subvec(0,1));
                    mat rTDD = XD.submat(0,j,1,j) - XD.submat(0,jj,1,jj);
                    rTDD = rTDD / arma::norm(rTDD);
                    mat tempM;
                    tempM = rTDD.t() * vA;
                    // cout << "3.1" << endl;
                    if (tempM(0,0) < 0)
                    {
                        rTDD = -rTDD;
                    }
                    
                    // rTDD.print("rTDD: ");
                    // vA_temp.print("vA_temp: ");
                    tempM = rTDD.t() * vA;
                    // cout << "countAPS: " << countAPS << endl;
                    // vAProjS0.print("vAProjS0: ");
                    vAProjS0.insert_cols(countAPS-1, tempM(0,0) * rTDD);
                    // vAProjS0.print("vAProjS0: ");
                    // cout<< "6.1" << endl;
                    vec tempV;
                    tempV = XD.submat(2,j,3,j) + 
                               (XD.submat(2,jj,3,jj) - XD.submat(2,j,3,j)) 
                               * arma::norm(rAProjS.submat(0,countAPS-1,1,countAPS-1) - XD.submat(0,j,1,j))
                               / arma::norm(XD.submat(0,jj,1,jj) - XD.submat(0,j,1,j));
                    vAProjS1.insert_cols(countAPS-1, tempV);
                    tempV = vAProjS0.col(countAPS-1) + vAProjS1.col(countAPS-1);
                    vAProjS.insert_cols(countAPS-1, tempV);
                    WDString_temp(j,jj) = 0;
                    WDString_temp(jj,j) = 0;
                }
            }
            
        }
        // rAProjS.print("rAProjS: ");
        // vAProjS.print("vAProjS: ");
        // cout << "44444444444444444444conA" << endl;
        
        // Check for the nearby strings
        mat uA_AProjS(2,1,fill::zeros);
        if(countAPS > 0)
        {
            mat XAProjS = join_cols(rAProjS, vAProjS);
            mat sig_para(4,1);
            sig_para(0,0) = A_A_D;
            sig_para(1,0) = B_A_D;
            sig_para(2,0) = C_A_D;
            sig_para(3,0) = D_A_D;
            vec potentialControl_result = potentialControl(1e-5, XA.col(i), XAProjS, rho_c_A, sig_para, R_m_AD, R_bar_AD, R_u_AD, Rij0(0), kADr, kADv, alphaADv);
            uA_AProjS = potentialControl_result.subvec(0,1);
        }
        // uA_AProjS.print("uA_AProjS: ");
        // cout << "555555555555555555555555conA" << endl;
        mat uA_AProjC(2,1,fill::zeros);
        int ci = clusteridA(i);
        // cout << "ci: " << ci << endl;
        // rAcm.print("rAcm");
        // rAcm.col(ci).print("rAcm col ci: ");
        // cout << "rhoA_con ci: " << rhoA_con(ci) << endl;
        if(arma::norm(rA-rAcm.col(ci)) < rhoA_con(ci))
        {
            // cout << "enter if loop" << endl;
            double thetaAAcm = atan2(rA(1) - rAcm(1, ci), rA(0) - rAcm(0, ci));
            mat rAprojC(2,1);
            mat tempM(2,1);
            // cout << "5.0" << endl;
            tempM(0) = cos(thetaAAcm);
            tempM(1) = sin(thetaAAcm);
            rAprojC = rAcm.col(ci) + rhoA_con(ci) * tempM;
            // cout << "5.1" << endl;
            mat rTP(2,1);
            rTP(0) = cos(thetaAAcm + M_PI/2);
            rTP(1) = sin(thetaAAcm + M_PI/2);
            tempM = rTP.t()*vA;
            
            if (tempM(0,0) < 0)
            {
                rTP = -rTP;
            }
            mat vAProjC = tempM(0,0) * rTP + vAcm.col(ci);
            mat XAProjC = join_cols(rAprojC, vAProjC);
            mat sig_para(4,1);
            sig_para(0,0) = A_A_A;
            sig_para(1,0) = B_A_A;
            sig_para(2,0) = C_A_A;
            sig_para(3,0) = D_A_A;
            vec potentialControl_result = potentialControl(1e-5, XA.col(i), XAProjC, rho_c_A, sig_para, R_m_AA, R_bar_AA, R_u_AA, Rik00(0), kAOr2, kAOv2, alphaAOv);
            uA_AProjC = potentialControl_result.subvec(0,1);
        }
        // uA_AProjS.print("uA_AProjC: ");
        // cout << "666666666666666666666conA" << endl;
        mat duA_goal;
        if(i == leaderIDA(i))
        {
            // cout << "if loop" << endl;
            duA_goal = XA_goal_dot.submat(2,i,3,i) - (kAPr * (rA-rA_goal) + kAPv * (vA-vA_goal));
        } else 
        {
            // cout << "else loop" << endl;
            // XA_goal_dot.submat(2,i,3,i).print("XA_goal_dot:");
            // cout << "kAPr2: " <<kAPr2  << endl;
            // cout << "kAPv2: " <<kAPv2  << endl;
            // rA.print("rA:");
            // rA_goal.print("rA_goal:");
            // vA.print("vA:");
            // vA_goal.print("vA_goal: ");
            // mat tempM;
            // tempM = (kAPr2 * (rA-rA_goal));
            // tempM.print(" kAPr2 * (rA-rA_goal)");
            // tempM = (kAPv2 * (vA-vA_goal));
            // tempM.print(" (kAPv2 * (vA-vA_goal)");
            // tempM = XA_goal_dot.submat(2,i,3,i) - (kAPr2 * (rA-rA_goal)) + (kAPv2 * (vA-vA_goal));
            // tempM.print("final result: ");
            duA_goal = XA_goal_dot.submat(2,i,3,i) - (kAPr2 * (rA-rA_goal) + kAPv2 * (vA-vA_goal));
        }
        // duA_goal.print("duA_goal");
        double norm_duA_goal = norm(duA_goal);
        // cout << "777777777777777777777conA" << endl;
        // cout << "norm_duA_goal: " << norm_duA_goal << endl;
        if(norm_duA_goal > 1e-10)
        {
            // cout << "u_maxA: " << u_maxA(i) << endl;
            // duA_goal.print("duA_goal");
            // uAFv.print("uAFv");
            // uAFr.print("uAFr");
            // uAOv.print("uAOv");
            // uAOr.print("uAOr");
            // uAD_pot.print("uAD_pot");
            // uA_AProjS.print("uA_AProjS");
            // cout << "C_d: " << C_d << endl;
            // cout << "norm(vA): " << arma::norm(vA) << endl;
            // vA.print("vA");
            // cout << "i: " << i << endl;
            uA.insert_cols(i, min(.7*u_maxA(i),norm_duA_goal) * (duA_goal/norm_duA_goal)+uAFv+uAFr+uAOv+uAOr+uAD_pot+uA_AProjS+uA_AProjC+C_d*arma::norm(vA)*vA);
        } else
        {
            uA.insert_cols(i, uAFv+uAFr+uAOv+uAOr+uAD_pot+uA_AProjS+uA_AProjC+C_d*arma::norm(vA)*vA);
        }
        // uA.print("uA: ");
        // cout << "88888888888888888888888conA" << endl;
        uA0.col(i) = uA.col(i);
        double uA_bar;
        
        for (int ca = 0; ca <= clusteridA.max(); ca++)
        {
            if(clusteridA(i) == ca)
            {
                // cout << "enter loop when ca = " << ca << endl;
                uA_bar = u_maxA(i);
                if (i == leaderIDA(i)){
                    uA_bar = 0.7*u_maxA(i);
                    double R_underbar = 3*R_bar_AD2;
                    double R_bar = 3*R_u_AD2;
                    double R_m;
                    if(flagAEnclosed(i)){
                        R_m = R_m_AD;
                    } else
                    {
                        R_m = 5*R_m_AD2;
                    }
                    double norm_uA = arma::norm(uA.col(i));
                    uAD_pot = zeros<vec>(2);
                    if(flagAEnclosed(i)){
                        vec potentialControl_result = potentialControl(1e-5, XA.col(i),XD, 4*rho_c_A, sigma_parameters(R_underbar, R_bar), R_m, R_underbar, R_bar, R_bar+10, kADr, kADv, alphaADv);
                    } else
                    {
                        for (int c = 0; c < assign.size(); c++)
                        {
                            vector<Point> points;
                            for (int i = 0; i < assign[c].n_elem; i++)
                                {
                                    Point tempP;
                                    tempP.x = XD(0,assign[c](i));
                                    tempP.y = XD(1,assign[c](i));
                                    points.push_back(tempP);
                                }
                            vector<Point> ans = convex_hull(points);
                            Point rDc_P = compute2DPolygonCentroid(ans);
                            mat rDc(2,1);
                            rDc(0,0) = rDc_P.x;
                            rDc(1,0) = rDc_P.y;
                            mat tempM(2, assign[c].n_elem);
                            for(int i=0; i<assign[c].n_elem; i++){
                                tempM.col(i) = XD.submat(2,assign[c](i),3,assign[c](i));
                            }
                            mat vDc = mean(tempM,1);
                            double Rad = 0;
                            for(int dj = 0; dj < ans.size(); dj++){
                                mat tempM1(2,1);
                                tempM1(0,0) = ans[dj].x;
                                tempM1(1,0) = ans[dj].y;
                                double dis = arma::norm(rDc - tempM1);
                                if (dis > Rad)
                                {
                                    Rad = dis;
                                }
                            }
                            R_m = Rad + R_m_AD2;
                            R_underbar = R_m + R_bar_AD2;
                            R_bar = R_m + R_bar_AD2 + 5;
                            tempM = join_cols(rDc,vDc);
                            vec potentialControl_result = potentialControl(1e-5, XA.col(i), tempM, 4*rho_c_A, sigma_parameters(R_underbar, R_bar), R_m, R_underbar, R_bar, R_bar+10, kADr, 0, alphaADv);
                            uAD_pot = uAD_pot + potentialControl_result.subvec(0,1);
                        }
                        

                    }
                    // uAD_pot.print("uAD_pot: ");
                    uvec i1 = find(clusteridA == ca-1); //ca-1???
                    if (!i1.is_empty()){
                        vec potentialControl_result = potentialControl(1e-5, XA.col(i), XA.col(leaderIDA(i1(0))),4*rho_c_A, sigma_parameters(R_underbar, R_bar), 1.5*R_m, R_underbar, R_bar, R_bar+10, kADr, 0, alphaADv);
                        uAD_pot = uAD_pot + potentialControl_result.subvec(0,1);
                    }

                    double norm_uAD_pot = arma::norm(uAD_pot);
                    uA.col(i) = uA.col(i) + uAD_pot;
                }
                double norm_uA = arma::norm(uA.col(i));
                if(norm_uA > uA_bar) {
                    uA.col(i) = uA.col(i) * uA_bar / norm_uA;
                }
                break;
            }
            

        }
        // cout << "9999999999999999999conA" << endl;
        
        vA_dot.insert_cols(i, uA.col(i)) ;
        vA_Des.insert_cols(i, zeros<mat>(2,1));
        vA_Des_dot.insert_cols(i,zeros<mat>(2,1));
        SigmaProdD(i,0) = sigmaProdD;
    }
    
    control_attacker_t control_result;
    control_result.uA = uA;
    control_result.uA0 = uA0;
    control_result.vA_des = vA_Des;
    control_result.vA_des_dot = vA_Des_dot;
    control_result.R_AO_min = R_AO_min;
    control_result.R_AAProjS_min = R_AAProjS_min;
    control_result.F_A = F_A;
    control_result.F_A_dot = F_A_dot;
    control_result.SigmaProdD = SigmaProdD;
    return control_result;
}

// int main(){
    
//     calAllparametersExperiment(18,18);
    
//     mat XA, XA_goal, XA_goal_dot, leaderIDA, XD, WA, R_tilde_AA, WDString, flagAEnclosed, rhoA_con;
//     vec assign1, NAinClusterA;
//     XA.load("/home/damon/Downloads/multi_swarm/controlA/XA.txt");
//     XA_goal.load("/home/damon/Downloads/multi_swarm/controlA/XA_goal.txt");
//     XA_goal_dot.load("/home/damon/Downloads/multi_swarm/controlA/XA_goal_dot.txt");
//     leaderIDA.load("/home/damon/Downloads/multi_swarm/controlA/leaderIdA.txt");
//     leaderIDA = leaderIDA -1;
//     leaderIDA.print("leaderIDA:");
//     XD.load("/home/damon/Downloads/multi_swarm/controlA/XD.txt");
//     WA.load("/home/damon/Downloads/multi_swarm/controlA/WA.txt");
//     R_tilde_AA.load("/home/damon/Downloads/multi_swarm/controlA/R_tilde_AA.txt");
//     WDString.load("/home/damon/Downloads/multi_swarm/controlA/WDString.txt");
//     clusteridA.load("/home/damon/Downloads/multi_swarm/controlA/clusterIdA.txt");
//     clusteridA = clusteridA - 1;
//     clusteridA.print("clusteridA");
//     indAinClusterA.print("indAinClusterA:");
//     // indAinClusterA.reset();
//     // indAinClusterA.set_size(3,1);
//     // vec tempV;
//     // tempV = {9,10,11,6,7,8};
//     // tempV = tempV-1;
//     // indAinClusterA(0) = tempV;
//     // tempV = {4,5,1,2,3};
//     // tempV = tempV-1;
//     // indAinClusterA(1) = tempV;
//     // tempV = {15,12,13,14,16,17,18};
//     // tempV = tempV-1;
//     // indAinClusterA(2) = tempV;
//     // indAinClusterA.print("indAinClusterA");
//     NAinClusterA.load("/home/damon/Downloads/multi_swarm/controlA/NAinClusterA.txt");
//     NAinClusterA = {18};
//     NAinClusterA.print("NAinClusterA:");
//     // NAinClusterA.print("NAinClusterA")
//     assign1.load("/home/damon/Downloads/multi_swarm/controlA/assign1.txt");
//     assign1 = assign1 - 1;
//     flagAEnclosed.load("/home/damon/Downloads/multi_swarm/controlA/flagAEnclosed.txt");
//     vector<vec> assign;
//     assign.push_back(assign1);
//     assign[0].print("assign1:");
//     rhoA_con.load("/home/damon/Downloads/multi_swarm/controlA/rhoA_con.txt");
//     // WDString.print("WDString:");
//     // cout << "rows: " << WDString.n_rows << " cols: " << WDString.n_cols << endl;
//     control_attacker_t a = controlAttacker4(XA, XA_goal, XA_goal_dot, leaderIDA, XD, WA, R_tilde_AA, WDString, 18, 18, clusteridA, indAinClusterA, NAinClusterA, rhoA_con, flagAEnclosed, assign, 0,100);
//     a.uA.print("uA:");
//     a.uA0.print("uA0:");
//     return 0;
// }

// int main() {
// //     double R_bar = 10.4177;
// //     double R_underbar = 7.4177;
// //     mat a = sigma_parameters(R_underbar, R_bar);
// //     a.print("a:");

// //     double alphav_12 = 0.5;
// //     double kr_12 = 0.5;
// //     double kv_12 = 0.5;
// //     double R_m_12 = 5.4177;
// //     double rho_sens_1 = 20;
// //     double R_tilde = 20.4176960858139;
// //     mat X1(4,1);
// //     X1(0,0) = 6.15480000000000;
// //     X1(1,0) = -33.0553000000000;
// //     X1(2,0) = 0;
// //     X1(3,0) = 0;
// //     mat X2(4,2);
// //     X2 = {{9.37240000000000,16.1838000000000},
// //          {-7.12780000000000,-5.84110000000000},
// //          {0,0},
// //          {0,0}};
// //     vec re = potentialControl(X1,X2, rho_sens_1, a.as_col(),R_m_12, R_underbar, R_bar,R_tilde, kr_12,kv_12,alphav_12);

// //     vec r = {6.15494381116162,
// //         -33.0547123910252};
// //     vec r1 = {9.37240000000000,
// // -7.12780000000000};
// //     vec r2 = {11.9907000000000,
// // -10.7580000000000};
// //     vec aa = projectionOnLine(r,r1,r2);
// //     aa.print("aa:");
//     calAllparametersExperiment();
//     mat XA;
//     mat XA_goal;
//     mat XA_goal_dot;
//     int flagEnclose;
//     int flagHerd;
//     mat BetaAv0;
//     mat SpeedA0;
//     mat mAv0;
//     mat cAv0;
//     mat XD;
//     mat WA;
//     mat WDString;
//     int ND = 3;
    
//     XA.load("../../../../../Downloads/swarm_matlab/controlA/XA.txt");
//     XA.print("XA:");
//     XA_goal.load("../../../../../Downloads/swarm_matlab/controlA/XA_goal.txt");
//     XA_goal.print("XA_goal:");
//     XA_goal_dot.load("../../../../../Downloads/swarm_matlab/controlA/XA_goal_dot.txt");
//     XA_goal_dot.print("XA_goal_dot:");
//     XD.load("../../../../../Downloads/swarm_matlab/controlA/XD.txt");
//     XD.print("XD:");
//     WA.load("../../../../../Downloads/swarm_matlab/controlA/WA.txt");
//     WA.print("WA:");
//     WDString.load("../../../../../Downloads/swarm_matlab/controlA/WDString.txt");
//     WDString.print("WDString:");
//     control_attacker_t a =  controlAttacker4(XA,XA_goal, XA_goal_dot, 0,0,XD,WA,WDString,NA,ND);
//     a.uA.print("uA:");
//     a.uA0.print("uA0:");

//     return 0;
// }