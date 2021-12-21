#pragma once
#include "AllParameters.cpp"
#include "helperFunction.cpp"
#include <armadillo>

using namespace arma;

mat controlDefenderFormation(mat XD, mat indD, mat XD_des, mat XD_des_dot, mat rhoA_con, mat XA, std::vector<uvec> indAinClusterAD, int NClusterAD, int flagAvoidAcon, int NA, int ND, mat& uD){
    // int nid = indDef.n_elem;
    int nd = indD.n_elem;
    // int na = XA.n_cols;
    // mat uD(2,nd,fill::zeros);
    mat uD0(2, nd, fill::zeros);
    mat tempM;
    // tempM = XD;
    // for (int i = 0; i < nd; i++)
    // {
    //     XD.col(i) = tempM.col(indD(i)-1);
    // }
    // mat rAcm = arma::sum(XA.submat(0,0,1,XA.n_cols-1),1)/NA;
    // mat vAcm = arma::sum(XA.submat(2,0,3,XA.n_cols-1),1)/NA;
    mat rAcm;
    mat vAcm;
    vec uD_Ac;
    for (int j = 0; j < nd; j++)
    {
        // cout << "j=========== " <<  j << endl;
        mat rD = XD.submat(0, indD(j), 1, indD(j));
        mat vD = XD.submat(2, indD(j), 3, indD(j));
        mat rD_des = XD_des.submat(0, j, 1, j);
        mat vD_des = XD_des.submat(2, j, 3, j);
        // rD.print("rD");
        // vD.print("vD");
        // rD_des.print("rD_des");
        // vD_des.print("vD_des");
        // cout << "kDFr2: " << kDFr2 << endl;
        // cout << "alphaDFr: " << alphaDFr << endl;
        // Avoid the boundary of the connectivity region of the attackers
        mat dr, dv;
        if(arma::norm(rD - rD_des) > 1e-5)
        {
            dr=-kDFr2*(rD-rD_des)*pow(arma::norm(rD-rD_des),(alphaDFr-1));
        } else
        {
            dr = zeros<mat>(2,1);
        }

        if(arma::norm(vD - vD_des) > 1e-5)
        {
            dv = -kDFv*(vD-vD_des)*pow(arma::norm(vD-vD_des),(alphaDFv-1));
        } else
        {
            dv = zeros<mat>(2,1);
        }

        mat XD_ind0_input;
        for(int i = 0 ; i < ND; i++)
        {
            if(i != indD(j))
            {
                XD_ind0_input.insert_cols(XD_ind0_input.n_cols, XD.col(i));
            }
        }

        mat sigma_para(1,4);
        sigma_para(0) = A_D_D;
        sigma_para(1) = B_D_D;
        sigma_para(2) = C_D_D;
        sigma_para(3) = D_D_D;
        vec uD_D_potential = potentialControl(1e-5, XD.col(indD(j)), XD_ind0_input, R_u_DD, sigma_para, R_m_DD, R_bar_DD, R_u_DD, Rjj0(j), kDDr, kDDv, alphaDDv);
        vec uD_D = uD_D_potential.subvec(0,1);
        if(flagAvoidAcon)
        {
            mat XDProj(4, NClusterAD);
            for (int ca = 0; ca < NClusterAD; ca++)
            {
                mat XA_input;
                for (int n_ca = 0; n_ca < indAinClusterAD[ca].n_elem; n_ca++)
                {
                    XA_input.insert_cols(XA_input.n_cols, XA.col(indAinClusterAD[ca](n_ca)));
                }
                rAcm = arma::sum(XA_input.submat(0,0,1,XA_input.n_cols-1),1) / indAinClusterAD[ca].n_elem;
                vAcm = arma::sum(XA_input.submat(2,0,3,XA_input.n_cols-1),1) / indAinClusterAD[ca].n_elem;
                // rAcm.print("rAcm: ");
                // vAcm.print("vAcm: ");
                double thetaDAcm = atan2(rD(1) - rAcm(1), rD(0) - rAcm(0));
                tempM.resize(2,1);
                tempM(0) = cos(thetaDAcm);
                tempM(1) = sin(thetaDAcm);
                mat rDProj = rAcm + rhoA_con(ca) * tempM;
                mat rTP(2,1);
                rTP(0) = cos(thetaDAcm + M_PI/2);
                rTP(1) = sin(thetaDAcm + M_PI/2);
                if(j < nd/2)
                {
                    rTP = -rTP;
                }
                tempM = rTP.t() * vD;
                // tempM.print("tempM: ");
                mat vDProj = tempM(0,0)*rTP + vAcm;
                mat tempV;
                tempV = join_cols(rDProj, vDProj);
                XDProj.col(ca) = tempV;
            }
            mat sig_para(1,4);
            sig_para = A_A_D;
            sig_para = B_A_D;
            sig_para = C_A_D;
            sig_para = D_A_D;
            vec uD_Ac_potential = potentialControl(1e-5, XD.col(indD(j)), XDProj, R_u_AD ,sig_para, R_m_DD, R_bar_AD, R_u_DD, Rjk00(0), kDDr, kDDv, alphaDDv);
            uD_Ac = uD_Ac_potential.subvec(0,1);
            
        }
        // cout <<"setting uD0 col j:" << j << endl;
        // uD_D.print("uD_D: ");
        // uD_Ac.print("uD_Ac");
        // dr.print("dr:");
        // dv.print("dv:");
        uD0.col(j) = XD_des_dot.submat(2, j, 3, j) + dr + dv + uD_D + uD_Ac;
        // cout << "finish seeting uD0 col j" << endl;
    }
    // cout << "saturation" << endl;
    for(int j = 0; j < nd; j++)
    {
        double norm_uD = arma::norm(uD0.col(j));
        if(norm_uD > u_maxD(j))
        {
            uD0.col(j) = uD0.col(j) * u_maxD(j) / norm_uD;
        }
    }
    // cout << "finish saturation" << endl;
    // tempM.reset();
    // tempM = uD;
    // for (int i = 0; i < uD.n_cols; i++)
    // {
    //     uD.col(indD(i)-1) = tempM.col(i);
    // }

    for(int i = 0; i < uD0.n_cols; i++)
    {
        uD.col(indD(i)) = uD0.col(i);
    }


    return uD;


    // for (int j = 0; j < nd; j++)
    // {
    //     mat rD = XD.submat(0,j,1,j);
    //     mat rD_des = XD_des.submat(0,j,1,j);
    //     mat vD = XD.submat(2,j,3,j);
    //     mat vD_des = XD_des.submat(2,j,3,j);
    //     mat dr, dv;
    //     if(arma::norm(rD - rD_des) > 1e-5)
    //     {
    //         dr=-kDFr2*(rD-rD_des)*pow(arma::norm(rD-rD_des),(alphaDFr-1));
    //     } else
    //     {
    //         dr = zeros<mat>(2,1);
    //     }

    //     if(arma::norm(vD - vD_des) > 1e-5)
    //     {
    //         dv = -kDFv*(vD-vD_des)*pow(arma::norm(vD-vD_des),(alphaDFv-1));
    //     } else
    //     {
    //         dv = zeros<mat>(2,1);
    //     }

    //     // check for nearby defenders
    //     std::vector<int> ind0;
    //     for(int i = 0; i < XD.n_cols; i++)
    //     {
    //         if (i != indD(j))
    //         {
    //             ind0.push_back(i);
    //         }
    //     }
    //     mat XD_sub;
    //     for(int i = 0; i < ind0.size(); i++){
    //         XD_sub.insert_cols(XD_sub.n_cols, XD.col(ind0[i]));
    //     }
    //     mat sigma_para(1,4);
    //     sigma_para(0, 0) = A_D_D;
    //     sigma_para(0, 1) = B_D_D;
    //     sigma_para(0, 2) = C_D_D;
    //     sigma_para(0, 3) = D_D_D;
    //     vec uD_D = potentialControl(1e-5 ,XD.col(indD(j)), XD_sub, R_u_DD, sigma_para, R_m_DD, R_bar_DD, R_u_DD, Rjj0(j), kDDr, kDDv, alphaDDv);

    //     mat uD_Ac(2, 1, fill::zeros);
    //     mat rAcm(2, 1, fill::zeros);
    //     mat vAcm(2, 1, fill::zeros);
    //     if (flagAvoidAcon)
    //     {
    //         for(int ca = 0; ca < NClusterAD; ca++)
    //         {
    //             for (int ca_i = 0; ca_i < indAinClusterAD[ca].n_elem; ca_i++)
    //             {
    //                 rAcm = rAcm + XA.submat(0,indAinClusterAD[ca](ca_i),1, indAinClusterAD[ca](ca_i));
    //                 vAcm = vAcm + XA.submat(2,indAinClusterAD[ca](ca_i),3, indAinClusterAD[ca](ca_i));
    //             }
    //             rAcm = rAcm / indAinClusterAD[ca].n_elem;
    //             vAcm = vAcm / indAinClusterAD[ca].n_elem;
    //             double thetaDAcm = atan2(rD(1)-rAcm(1),rD(0)-rAcm(0));
    //             tempM.resize(2,1);
    //             tempM(0,0) = cos(thetaDAcm);
    //             tempM(1,0) = sin(thetaDAcm);
    //             mat rDProj = rAcm + rhoA_con(ca) * tempM;
    //             mat rTP(2,1);
    //             rTP(0,0) = cos(thetaDAcm+M_PI/2);
    //             rTP(1,0) = sin(thetaDAcm+M_PI/2);
    //             tempM.reset();
    //             tempM = rTP.t() * vD;
    //             if (j < nd/2 )
    //             {
    //                 rTP = -rTP;
    //             }
    //             tempM.reset();
    //             tempM = rTP.t()*vD;
    //             mat vDProj = tempM(0,0) * rTP + vAcm;
    //             double R_underbar = 4;
    //             double R_bar = 3;
    //             double R_m = 1.5;
    //             mat XDproj;
    //             XDproj = arma::join_cols(rDProj, vDProj);
    //             // rvDproj.print("rvDproj");
    //             vec potentialControl_result = potentialControl(1e-5,XD.col(indD(j)),XDproj,R_u_DD,sigma_parameters(R_underbar,R_bar),R_m,R_underbar,R_bar, R_bar+10,kDOr,kDOv2, alphaDOv);
    //             uD_Ac = potentialControl_result.subvec(0,1);
    //         }
            
            
    //     }

    //     mat uDOv(2,1,fill::zeros);
    //     mat uDOr(2,1,fill::zeros);
    //     // mat uD1,uD2;
    //     // double norm_uD1, norm_uD2;
    //     // uD1=-kDDesr*(rD-XD_des.submat(0,j,1,j))-kDDesv*(vD-XD_des.submat(2,j,3,j)) + uAcon+uDOv+uDOr;
    //     // norm_uD1=arma::norm(uD1);
    //     // // cout << "norm_uD1: " << norm_uD1 << endl;
    //     // uD2=-kDDesv*(vD-XD_des.submat(2,j,3,j))+C_d*vD*arma::norm(vD);
    //     // norm_uD2=arma::norm(uD2);
    //     // // cout << "norm_uD2: " << norm_uD2 << endl;
    //     // uD.col(j)=XD_des_dot.submat(2,j,3,j) + std::min(umd_e1,norm_uD1)*uD1/norm_uD1 + std::min(umd_e2,norm_uD2)*uD2/norm_uD2;
    //     uD0.col(j) = XD_des_dot.submat(2,j,3,j) + dr + dv + uD_D + uD_Ac + uDOv + uDOr;
    //     double norm_uD0 = arma::norm(uD0.col(j));
    //     if (norm_uD0 > u_maxD(j))
    //     {
    //         uD0.col(j) = uD0.col(j) * u_maxD(j) / norm_uD0;
    //     }
        
    // }
    // // uD.col(ND) = uDFc_trans;
    // // uD.print("uD after all cal before switch:");
    // tempM.reset();
    // tempM = uD;
    // for (int i = 0; i < uD.n_cols; i++)
    // {
    //     uD.col(indD(i)-1) = tempM.col(i);
    // }

    // return uD;
}

// int deleted_index = indD(j);
        // vec ind1 = linspace(0,deleted_index-1, deleted_index);
        // vec ind2 = linspace(deleted_index+1, NA-1, NA - deleted_index-1);
        // vec ind0 = join_cols(ind1, ind2);
        // mat XD_input;
        // for (int ind_i = 0; ind_i < ind0.n_elem; ind_i++)
        // {
        //     XD_input.insert_cols(XD_input.n_cols, XD.col(ind0(ind_i)));
        // }
        // mat sigma_para(1,4);
        // sigma_para(0) = A_D_D;
        // sigma_para(1) = B_D_D;
        // sigma_para(2) = C_D_D;
        // sigma_para(3) = D_D_D;
        // mat uD_D = potentialControl(1e-5, XD.col(indD(j)), XD_input, R_u_DD, sigma_para, R_m_DD, R_bar_DD, R_u_DD, Rjj0(j), kDDr, kDDv, alphaDDv);
        // uD0.col(j) = XD_des_dot.submat(2,j,3,j) - kDDesr * (rD - XD_des.submat(0,j,1,j)) - kDDesv * (vD - XD_des.submat(2,j,3,j)) + uD_D + uD_Ac
        //              + C_d * arma::norm(vD) * vD;

// int main(){
    
//     calAllparametersExperiment(18,18);
//     mat XD, indD, XD_des, XD_des_dot, rhoA_con, uD, YA;
//     int NClusterAD_ti = 3;
//     std::vector<uvec> indAinClusterAD(3);
//     uvec tempV = {0,1,2,3,4};
//     indAinClusterAD[0] = tempV;
//     tempV = {5,6,7,8,9,10};
//     indAinClusterAD[1] = tempV;
//     tempV = {11,12,13,14,15,16,17};
//     indAinClusterAD[2] = tempV;
//     XD.load("/home/damon/Downloads/multi_swarm/Formation/XD.txt");
//     indD.load("/home/damon/Downloads/multi_swarm/Formation/indD.txt");
//     indD = indD - 1;
//     XD_des.load("/home/damon/Downloads/multi_swarm/Formation/XD_des.txt");
//     XD_des_dot.load("/home/damon/Downloads/multi_swarm/Formation/XD_des_dot.txt");
//     rhoA_con.load("/home/damon/Downloads/multi_swarm/Formation/rhoA_con.txt");
//     uD.load("/home/damon/Downloads/multi_swarm/Formation/uD.txt");
//     YA.load("/home/damon/Downloads/multi_swarm/Formation/XA.txt");

//     XD.print("XD: ");
//     indD.print("indD: ");
//     XD_des.print("XD_des: ");
//     XD_des_dot.print("XD_des_dot: ");
//     rhoA_con.print("rhoA_con: ");
//     uD.print("uD: ");
//     YA.print("YA: ");


//     uD = controlDefenderFormation(XD, indD, XD_des, XD_des_dot, rhoA_con, YA, indAinClusterAD, NClusterAD_ti, 1, 18, 18, uD);
//     uD.print("uD: ");
//     return 0;
// }