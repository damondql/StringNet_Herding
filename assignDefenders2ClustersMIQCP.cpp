#pragma once

#include "gurobi_c++.h"
#include <armadillo>
// #include <chrono>
using namespace std;
using namespace arma;
using namespace std::chrono;

struct new_assign_elem
{
    std::vector<vec>assign;
    std::vector<vec> indDinSwarm;
    double cost;
};



new_assign_elem assignDefenders2ClustersMIQCP(mat rAcm, mat rD, int ND, vec NAinClusterAD, vec assign0, vec indDinSwarm0, double R_DD_string){
    new_assign_elem assign_result;
    int NClusterAD = rAcm.row(0).n_elem;
    mat Dist(ND, NClusterAD);
    for (int i0 = 0; i0 < ND; i0++)
    {
        int i = assign0(i0);
        // cout << "i: " << i << endl;
        for (int c = 0; c < NClusterAD ; c++)
        {
            // cout << "c: " << c << endl;
            Dist(i0,c) = arma::norm(rD.col(i)- rAcm.col(c));
        }
        
    }
    // Dist.print("Dist:");
    mat Dist_col = Dist.as_col();
    int N = Dist_col.n_elem;

    int NV = ND * NClusterAD;
    mat A = zeros(ND+NClusterAD+1, NV);
    // cout << "A.col: " << A.n_cols << endl;
    mat b = zeros<mat>(ND+NClusterAD+1,1);
    // A(17,53) = 1;
    for (int j = 0; j < ND; j++)
    {
        vec tempV = regspace(j,ND,NV-1);
        // tempV.print("tempV);
        for (int i = 0; i < tempV.n_elem; i++)
        {
            // cout<<"j: " << j << endl;
            // cout << "tempV(i): " << tempV(i) << endl;
            A(j,tempV(i)) = 1;
        }
        b(j,0) = 1;
    }
    // cout << "finish first for loop" << endl;
    for (int c = 0; c < NClusterAD; c++)
    {
        vec tempV = regspace(ND*(c),1,ND*(c+1)-1);
        for (int i = 0; i < tempV.n_elem; i++)
        {
            A(ND+c, tempV(i)) = 1;
        }
        b(ND+c,0) = NAinClusterAD(c);
        
    }
    A.row(ND+NClusterAD) = ones<mat>(1,NV);
    b(ND+NClusterAD,0) = ND;
    // A.print("A: ");
    // b.print("b: ");
    
    // A.save("result_A.txt",raw_ascii);
    // b.save("result_b.txt",raw_ascii);
    mat result(NV,1);
    double cost;
    try {
        GRBEnv* env = 0;
        GRBVar* open = 0;


        // Model
        env = new GRBEnv();
        GRBModel model = GRBModel(*env);
        model.set(GRB_StringAttr_ModelName, "grbTest");

        // Plant open decision variables: open[p] == 1 if plant p is open.
        open = model.addVars(N, GRB_BINARY);

        int p;
        for (p = 0; p < N; ++p)
        {
            ostringstream vname;
            vname << "var" << p;
            open[p].set(GRB_StringAttr_VarName, vname.str());
        }
        
        // Set objective
        GRBQuadExpr obj = 0;
        // obj += Pbar(0,0)*a*a + Pbar(0,1)*a*b + Pbar(0,2)*a*c + Pbar(0,3)*a*d + Pbar(0,4)*a*e + Pbar(0,5)*a*f + Pbar(0,6)*a*g + Pbar(0,7)*a*h + Pbar(0,8)*
        // for (int i = 0; i < Pbar.n_rows; i++)
        // {
        //     for (int j = 0; j < Pbar.n_cols; j++)
        //     {
        //         obj += Pbar(i,j)* open[i] * open[j];
        //     }
        // }

        for (int i = 0; i < Dist_col.n_rows; i++)
        {
            obj += Dist_col(i,0) * open[i];
        }

        // cout << "finish obj" << endl;

        model.setObjective(obj);

        for (int p = 0; p < ND+NClusterAD+1; p++)
        {
            GRBLinExpr ptot = 0;
            for (int w = 0; w < NV; ++w)
            {
            ptot += A(p,w)*open[w];
            }
            ostringstream cname;
            cname << "Constriain" << p;
            model.addConstr(ptot == b(p), cname.str());
        }


        for (int c = 0; c < NClusterAD; c++)
        {
            mat tempM = zeros<mat>(NV,NV);
            tempM.submat(ND*c, ND*c+1 ,ND*(c+1)-2,ND*(c+1)-1) = -arma::eye(ND-1, ND-1);
            GRBQuadExpr QExp = 0;
            for (int j = 0; j < NV; j++)
            {
                for (int k = 0; k < NV; k++)
                {
                    QExp += tempM(j,k) * open[j] * open[k];
                }
            }
            ostringstream cname;
            cname << "QConstriain" << c;
            model.addQConstr(QExp <= -NAinClusterAD(c)+1, cname.str());
        }
        
        // cout << "finished adding constrains" << endl;
        model.optimize();

        
        cost = model.get(GRB_DoubleAttr_ObjVal);



        for (int i = 0; i < NV; i++)
        {
            result(i) = open[i].get(GRB_DoubleAttr_X);
        }
        
        


        // cal_result.result.resize(ND*ND,1);
        // for (int i = 0; i < ND*ND; i++)
        // {
        //     cal_result.result(i) = open[i].get(GRB_DoubleAttr_X);
        // }

        // cal_result.obj = model.get(GRB_DoubleAttr_ObjVal);

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    // result.print("result: ");
    // Dist_col.print("obj:");
    mat cost_mat = Dist_col.t() * result;
    // double cost = cost_mat(0);
    // cout << "cost: " << cost << endl;
    assign_result.cost = cost;

    uvec tempInd = find(result >= 0.9999);
    int tempSum = 0;
    // tempInd.print("tempInd:");
    for (int c = 0; c < NClusterAD; c++)
    {
        if(c > 0)
        {
            tempSum += NAinClusterAD(c-1);
        }
        uvec tempV;
        // cout  << "tempSum " << tempSum << endl;
        // cout << "tempSum+NAinClusterAD(c): " << tempSum+NAinClusterAD(c)-1 << endl;
        // tempInd.subvec(tempSum, tempSum+NAinClusterAD(c)-1).print("tempIND");
        tempV = tempInd.subvec(tempSum, tempSum+NAinClusterAD(c)-1);
        // tempInd.subvec(tempSum,tempSum+NAinClusterAD(c)-1);
        tempV = tempV - ND*c;
        // tempV.print("tempV: ");
        vec tempV2(tempV.n_elem);
        for (int i = 0; i < tempV2.n_elem; i++)
        {
            tempV2(i) = assign0(tempV(i));
        }
        tempV2.print("assign: ");
        assign_result.assign.push_back(tempV2);
        vec tempV3(tempV.n_elem);
        for (int i = 0; i < tempV2.n_elem; i++)
        {
            tempV3(i) = indDinSwarm0(tempV(i));
        }
        tempV3.print("indDinSwarm");
        assign_result.indDinSwarm.push_back(tempV3);
        
    }
    
    return assign_result;
    
}

// int main(){
//     mat rAcm;
//     mat rD;
//     mat NAinClusterAD_mat; 
//     mat assign0_mat; 
//     mat indDinSwarm0_mat;
//     vec NAinClusterAD; 
//     vec assign0; 
//     vec indDinSwarm0;
//     rAcm.load("/home/damon/Downloads/multi_swarm/new_assign/rAcm.txt");
//     rD.load("/home/damon/Downloads/multi_swarm/new_assign/rD.txt");
//     NAinClusterAD_mat.load("/home/damon/Downloads/multi_swarm/new_assign/NAinClusterAD.txt");
//     assign0_mat.load("/home/damon/Downloads/multi_swarm/new_assign/assign0.txt");
//     indDinSwarm0_mat.load("/home/damon/Downloads/multi_swarm/new_assign/indDinSwarm0.txt");
//     NAinClusterAD = NAinClusterAD_mat.as_col();
//     assign0 = assign0_mat.as_col();
//     indDinSwarm0 = indDinSwarm0_mat.as_col();
//     double R_DD_String = 80;
//     int ND = 18;
//     assign0 = assign0 -1;
//     indDinSwarm0 = indDinSwarm0 - 1;
//     assign0.print("assign0: ");
//     indDinSwarm0.print("indDinSwarm0:");
//     // NAinClusterAD.print("NAinClusterAD");
//     // assign0.print("assign:");
//     // indDinSwarm0.print("indinSwarm");
//     new_assign_elem a = assignDefenders2ClustersMIQCP(rAcm, rD, ND, NAinClusterAD, assign0, indDinSwarm0, R_DD_String);
//     cout << "cost: " << a.cost <<endl;
//     for (int i = 0; i < a.assign.size(); i++)
//     {
//         a.assign[i].print("assign");
//         a.indDinSwarm[i].print("indDinSwarm:");
//     }
    
// }