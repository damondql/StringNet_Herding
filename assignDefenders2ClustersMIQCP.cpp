#pragma once

#include "gurobi_c++.h"
#include <armadillo>
// #include <chrono>
using namespace std;
using namespace arma;
using namespace std::chrono;

void assignDefenders2ClustersMIQCP(mat rAcm, mat rD, int ND, vec NAinClusterAD, vec assign0, vec indDinSwarm0, double R_DD_string){
    int NClusterAD = rAcm.row(0).n_elem;
    mat Dist(ND, NClusterAD);
    for (int i0 = 0; i < ND; i++)
    {
        int i = assign0(i0);
        for (int c = 1; c < NClusterAD ; c++)
        {
            Dist(i0,c) = arma::norm(rD.col(i)- rAcm.col(c));
        }
        
    }

    mat Dist_col = Dist.as_col();
    int N = Dist_col.n_elem;

    int NV = ND * NClusterAD;
    mat A = zeros(NE+NClusterAD+1, NV);

    for (int j = 0; j < ND; i++)
    {
        vec tempV = regspace(j,ND,NV);
        for (int i = 0; i < count; i++)
        {
            A(j,tempV(i)) = 1;
        }
        
    }
    mat b = ones<mat>(ND,1);
    
    // try {
    //     GRBEnv* env = 0;
    //     GRBVar* open = 0;


    //     // Model
    //     env = new GRBEnv();
    //     GRBModel model = GRBModel(*env);
    //     model.set(GRB_StringAttr_ModelName, "grbTest");

    //     // Plant open decision variables: open[p] == 1 if plant p is open.
    //     open = model.addVars(N, GRB_BINARY);

    //     int p;
    //     for (p = 0; p < N; ++p)
    //     {
    //         ostringstream vname;
    //         vname << "var" << p;
    //         open[p].set(GRB_StringAttr_VarName, vname.str());
    //     }
        
    //     // Set objective
    //     GRBQuadExpr obj = 0;
    //     // obj += Pbar(0,0)*a*a + Pbar(0,1)*a*b + Pbar(0,2)*a*c + Pbar(0,3)*a*d + Pbar(0,4)*a*e + Pbar(0,5)*a*f + Pbar(0,6)*a*g + Pbar(0,7)*a*h + Pbar(0,8)*
    //     // for (int i = 0; i < Pbar.n_rows; i++)
    //     // {
    //     //     for (int j = 0; j < Pbar.n_cols; j++)
    //     //     {
    //     //         obj += Pbar(i,j)* open[i] * open[j];
    //     //     }
    //     // }

    //     for (int i = 0; i < Dist_col.n_rows; i++)
    //     {
    //         obj += Dist_col(i,0) * open[i];
    //     }

    //     // cout << "finish obj" << endl;

    //     model.setObjective(obj);
    //     // cout << "start setting up matrix A" << endl;
    //     mat A(2*ND,ND*ND,fill::zeros);
    //     for (int i = 1; i <= ND; i++)
    //     {
    //         A.submat(i-1, ND*(i-1), i-1, ND*i-1) = ones<mat>(1,ND);
    //         vec tempV = regspace(i, ND, ND*ND);
    //         for (int j = 0; j < tempV.n_rows; j++)
    //         {
    //             A(ND+i-1, tempV(j)-1) = 1;
    //         }
            
    //     }

    //     // cout << "finished setting up matrix A" << endl;
    //     for (int p = 0; p < ND*2; p++)
    //     {
    //         GRBLinExpr ptot = 0;
    //         for (int w = 0; w < ND*ND; ++w)
    //         {
    //         ptot += A(p,w)*open[w];
    //         }
    //         ostringstream cname;
    //         cname << "Constriain" << p;
    //         model.addConstr(ptot == 1, cname.str());
    //     }
    //     // cout << "finished adding constrains" << endl;
    //     model.optimize();


    //     cal_result.result.resize(ND*ND,1);
    //     for (int i = 0; i < ND*ND; i++)
    //     {
    //         cal_result.result(i) = open[i].get(GRB_DoubleAttr_X);
    //     }

    //     cal_result.obj = model.get(GRB_DoubleAttr_ObjVal);

    // } catch(GRBException e) {
    //     cout << "Error code = " << e.getErrorCode() << endl;
    //     cout << e.getMessage() << endl;
    // } catch(...) {
    //     cout << "Exception during optimization" << endl;
    // }


}

int main(){
    mat rAcm;
    mat rD;
    int ND;
    vec NAinClusterAD; 
    vec assign0; 
    vec indDinSwarm0;
}