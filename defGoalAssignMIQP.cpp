#pragma once

#include "gurobi_c++.h"
#include <armadillo>
// #include <chrono>
using namespace std;
using namespace arma;
using namespace std::chrono;

struct goal_assign
{
     vec assign;
     double cost;
     double maxOptT;
};


struct gurobi_result
{
     mat result;
     double obj;
};

goal_assign defGoalAssignMIQP(mat optT, mat Pbar, int ND)
{
     goal_assign result;
     gurobi_result cal_result;
  try {
     GRBEnv env = GRBEnv();
     env.set("NonConvex","2");
     GRBModel model = GRBModel(env);

     // Create variables

     GRBVar a = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "a");
     GRBVar b = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "b");
     GRBVar c = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "c");
     GRBVar d = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "d");
     GRBVar e = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "e");
     GRBVar f = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "f");
     GRBVar g = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "g");
     GRBVar h = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "h");
     GRBVar i = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "i");

    //  cout << "finish variable" <<endl;
     std::vector<GRBVar> variables;
     
     variables = {a,b,c,d,e,f,g,h,i};
    //  cout <<"make varibles vector" << endl;
     // Set objective
     GRBQuadExpr obj = 0;
     // obj += Pbar(0,0)*a*a + Pbar(0,1)*a*b + Pbar(0,2)*a*c + Pbar(0,3)*a*d + Pbar(0,4)*a*e + Pbar(0,5)*a*f + Pbar(0,6)*a*g + Pbar(0,7)*a*h + Pbar(0,8)*
     for (int i = 0; i < Pbar.n_rows; i++)
     {
          for (int j = 0; j < Pbar.n_cols; j++)
          {
               obj += Pbar(i,j)* variables[i] * variables[j];
          }
     }
     
     for (int i = 0; i < optT.n_rows; i++)
     {
          obj += optT(i,0) * variables[i];
     }
    //  cout << "finish obj" << endl;


     model.setObjective(obj);

     // Add constraint: x + 2 y + 3 z >= 4
     mat A(2*ND,ND*ND,fill::zeros);
     for (int i = 1; i <= ND; i++)
     {
          A.submat(i-1, ND*(i-1), i-1, ND*i-1) = ones<mat>(1,ND);
          vec tempV = regspace(i, ND, ND*ND);
          for (int j = 0; j < tempV.n_rows; j++)
          {
               A(ND+i-1, tempV(j)-1) = 1;
          }
          
     }


     model.addConstr(A(0,0)*a + A(0,1)*b + A(0,2)*c + A(0,3)*d + A(0,4)*e + A(0,5)*f + A(0,6)*g + A(0,7)*h + A(0,8)*i == 1, "c0");
     model.addConstr(A(1,0)*a + A(1,1)*b + A(1,2)*c + A(1,3)*d + A(1,4)*e + A(1,5)*f + A(1,6)*g + A(1,7)*h + A(1,8)*i == 1, "c1");
     model.addConstr(A(2,0)*a + A(2,1)*b + A(2,2)*c + A(2,3)*d + A(2,4)*e + A(2,5)*f + A(2,6)*g + A(2,7)*h + A(2,8)*i == 1, "c2");
     model.addConstr(A(3,0)*a + A(3,1)*b + A(3,2)*c + A(3,3)*d + A(3,4)*e + A(3,5)*f + A(3,6)*g + A(3,7)*h + A(3,8)*i == 1, "c3");
     model.addConstr(A(4,0)*a + A(4,1)*b + A(4,2)*c + A(4,3)*d + A(4,4)*e + A(4,5)*f + A(4,6)*g + A(4,7)*h + A(4,8)*i == 1, "c4");
     model.addConstr(A(5,0)*a + A(5,1)*b + A(5,2)*c + A(5,3)*d + A(5,4)*e + A(5,5)*f + A(5,6)*g + A(5,7)*h + A(5,8)*i == 1, "c5");

     // model.addConstr(a+b+c == 1, "c0");
     // model.addConstr(d+e+f == 1, "c1");
     // model.addConstr(g+h+i == 1, "c2");
     // model.addConstr(a+d+g == 1, "c3");
     // model.addConstr(b+e+h == 1, "c4");
     // model.addConstr(c+f+i == 1, "c5");

     // Add constraint: x + y >= 1


     // Optimize model

     model.optimize();

     // cout << a.get(GRB_StringAttr_VarName) << " "
     //      << a.get(GRB_DoubleAttr_X) << endl;
     // cout << b.get(GRB_StringAttr_VarName) << " "
     //      << b.get(GRB_DoubleAttr_X) << endl;
     // cout << c.get(GRB_StringAttr_VarName) << " "
     //      << c.get(GRB_DoubleAttr_X) << endl;
     // cout << d.get(GRB_StringAttr_VarName) << " "
     //      << d.get(GRB_DoubleAttr_X) << endl;
     // cout << e.get(GRB_StringAttr_VarName) << " "
     //      << e.get(GRB_DoubleAttr_X) << endl;
     // cout << f.get(GRB_StringAttr_VarName) << " "
     //      << f.get(GRB_DoubleAttr_X) << endl;
     // cout << g.get(GRB_StringAttr_VarName) << " "
     //      << g.get(GRB_DoubleAttr_X) << endl;
     // cout << h.get(GRB_StringAttr_VarName) << " "
     //      << h.get(GRB_DoubleAttr_X) << endl;
     // cout << i.get(GRB_StringAttr_VarName) << " "
     //      << i.get(GRB_DoubleAttr_X) << endl;

     // cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

     double a1 = a.get(GRB_DoubleAttr_X);
     double b1 = b.get(GRB_DoubleAttr_X);
     double c1 = c.get(GRB_DoubleAttr_X);
     double d1 = d.get(GRB_DoubleAttr_X);
     double e1 = e.get(GRB_DoubleAttr_X);
     double f1 = f.get(GRB_DoubleAttr_X);
     double g1 = g.get(GRB_DoubleAttr_X);
     double h1 = h.get(GRB_DoubleAttr_X);
     double i1 = i.get(GRB_DoubleAttr_X);

     cal_result.result = {a1,b1,c1,d1,e1,f1,g1,h1,i1};
     cal_result.result = cal_result.result.t();

     cal_result.obj = model.get(GRB_DoubleAttr_ObjVal);
     // cost.print("cost: ");
     // Change variable types to integer

     // x.set(GRB_CharAttr_VType, GRB_INTEGER);
     // y.set(GRB_CharAttr_VType, GRB_INTEGER);
     // z.set(GRB_CharAttr_VType, GRB_INTEGER);

     // // Optimize model

     // model.optimize();

     // cout << x.get(GRB_StringAttr_VarName) << " "
     //      << x.get(GRB_DoubleAttr_X) << endl;
     // cout << y.get(GRB_StringAttr_VarName) << " "
     //      << y.get(GRB_DoubleAttr_X) << endl;
     // cout << z.get(GRB_StringAttr_VarName) << " "
     //      << z.get(GRB_DoubleAttr_X) << endl;

     // cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

     } catch(GRBException e) {
     cout << "Error code = " << e.getErrorCode() << endl;
     cout << e.getMessage() << endl;
     } catch(...) {
     cout << "Exception during optimization" << endl;
     }
     // cal_result.result.print("get result: ");
     // cout << "obj: " << cal_result.obj << endl;
     mat cost_mat;
     cost_mat = optT.t()*cal_result.result + cal_result.result.t()*Pbar*cal_result.result;
     result.cost = cost_mat(0,0);
     // cost.print("cost: ");
     mat cal_OptT;
     cal_OptT.copy_size(cal_result.result);
     for (int i = 0; i < cal_OptT.n_rows; i++)
     {
          cal_OptT(i,0) = optT(i,0) * cal_result.result(i,0);
     }
     // cal_OptT.print("cal_OptT:");
     result.maxOptT = cal_OptT.max();
     // cout << "maxOptT: " << maxOptT << endl;
     uvec tempInd = find(cal_result.result == 1);
     
     tempInd += 1;
     
     result.assign.copy_size(tempInd);
     for (int i = 0; i < result.assign.n_rows; i++)
     {
          result.assign(i) = tempInd(i)%ND;
          if (result.assign(i) == 0)
          {
               result.assign(i) = ND;
          }
          
     }
     // assign.print("assign: ");

     return result;
}

goal_assign defGoalAssignMIQP_new(mat optT, mat Pbar, int ND)
{
     goal_assign result;
     gurobi_result cal_result;
  try {
     GRBEnv* env = 0;
     GRBVar* open = 0;

     // Model
     env = new GRBEnv();
     GRBModel model = GRBModel(*env);
     model.set(GRB_StringAttr_ModelName, "grbTest");

     // Plant open decision variables: open[p] == 1 if plant p is open.
     open = model.addVars(ND*ND, GRB_BINARY);

     int p;
     for (p = 0; p < ND*ND; ++p)
     {
          ostringstream vname;
          vname << "var" << p;
          open[p].set(GRB_StringAttr_VarName, vname.str());
     }
     
     // Set objective
     GRBQuadExpr obj = 0;
     // obj += Pbar(0,0)*a*a + Pbar(0,1)*a*b + Pbar(0,2)*a*c + Pbar(0,3)*a*d + Pbar(0,4)*a*e + Pbar(0,5)*a*f + Pbar(0,6)*a*g + Pbar(0,7)*a*h + Pbar(0,8)*
     for (int i = 0; i < Pbar.n_rows; i++)
     {
          for (int j = 0; j < Pbar.n_cols; j++)
          {
               obj += Pbar(i,j)* open[i] * open[j];
          }
     }

     for (int i = 0; i < optT.n_rows; i++)
     {
          obj += optT(i,0) * open[i];
     }

     cout << "finish obj" << endl;

     model.setObjective(obj);
     cout << "start setting up matrix A" << endl;
     mat A(2*ND,ND*ND,fill::zeros);
     for (int i = 1; i <= ND; i++)
     {
          A.submat(i-1, ND*(i-1), i-1, ND*i-1) = ones<mat>(1,ND);
          vec tempV = regspace(i, ND, ND*ND);
          for (int j = 0; j < tempV.n_rows; j++)
          {
               A(ND+i-1, tempV(j)-1) = 1;
          }
          
     }

     cout << "finished setting up matrix A" << endl;
     for (int p = 0; p < ND*2; p++)
     {
          GRBLinExpr ptot = 0;
          for (int w = 0; w < ND*ND; ++w)
          {
          ptot += A(p,w)*open[w];
          }
          ostringstream cname;
          cname << "Constriain" << p;
          model.addConstr(ptot == 1, cname.str());
     }
     cout << "finished adding constrains" << endl;
     model.optimize();
    

     cal_result.result.resize(ND*ND,1);
     for (int i = 0; i < ND*ND; i++)
     {
          cal_result.result(i) = open[i].get(GRB_DoubleAttr_X);
     }

     cal_result.obj = model.get(GRB_DoubleAttr_ObjVal);

     } catch(GRBException e) {
     cout << "Error code = " << e.getErrorCode() << endl;
     cout << e.getMessage() << endl;
     } catch(...) {
     cout << "Exception during optimization" << endl;
     }
     // cal_result.result.print("get result: ");
     // cout << "obj: " << cal_result.obj << endl;
     mat cost_mat;
     cost_mat = optT.t()*cal_result.result + cal_result.result.t()*Pbar*cal_result.result;
     result.cost = cost_mat(0,0);
     // cost.print("cost: ");
     mat cal_OptT;
     cal_OptT.copy_size(cal_result.result);
     for (int i = 0; i < cal_OptT.n_rows; i++)
     {
          cal_OptT(i,0) = optT(i,0) * cal_result.result(i,0);
     }
     // cal_OptT.print("cal_OptT:");
     result.maxOptT = cal_OptT.max();
     // cout << "maxOptT: " << maxOptT << endl;
     uvec tempInd = find(cal_result.result == 1);
     
     tempInd += 1;
     
     result.assign.copy_size(tempInd);
     for (int i = 0; i < result.assign.n_rows; i++)
     {
          result.assign(i) = tempInd(i)%ND;
          if (result.assign(i) == 0)
          {
               result.assign(i) = ND;
          }
          
     }
     // assign.print("assign: ");

     return result;
}

// int main() {
//      int ND = 3;
//      mat optT;
//      optT.load("../../../../../Downloads/swarm_matlab/defV/optT.txt");
//      mat Pbar;
//      Pbar.load("../../../../../Downloads/swarm_matlab/defV/Pbar.txt");
//      optT.print("optT: ");
//      Pbar.print("Pbar: ");
//      auto start = high_resolution_clock::now();
//      goal_assign a = defGoalAssignMIQP(optT, Pbar, ND);
//      auto stop = high_resolution_clock::now();
//      auto duration = duration_cast<microseconds>(stop - start);
//      cout << "It takes micro-s: " << duration.count() <<endl;
//      a.assign.print("assign: ");
//      cout << "cost: " << a.cost << endl;
//      cout << "maxOptT: " << a.maxOptT << endl;


// }