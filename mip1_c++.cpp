/* Copyright 2021, Gurobi Optimization, LLC */

/* This example formulates and solves the following simple QP model:

     minimize    x^2 + x*y + y^2 + y*z + z^2 + 2 x
     subject to  x + 2 y + 3 z >= 4
                 x +   y       >= 1
                 x, y, z non-negative

   It solves it once as a continuous model, and once as an integer model.
*/

#include "gurobi_c++.h"
#include <armadillo>

using namespace std;
using namespace arma;


struct goal_assign
{
     vec assign;
     double cost;
     double maxOpt;
};


struct gurobi_result
{
     mat result;
     double obj;
};

gurobi_result defGoalAssignMIQP(mat optT, mat Pbar, int ND)
{
     gurobi_result result;
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

     cout << "finish variable" <<endl;
     std::vector<GRBVar> variables;
     
     variables = {a,b,c,d,e,f,g,h,i};
     cout <<"make varibles vector" << endl;
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
     cout << "finish obj" << endl;


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

     mat cal_result = {a1,b1,c1,d1,e1,f1,g1,h1,i1};
     result.result = cal_result.t();

     result.obj = model.get(GRB_DoubleAttr_ObjVal);
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


     return result;
}

int main() {
     int ND = 3;
     mat optT;
     optT.load("../defV/optT.txt");
     mat Pbar;
     Pbar.load("../defV/Pbar.txt");
     optT.print("optT: ");
     Pbar.print("Pbar: ");
     gurobi_result a = defGoalAssignMIQP(optT, Pbar, ND);
     a.result.print("get result: ");
     cout << "obj: " << a.obj << endl;
     mat cost;
     cost = optT.t()*a.result + a.result.t()*Pbar*a.result;
     cost.print("cost: ");
     mat cal_OptT;
     cal_OptT.copy_size(a.result);
     for (int i = 0; i < cal_OptT.n_rows; i++)
     {
          cal_OptT(i,0) = optT(i,0) * a.result(i,0);
     }
     cal_OptT.print("cal_OptT:");
     double maxOptT;
     maxOptT = cal_OptT.max();
     cout << "maxOptT: " << maxOptT << endl;
     uvec tempInd = find(a.result == 1);
     tempInd.print("tempInd: ");
     tempInd += 1;
     tempInd.print("tempInd: ");
     vec assign;
     assign.copy_size(tempInd);
     for (int i = 0; i < assign.n_rows; i++)
     {
          assign(i) = tempInd(i)%ND;
          if (assign(i) == 0)
          {
               assign(i) = ND;
          }
          
     }
     assign.print("assign: ");

}