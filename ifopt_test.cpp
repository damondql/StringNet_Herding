#include <iostream>
#include <math.h>
#include <vector>
#include <typeinfo>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <armadillo>

using namespace Eigen;
int main()
{
// First step, defines symbols:
struct g_tag {};  static const symbolic::SymbolExpr<g_tag> g;


double gk1 = 2;
double gk2 = 5;
arma::mat a(3,3,arma::fill::zeros);
a.print("a:");


}
