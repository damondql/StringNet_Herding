#include <armadillo>
#include "AllParameters.cpp"
using namespace arma;

vec lambdaInterSecCircLine(vec rc, double rhoc, vec r1, double mL, double cL, double drx, double dry) {
    double lambda01, lambda02;
    if(mL < 1e16) {
        double b = -2*(rc(0)-mL*cL+mL*rc(1));
        double a=(1+pow(mL,2));
        double c=pow(rc(0),2)+pow(cL,2)+pow(rc(1),2)-2*cL*rc(1)-pow(rhoc,2);
        double x_int01 = (-b-sqrt(pow(b,2)-4*a*c))/(2*a);
        lambda01=(x_int01-r1(0))/drx;
        double x_int02 = (-b+sqrt(pow(b,2)-4*a*c))/(2*a);
        lambda02=(x_int02-r1(0))/drx;
    } else {
        double x_int01 = r1(0);
        double y_int01 = rc(1)+sqrt(pow(rhoc,2)-pow((x_int01-rc(0)),2));
        double y_int02 = rc(1)-sqrt(pow(rhoc,2)-pow((x_int01-rc(0)),2));
        lambda01 = (y_int01-r1(1))/dry;
        lambda02 = (y_int02-r1(1))/dry;
    }
    vec lambda = {lambda01, lambda02};
    return lambda;
}