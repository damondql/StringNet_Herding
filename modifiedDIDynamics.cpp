#include <armadillo>

using namespace arma;

mat Xdot(mat X, mat U, double C_d){
    int N = (int)X.n_rows/4;
    mat f  = zeros(N*4,1);
    for (int i = 1; i <= N; i++)
    {
        f.submat(4*i-4,0,4*i-3,0) = X.submat(4*i-2,0,4*i-1,0);
        f.submat(4*i-2,0,4*i-1,0) = U.submat(2*i-2,0,2*i-1,0) - C_d*arma::norm(X.submat(4*i-2,0,4*i-1,0))*X.submat(4*i-2,0,4*i-1,0);
    }
    return f;
}

mat modifiedDIDynamics(mat X0, mat U, double dt, double C_d){
    
    mat k1 = Xdot(X0,U, C_d);
    mat k2 = Xdot(X0 + dt/2*k1,U, C_d);
    mat k3 = Xdot(X0 + dt/2*k2,U, C_d);
    mat k4 = Xdot(X0 + dt*k3, U, C_d);
    mat X1 = X0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    return X1;
}