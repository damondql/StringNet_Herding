#include <armadillo>

using namespace std;
using namespace arma;

struct CommGraph {
    mat W;
    mat Rij_tilde;
};

CommGraph findCommGraphAndFormDist(double N, double shapeld, double R0);

CommGraph findCommGraphAndFormDist(double N, double shapeld, double R0) {
    CommGraph result;
    if (shapeld == 0)
    {
        double Rij0 = R0;
        double Rij1 = (N-1) * R0;
        if (N > 1)
        {
            result.W = ones<mat>(N, N);
            result.Rij_tilde = zeros<mat>(N,N);
            for (int i = 0; i < N; i++)
            {
                result.W(i,i) = 0;
                for (int ii = 0; ii < N; ii++)
                {
                    result.Rij_tilde(i,ii) = abs(i - ii) * R0;
                }
            }
        } else {
            result.W = {0};
        }
    } else if (shapeld == 1) {
      result.Rij_tilde = zeros<mat>(N,N);
        if (N > 1) {
            result.W = zeros<mat> (N,N);
            int NC = N;
            double Rij0 = R0 * sqrt(2 * (1-cos(2*M_PI / NC)));
            double Rij1 = Rij0 * sqrt(2* (1-cos(M_PI - 2*M_PI/NC)));
            mat temp(1,NC*2);
            vec temp_v;
            temp_v = arma::regspace(0,NC-1);
            for (size_t i = 0; i < NC; i++)
            {
                temp(0,i) = temp_v(i);
                temp(0,i+NC) = temp_v(i);
            }
            for (int i = 0; i < NC; i++)
            {
                result.W(i,temp(i+1)) = 1;
                result.W(i,temp(i+2)) = 1;
                result.W(i,temp(i+NC-1)) = 1;
                result.Rij_tilde(i,temp(0,i+1)) = Rij0;
                result.Rij_tilde(i,temp(0,i+NC-1)) = Rij0;
                result.Rij_tilde(i,temp(0,temp(i+2))) = Rij1;
                result.Rij_tilde(i,i) = 0;
            }
        } else {
            result.W = {0};
        }
    } else if (shapeld ==2) {
        result.Rij_tilde = zeros<mat>(N,N);
        if (N > 1) {
            result.W = zeros<mat> (N,N);
            int NC = N-1;
            double Rij0=R0*sqrt(2*(1-cos(2*M_PI/NC)));
            double Rij1=Rij0*sqrt(2*(1-cos(M_PI-2*M_PI/NC)));
            mat temp(1,NC*2);
            vec temp_v;
            temp_v = arma::regspace(0,NC-1);
            for (size_t i = 0; i < NC; i++)
            {
                temp(0,i) = temp_v(i);
                temp(0,i+NC) = temp_v(i);
            }
            for (int i = 0; i < NC; i++)
            {
                result.W(i,temp(i+1)) = 1;
                result.W(i,temp(i+2)) = 1;
                result.W(i,temp(i+NC-1)) = 1;
                result.Rij_tilde(i,temp(0,i+1)) = Rij0;
                result.Rij_tilde(i,temp(0,i+NC-1)) = Rij0;
                result.Rij_tilde(i,temp(0,temp(i+2))) = Rij1;
                result.Rij_tilde(i,i) = 0;
            }
            for (int i = 0; i < N-1; i++) {
                result.W(N-1, i) = 1;
                result.W(i, N-1) = 1;
                result.Rij_tilde(N-1,i) = R0;
                result.Rij_tilde(i, N-1) = R0;
            }
        } else {
            result.W = {0};
        }
    } else if (shapeld == 3) {
      result.Rij_tilde = zeros<mat>(N,N);
        if (N > 1) {
            result.W = zeros<mat> (N,N);
            int NC = N;
            double Rij0 = R0*sqrt(2*(1-cos(M_PI/(NC-1))));
            double Rij1 = Rij0*sqrt(2*(1-cos(M_PI-M_PI/(NC-1))));
            mat temp(1,NC*2);
            vec temp_v;
            temp_v = arma::regspace(0,NC-1);
            for (size_t i = 0; i < NC; i++)
            {
                temp(0,i) = temp_v(i);
                temp(0,i+NC) = temp_v(i);
            }
            for (int i = 0; i < NC; i++)
            {
                result.W(i,temp(i+1)) = 1;
                result.W(i,temp(i+2)) = 1;
                result.W(i,temp(i+NC-1)) = 1;
                result.Rij_tilde(i,temp(0,i+1)) = Rij0;
                result.Rij_tilde(i,temp(0,i+NC-1)) = Rij0;
                result.Rij_tilde(i,temp(0,temp(i+2))) = Rij1;
                result.Rij_tilde(i,i) = 0;
            }
            result.Rij_tilde(0, N-1) = 0;
            result.Rij_tilde(N-1, 0) = 0;
            result.W(0, N-1) = 0;
            result.W(N-1, 0) = 0;
        } else {
            result.W = {0};
        }

    } else if (shapeld == 4) {
        result.Rij_tilde = zeros<mat>(N,N);
        if (N > 1) {
            result.W = zeros<mat> (N,N);
            int NC = N-1;
            double Rij0 = R0*sqrt(2*(1-cos(M_PI/(NC-1))));
            double Rij1 = Rij0*sqrt(2*(1-cos(M_PI-M_PI/(NC-1))));
            mat temp(1,NC*2);
            vec temp_v;
            temp_v = arma::regspace(0,NC-1);
            for (size_t i = 0; i < NC; i++)
            {
                temp(0,i) = temp_v(i);
                temp(0,i+NC) = temp_v(i);
            }
            for (int i = 0; i < NC; i++)
            {
                result.W(i,temp(i+1)) = 1;
                result.W(i,temp(i+2)) = 1;
                result.W(i,temp(i+NC-1)) = 1;
                result.Rij_tilde(i,temp(0,i+1)) = Rij0;
                result.Rij_tilde(i,temp(0,i+NC-1)) = Rij0;
                result.Rij_tilde(i,temp(0,temp(i+2))) = Rij1;
                result.Rij_tilde(i,i) = 0;
            }
            for (int i = 0; i < N-1; i++) {
                result.W(N-1, i) = 1;
                result.W(i, N-1) = 1;
                result.Rij_tilde(N-1,i) = R0;
                result.Rij_tilde(i, N-1) = R0;
            }
            result.Rij_tilde(0, N-2) = 0;
            result.Rij_tilde(N-2, 0) = 0;
            result.W(0, N-2) = 0;
            result.W(N-2, 0) = 0;
        } else {
            result.W = {0};
        }
    } else if (shapeld == 5)
    {
        result.W = ones<mat>(N,N);
        result.Rij_tilde = zeros<mat>(N,N);
        for (int i = 0; i < N; i++)
        {
            result.W(i,i) = 0;
            for (int j = 0; j < N; j++)
            {
                result.Rij_tilde(i,j) = R0 * abs(i-j);
            }
            
        }
        
    } else if (shapeld == 6)
    {
        result.W = ones<mat>(N,N);
        result.Rij_tilde = zeros<mat>(N,N);
        for (int i = 0; i < N; i++)
        {
            result.W(i,i) = 0;
            if (i == 0)
            {
                cout<< "0" <<endl;
                result.Rij_tilde(0,1) = R0;
                result.Rij_tilde(0,2) = R0;
                for (int ii = 3; ii < N; ii++)
                {
                    result.Rij_tilde(0,ii) = sqrt(2) * R0 + (ii-3) * R0;
                }
                

                
            } else if (i == 1)
            {
                cout<< "1" <<endl;
                result.Rij_tilde(1,0) = R0;
                result.Rij_tilde(1,2) = sqrt(2) * R0;
                for (int ii = 3; ii < N; ii++)
                {
                    result.Rij_tilde(1,ii) = sqrt(pow(R0,2)+pow(R0*(ii-3),2)-2*R0*(R0*(ii-3))*cos(3*M_PI/4));
                }
            } else if (i == 2)
            {
                cout<< "2" <<endl;
                result.Rij_tilde(2,0) = R0;
                result.Rij_tilde(2,1) = sqrt(2) * R0;
                for (int ii = 3; ii < N; ii++)
                {
                    result.Rij_tilde(2,ii) =sqrt(pow(R0,2)+pow((R0*(ii-3)),2)-2*R0*(R0*(ii-3))*cos(3*M_PI/4));
                }
                
            } else {
                cout<< "else" <<endl;
                result.Rij_tilde(i,0) = result.Rij_tilde(0,i);
                result.Rij_tilde(i,1) = result.Rij_tilde(1,i);
                result.Rij_tilde(i,2) = result.Rij_tilde(2,i);
                for (int ii = 4; ii < N; ii++)
                {
                    result.Rij_tilde(i,ii) = (ii - i) * R0;
                }
            }
        }
        
    } else {
        result.W = {0};
        cout << "No valid shape specified!" << endl;
    }
    
  return result;
}

