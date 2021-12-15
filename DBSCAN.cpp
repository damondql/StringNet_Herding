#include <armadillo>
#include <iostream>
using namespace std;
using namespace arma;


vec DBSCAN(mat X, double epsilon, int MinPts)
{
    int C = -1;
    int n = X.n_cols;
    vec IDX = zeros<vec>(n,1);

    mat D(n,n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            D(i,j) = sqrt(pow(X(0,i)-X(0,j),2) + pow(X(1,i)-X(1,j),2) );
        }
        
    }
    mat visited = zeros<mat>(n,1);
    mat isnoise = zeros<mat>(n,1);

    for (int i = 0; i < n; i++)
    {
        if(!visited(i))
        {
            visited(i) = 1;
            uvec Neighbors = arma::find(D.row(i) <= epsilon);
            if(Neighbors.n_elem -1 < MinPts)
            {
                isnoise(i) = 1;
            } else
            {
                C++;
                //Expand Cluster
                IDX(i) = C;
                int k = 0;
                while (1)
                {
                    double j = Neighbors(k);
                    if(!visited(j))
                    {
                        visited(j) = 1;
                        uvec Neighbors2 = arma::find(D.row(j) <= epsilon);
                        if(Neighbors2.n_elem -1 >= MinPts)
                        {
                            Neighbors = arma::join_cols(Neighbors, Neighbors2);
                        }
                    }
                    if(IDX(j) == 0)
                    {
                        IDX(j) = C;
                    }
                    k++;
                    
                    if(k > Neighbors.n_elem-1)
                    {
                        break;
                    }
                }
                
            }
        }
    }
    return IDX;
    
}


// int main(){
//     mat X;
//     X.load("../dbs/XA_submat.csv");
//     // X = X.t();
//     X.print("X");
//     double epsilon = 31.2254;
//     vec a = DBSCAN(X,epsilon, 4);
//     a.print("a");

// }