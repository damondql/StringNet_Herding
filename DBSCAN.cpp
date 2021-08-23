#include <armadillo>
#include <iostream>
using namespace std;
using namespace arma;





vec DBSCAN(mat X, double epsilon, int MinPts)
{
    int C = 0;
    int n = X.n_rows;
    vec IDX = zeros<vec>(n,1);

    mat D(n,n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            D(i,j) = sqrt(pow(X(i,0)-X(j,0),2) + pow(X(i,1)-X(j,1),2) );
        }
        
    }
    mat visited = zeros<mat>(n,1);
    mat isnoise = zeros<mat>(n,1);

    for (int i = 0; i < n; i++)
    {
        if(!visited(i))
        {
            visited(i) = 1;
            uvec Nieghbors = arma::find(D.row(i) <= epsilon);
            if(Nieghbors.n_elem -1 < MinPts)
            {
                isnoise(i) = 1;
            } else
            {
                C++;
                //Expand Cluster
            }
        }
    }
    
    
}