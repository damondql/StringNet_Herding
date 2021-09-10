#pragma once
#include <armadillo>
#include <bits/stdc++.h>
using namespace std;
using namespace arma;

 
struct Point {

	double x, y;

	bool operator<(Point p)
	{
		return x < p.x || (x == p.x && y < p.y);
	}
};
 
bool line_intersects(double p1_x1, double p1_y1, double p1_x2, double p1_y2,
                     double p2_x1, double p2_y1, double p2_x2, double p2_y2) {
    double z1 = (p1_x2 - p1_x1) * (p2_y1 - p1_y1) - (p2_x1 - p1_x1) * (p1_y2 - p1_y1);
    double z2 = (p1_x2 - p1_x1) * (p2_y2 - p1_y1) - (p2_x2 - p1_x1) * (p1_y2 - p1_y1);
    int os;
    if ((z1 > 0 && z2 < 0) || (z1 < 0 && z2 > 0)) {
        os = 1;
    } else {
        os = 0;
    }
    double z3 = (p2_x2 - p2_x1) * (p1_y1 - p2_y1) - (p1_x1 - p2_x1) * (p2_y2 - p2_y1);
    double z4 = (p2_x2 - p2_x1) * (p1_y2 - p2_y1) - (p1_x2 - p2_x1) * (p2_y2 - p2_y1);
    if ((z3 > 0 && z4 < 0) || (z3 < 0 && z4 > 0)) {
        os++;
    }
    int colls;
    if (os == 2) {
        colls = 1;
    } else {
        colls = 0;
    }
    /*if (colls == 1) {
        printf("line intersection between (%.5f, %.5f), (%.5f, %.5f) ",
               p1_x1, p1_y1, p1_x2, p1_y2);
        printf("and (%.5f, %.5f), (%.5f, %.5f)\n", p2_x1, p2_y1, p2_x2, p2_y2);
    } else {
        printf("NO line intersection between (%.5f, %.5f), (%.5f, %.5f) ",
               p1_x1, p1_y1, p1_x2, p1_y2);
        printf("and (%.5f, %.5f), (%.5f, %.5f)\n", p2_x1, p2_y1, p2_x2, p2_y2);
    }*/
    /*bool inters = false;
    if (colls == 1) {
        inters = true;
    }*/
    return colls;
}

int poly_inters(mat poly1, mat poly2) {
    int if_inters = 0;
    for (int i = 0; i < poly1.n_cols; i++) {
        int a = (i + 1) % (int)poly1.n_cols;
        for (int j = 0; j < poly2.n_cols; j++) {
            int b = (j + 1) % (int)poly2.n_cols;
            if_inters = if_inters + line_intersects(poly1(0,i), poly1(1,i),
                                                    poly1(0,a), poly1(1,a),
                                                    poly2(0,j), poly2(1,j),
                                                    poly2(0,b), poly2(1,b));
        }
    }
    return if_inters;
}

double point_contain(double p1_x1, double p1_y1, double p1_x2, double p1_y2,
                     double p2_x1, double p2_y1) {
    double z1 = (p1_x2 - p1_x1) * (p2_y1 - p1_y1) - (p2_x1 - p1_x1) * (p1_y2 - p1_y1);
    return z1;
}

int poly_contain(mat poly1, mat poly2) {
    int if_contain = 0;
    double contain_result;
    for (int i = 0; i < poly2.n_cols; i++) {
        int pos = 0;
        int neg = 0;
        int zero = 0;
        for (int j = 0; j < poly1.n_cols; j++) {
            int a = (j + 1) % (int)poly1.n_cols;
            contain_result = point_contain(poly1(0,j), poly1(1,j),
                                           poly1(0,a), poly1(1,a),
                                           poly2(0,i), poly2(1,i));
            if (contain_result > 0) {
                pos++;
            } else if (contain_result < 0) {
                neg++;
            } else {
                zero++;
            }
        }
        if ((pos + zero == poly1.n_cols) || (neg + zero == poly1.n_cols)) {
            if_contain++;
        }
    }
    return if_contain;
}



// Cross product of two vectors OA and OB
// returns positive for counter clockwise
// turn and negative for clockwise turn
double cross_product(Point O, Point A, Point B)
{
	return (A.x - O.x) * (B.y - O.y)
		- (A.y - O.y) * (B.x - O.x);
}

// Returns a list of points on the convex hull
// in counter-clockwise order
vector<Point> convex_hull(vector<Point> A)
{
	int n = A.size(), k = 0;

	if (n <= 3)
		return A;

	vector<Point> ans(2 * n);

	// Sort points lexicographically
	sort(A.begin(), A.end());

	// Build lower hull
	for (int i = 0; i < n; ++i) {

		// If the point at K-1 position is not a part
		// of hull as vector from ans[k-2] to ans[k-1]
		// and ans[k-2] to A[i] has a clockwise turn
		while (k >= 2 && cross_product(ans[k - 2],
						ans[k - 1], A[i]) <= 0)
			k--;
		ans[k++] = A[i];
	}

	// Build upper hull
	for (size_t i = n - 1, t = k + 1; i > 0; --i) {

		// If the point at K-1 position is not a part
		// of hull as vector from ans[k-2] to ans[k-1]
		// and ans[k-2] to A[i] has a clockwise turn
		while (k >= t && cross_product(ans[k - 2],
						ans[k - 1], A[i - 1]) <= 0)
			k--;
		ans[k++] = A[i - 1];
	}

	// Resize the array to desired size
	ans.resize(k - 1);

	return ans;
}


Point compute2DPolygonCentroid(vector<Point> vertices)
{
    Point centroid = {0, 0};
    double signedArea = 0.0;
    double x0 = 0.0; // Current vertex X
    double y0 = 0.0; // Current vertex Y
    double x1 = 0.0; // Next vertex X
    double y1 = 0.0; // Next vertex Y
    double a = 0.0;  // Partial signed area

    // For all vertices except last
    int i=0;
    for (i=0; i<vertices.size()-1; ++i)
    {
        x0 = vertices[i].x;
        y0 = vertices[i].y;
        x1 = vertices[i+1].x;
        y1 = vertices[i+1].y;
        a = x0*y1 - x1*y0;
        signedArea += a;
        centroid.x += (x0 + x1)*a;
        centroid.y += (y0 + y1)*a;
    }

    // Do last vertex separately to avoid performing an expensive
    // modulus operation in each iteration.
    x0 = vertices[i].x;
    y0 = vertices[i].y;
    x1 = vertices[0].x;
    y1 = vertices[0].y;
    a = x0*y1 - x1*y0;
    signedArea += a;
    centroid.x += (x0 + x1)*a;
    centroid.y += (y0 + y1)*a;

    signedArea *= 0.5;
    centroid.x /= (6.0*signedArea);
    centroid.y /= (6.0*signedArea);

    return centroid;
}
 
// Driver program to test above functions
// int main()
// {
//     mat XD;
//     XD.load("/home/damon/Downloads/multi_swarm/controlA/XD.txt");
//     XD.print("XD");
//     vector<Point> points;
//     for (int i = 0; i < XD.n_cols; i++)
//     {
//         Point tempP;
//         tempP.x = XD(0,i);
//         tempP.y = XD(1,i);
//         points.push_back(tempP);
//     }
//     int n = XD.n_cols;
    
//     // int n = sizeof(points)/sizeof(points[0]);
//     convexHull(points, n);
//     return 0;
// }


// int main(){
//     mat xA(2,1);
//     xA.col(0) = {4,1};
//     mat xD = {{ 1,4,6},{ 1,4,1}};
//     xA.print("xA: ");
//     xD.print("xD: ");
//     int a= poly_contain(xD, xA);
//     cout << "contain result: " << a << endl;
// }