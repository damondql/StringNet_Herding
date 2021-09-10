
// C++ implementation of the approach
#pragma once
#include <armadillo>
#include <bits/stdc++.h>
#define llu long long int
using namespace std;
using namespace arma;

struct Point {

	double x, y;

	bool operator<(Point p)
	{
		return x < p.x || (x == p.x && y < p.y);
	}
};

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

// int main()
// {
//     // Point polygon[] = {{0.0,0.0}, {0.0,10.0}, {10.0,10.0}, {10.0,0.0}};
//     vector<Point> polygon;
//     polygon.push_back({0,0});
//     polygon.push_back({0,10});
//     polygon.push_back({10,10});
//     polygon.push_back({10,0});
//     Point centroid = compute2DPolygonCentroid(polygon);
//     std::cout << "Centroid is (" << centroid.x << ", " << centroid.y << ")\n";
// }


// Driver code
int main()
{

	// Add points
  mat XD;
	XD.load("/home/damon/Downloads/multi_swarm/controlA/XD.txt");
  XD.print("XD");
  vector<Point> points;
  for (int i = 0; i < XD.n_cols; i++)
  {
      Point tempP;
      tempP.x = XD(0,i);
      tempP.y = XD(1,i);
      points.push_back(tempP);
  }

  cout << "number of points in vector: " << points.size() << endl;
	// Find the convex hull
	vector<Point> ans = convex_hull(points);

	// Print the convex hull
	for (int i = 0; i < ans.size(); i++)
		cout << "(" << ans[i].x << ", "
			<< ans[i].y << ")" << endl;

  Point centroid = compute2DPolygonCentroid(ans);
  std::cout << "Centroid is (" << centroid.x << ", " << centroid.y << ")\n";

	return 0;
}


