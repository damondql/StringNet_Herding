#pragma once
#include <armadillo>

using namespace arma;

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

// int main(){
//     mat xA(2,1);
//     xA.col(0) = {4,1};
//     mat xD = {{ 1,4,6},{ 1,4,1}};
//     xA.print("xA: ");
//     xD.print("xD: ");
//     int a= poly_contain(xD, xA);
//     cout << "contain result: " << a << endl;
// }