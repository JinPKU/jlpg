// constructing a lasso problem
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"

using namespace std;

int main(){
    Mat A = Mat::Random(3,3);
    Vec u = Vec::Random(3);
    cout << u << endl;
    Vec b = A*u;
    Vec x(3); x<<0,0,0;
    Problem<Vec> p(LS(A,b),Linf_BALL(0.1));
    p.mu = 0.01;
    Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    Outputs out;
    opts.setClassical();
    x = pgm(p,x,opts,out);
    cout << x << endl;
}
