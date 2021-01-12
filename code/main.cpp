#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"

//#include "./src/solver.cpp"
using namespace std;
//using namespace Solver;
Real f(Vec x){return x(0);}
Vec gradf(Vec x){return x;}



int main(){
    Mat A(3,2);
    A<<1,2,3,4,5,6;
    Mat x(2,2);  x<<1,2,3,4; 
    Mat b(3,2); b<<1,1,1,1,1,1;
    // Vec b(2);
    // b<<1,2;
    Problem<Mat> p(LS(A,b),NO_PROX);
    // p.mu = 1;
    // Vec x(3);
    // x << -2,0.5,1;
    // cout << p.f(x) << endl;
	Options opts(true, true, 10000, 1e-8, 1e-6, 1e-2);
    cout << pgm(p, x, opts) << endl;
}
