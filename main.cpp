#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"

//#include "./src/solver.cpp"
using namespace std;
//using namespace Solver;
Real f(Vec x){return x(0);}
Vec gradf(Vec x){return x;}



int main(){
    Mat A(2,2);
    A<<1,2,3,4;
    Vec b(2);
    b<<1,2;
    Problem<Vec> p(LS(A,b),L1_NORM);
    p.mu = 1;
    Vec x(2);
    x << -2,0.5;
    cout << p.proxh(x,1)<<endl;
    cout << p.value(x) << endl;
	Options opts(true, true, 10000, 1e-8, 1e-6, 1e-2);
    cout << pgm(p, x, opts);
}