#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"

//#include "./src/solver.cpp"
using namespace std;
//using namespace Solver;
Real f(Vec x){return x(0);}
Vec gradf(Vec x){return x;}



int main(){
    Mat A(3,2); A<<1,2,3,4,5,6;
    Vec x(2); x<<1,2;
    Vec b(3); b<<1,1,1;
    Problem<Vec> p(LS(A,b),L0_NORM);
	Options opts(1,1, 10000, 1e-8, 1e-6, 1e-2);
    cout << pgm(p, x, opts) << endl;
}
