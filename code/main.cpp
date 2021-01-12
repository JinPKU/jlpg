#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"
#include "example/example_gl.hpp"

//#include "./src/solver.cpp"
using namespace std;
//using namespace Solver;
Real f(Vec x){return x(0);}
Vec gradf(Vec x){return x;}



int main(){
    Mat x0(2,2);  x0<<1,2,3,4; 
    Problem<Mat> p;
    example_gl(p);
	Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    opts.setClassical();
    Mat x = pgm(p, x0, opts);
    cout << x << endl;
}
