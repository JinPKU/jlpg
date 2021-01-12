#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"
#include "example/example_gl.hpp"

//#include "./src/solver.cpp"
using namespace std;
//using namespace Solver;


int main(){
    Vec x(3); x<<999,997,995;
    Mat A(2,3); A<<0,0,0,0,0,0;
    Vec b(2); b << 0,0;
    Vec u(3); u << 1,1,1;
    Problem<Vec> p(LS(A,b), L1_NORM.shift(u));
    cout << p.h(x);
    p.mu = 0.1;
	Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    opts.setClassical();
    x = pgm(p, x, opts);
    cout << x << endl;
    p.mu = 0.001;
    x = pgm(p, x, opts);
    cout << x << endl;
}
