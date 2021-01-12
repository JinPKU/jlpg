#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"
#include "example/example_gl.hpp"

//#include "./src/solver.cpp"
using namespace std;
//using namespace Solver;


int main(){
<<<<<<< HEAD
//    Vec x(6); x<<3,1,2,4,5,6;
//    cout << L0_BALL(4).proxh(x,0) << endl;
    Mat x0(2,2);  x0<<1,2,3,4; 
    Problem<Mat> p;
    example_gl(p);
 	Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    opts.setClassical();
    Mat x = pgm(p, x0, opts);
=======
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
>>>>>>> d9713fb72db1ad6525d04dbb3bced37781f35332
    cout << x << endl;
}
