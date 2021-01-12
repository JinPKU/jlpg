// constructing a lasso problem
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"
#include "continuation.hpp"

using namespace std;

int main(){
    Mat A = Mat::Random(3,3);
    Vec u = Vec::Random(3);
    cout << u << endl;
    Vec b = A*u;
    Vec x(3); x<<0,0,0;
    Problem<Vec> p(LS(A,b),L1_NORM);
    p.mu = 0.0001;
    Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    opts.setClassical();
    ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);

    Outputs out;
    
    x = cont_pgm(p, x, con_opts, out);
    // x = pgm(p, x, opts, out);
    cout << x << endl;
    cout << "cputime=" << out.cputime << endl;
    cout << "f=" << out.F_cur << "; nrmG=" << out.nrmG << "; Flag=" << out.Flag << endl;
}
