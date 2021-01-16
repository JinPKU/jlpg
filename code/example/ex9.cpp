
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"

using namespace std;

int main(){
    Mat A = Mat::Random(10,10);
    Vec u = Vec::Random(10);
    cout << u << endl;
    Vec b = A*u;
    Vec x = Vec::Zero(10);
    Problem<Vec> p(LS(A,b), Linf_NORM);
    p.mu = 1e-2;
    Options opts(10000, 1e-8, 1e-6, 1e-2, 5e-1);
    opts.setClassical();
    // ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);

    Outputs out;
    
    // x = cont_pgm(p, x, con_opts, out);
    x = pgm(p, x, opts, out);
    cout << x << endl;
    cout << "cputime=" << out.cputime << endl;
    cout << "f=" << out.F_cur << "; nrmG=" << out.nrmG << "; Flag=" << out.Flag << endl;
}
