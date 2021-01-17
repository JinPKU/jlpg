
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
    // ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);

    Outputs out;
    
    Vec x1 = x, x2 = x, x3 = x;
    printf("type & cont. & iters & cputime & fval & optimality & flag \\\\\n");
    
    // Armijo
    opts.setArmijo(0.1);
    x1 = pgm(p, x1, opts, out);
    printf("Armijo & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Nonmonotone
    opts.setNonmonotone(0.1, 10);
    x2 = pgm(p, x2, opts, out);
    printf("Nonmonotone & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Classical
    opts.setClassical();
    x3 = pgm(p, x3, opts, out);
    printf("Classical & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);
}
