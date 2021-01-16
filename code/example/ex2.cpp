// Group Lasso problem

#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../include/jlpg.hpp"

using namespace std;

int main(){
    int m = 128, n = 256, l = 2;
    Mat A = Mat::Random(m, n);
    Mat u = Mat::Random(n, l);
    Mat u0 = 10 * (u.array().sign()) * (u.array().abs() - 0.9).max(0);
    Mat b = A * u0;
    Mat x = Mat::Zero(n, l);
    
    Problem<Mat> p(LS(A,b), L12_NORM);
    p.mu = 1e-2;

    Options opts(200, 1e-10, 1e-8, 1e-1, 5e-1);
    ContOptions con_opts(opts, 100, 0.1, 100, 1e6, 1e6, 1e-1, 1e-1);
    Outputs out;

    Mat x1 = x, x2 = x, x3 = x;
    printf("type & cont. & iters & cputime & fval & optimality & flag \\\\\n");
    
    // Armijo + continuation
    opts.setArmijo(0.9);
    x1 = cont_pgm(p, x1, con_opts, out);
    printf("Armijo & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Nonmonotone + continuation
    opts.setNonmonotone(0.9, 5);
    x2 = cont_pgm(p, x2, con_opts, out);
    printf("Nonmonotone & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Classical + continuation
    opts.setClassical();
    x3 = cont_pgm(p, x3, con_opts, out);
    printf("Classical & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    
}
