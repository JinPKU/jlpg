// logistic regression with sparsity 
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
    Problem<Vec> p(LOGISTIC(A,b), L1_NORM);
    p.mu = 0.0001;
    Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);
    Outputs out;
    
    Vec x1 = x, x2 = x, x3 = x;
    printf("type & cont. & iters & cputime & fval & optimality & flag \\\\\n");
    
    // Armijo + continuation
    opts.setArmijo(0.1);
    x1 = cont_pgm(p, x1, con_opts, out);
    printf("Armijo & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Nonmonotone + continuation
    opts.setNonmonotone(0.1, 10);
    x2 = cont_pgm(p, x2, con_opts, out);
    printf("Nonmonotone & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);

    // Classical + continuation
    opts.setClassical();
    x3 = cont_pgm(p, x3, con_opts, out);
    printf("Classical & 1 & %d & %g & %g & %g & %d \\\\\n", out.iter, out.cputime, out.F_cur, out.nrmG, out.Flag);
}
