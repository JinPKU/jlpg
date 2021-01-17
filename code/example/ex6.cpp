// simple box constrainted optimization.
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../include/jlpg.hpp"
#include <cmath>
using namespace std;

Real f(Vec x){
    return exp(x[0]*x[3] - x[1]*x[2]) + exp(x[2]+x[3]);
}
Vec gradf(Vec x){
    Vec u(4);
    u[0] = x[3]*exp(x[0]*x[3] - x[1]*x[2]);
    u[1] = -x[2]*exp(x[0]*x[3] - x[1]*x[2]);
    u[2] = -x[1]*exp(x[0]*x[3] - x[1]*x[2]) + exp(x[2]+x[3]);
    u[3] = x[0]*exp(x[0]*x[3] - x[1]*x[2]) + exp(x[2]+x[3]);
    return u;
    
}
int main(){
    Vec lx(4); lx<<0,0,0,0;
    Vec ux(4); ux<<3,2,2,2;
    Vec x(4); x<<1,1,1,1;
    grad_pair<Vec> my_grad(f,gradf);
    Problem<Vec> p(my_grad,BOX(lx,ux));
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