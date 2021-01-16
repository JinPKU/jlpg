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
    opts.setClassical();
    ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);

    Outputs out;
    
    x = cont_pgm(p, x, con_opts, out);
    // x = pgm(p, x, opts, out);
    cout << x << endl;
    cout << "cputime=" << out.cputime << endl;
    cout << "f=" << out.F_cur << "; nrmG=" << out.nrmG << "; Flag=" << out.Flag << endl;
}