// logistic regression with sparsity 
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"
#include "continuation.hpp"

using namespace std;

Real h(Vec x){
    return L2_BALL_H(x.segment<10>(0), 1) + Linf_BALL_H(x.segment<10>(10),1);
}
Vec proxh(Vec x, Real t){
    Vec u  = x; 
    u.segment<10>(0) = L2_BALL_PROXH(x.segment<10>(0),1,t);
    u.segment<10>(10) = Linf_BALL_PROXH(x.segment<10>(10),1,t);
    return u;
}
int main(){
    Mat A = Mat::Random(20,20);
    Vec u = Vec::Random(20);
    cout << u << endl;
    Vec b = A*u;
    Vec x = Vec::Zero(20);
    prox_pair<Vec> my_prox(h,proxh);
    Problem<Vec> p(LS(A,b),my_prox);
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
