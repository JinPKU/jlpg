// constructing a blocked proximal problem

#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"

using namespace std;

int main(){
    Mat A = Mat::Random(500,1000);
    Vec u = Vec::Random(1000);
    Vec u0 = 10*(u.array().sign())*(u.array().abs() - 0.9).max(0);
    Vec b = A*u0;
    Vec x = Vec::Zero(1000);
    Problem<Vec> p(LS(A,b),L1_NORM);
    p.mu = 0.01;
    cout << p.value(u0) << endl;
    p.mu*=(1<<15);
    Options opts(10000, 1e-8*(1<<15), 1e-6*(1<<15), 1e-0, 5e-1);
    opts.setClassical();

    for(int i = 0; i < 15; i ++){
        x = pgm(p,x,opts);
        p.mu/=2, opts.ftol/=2,opts.gtol/=2;
    }
}
