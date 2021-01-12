// constructing a blocked proximal problem

#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"

using namespace std;

int main(){
    Mat A = Mat::Random(256,512);
    Vec u = Vec::Random(512);
    Vec u0 = 10*(u.array().sign())*(u.array().abs() - 0.9).max(0);
    Vec b = A*u0;
    Vec x = Vec::Zero(512);
    Problem<Vec> p(LS(A,b),L1_NORM);
<<<<<<< HEAD
    Options opts(10000, 1e-10, 1e-8, 1e-1, 5e-1);
    opts.setClassical();

    p.mu = 1e2;
    x = pgm(p,x,opts);

    p.mu = 1e0;
    x = pgm(p,x,opts);

    p.mu = 1e-2;
    x = pgm(p,x,opts);

    p.mu = 1e-4;
    x = pgm(p,x,opts);


    int cnt = 0;
    for(int i = 0; i < 512; i++){
        if(abs(x[i])<1e-4) cnt ++ ;
    }
    cout << cnt / double(512) << endl;
    cout << (x-u0).lpNorm<1>()/u0.lpNorm<1>() << endl;
=======
    p.mu = 0.01;
    cout << p.value(u0) << endl;
    p.mu*=(1<<15);
    Options opts(10000, 1e-8*(1<<15), 1e-6*(1<<15), 1e-0, 5e-1);
    opts.setClassical();

    for(int i = 0; i < 15; i ++){
        x = pgm(p,x,opts);
        p.mu/=2, opts.ftol/=2,opts.gtol/=2;
    }
>>>>>>> f2048ea3be15149c6cd4e01b4d7fc22bf04093b8
}
