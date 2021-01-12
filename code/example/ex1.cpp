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
    Options opts(10000, 1e-10, 1e-8, 1e-1, 5e-1);
    opts.setNonmonotone(0.9, 10);

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
}
