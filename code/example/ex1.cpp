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
    p.mu = 10;
    Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    opts.setClassical();

    x = pgm(p,x,opts);
    int cnt = 0;
    for(int i = 0; i < 1000; i++){
        if(abs(x[i])<1e-4) cnt ++ ;
    }
    cout << cnt << endl;
    cout << (x-u0).lpNorm<1>() << endl;
}
