// nuclear norm 
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "jlpg.hpp"
#include <cmath>
using namespace std;

int main(){
    Mat A = Mat::Random(200,100);
    Mat B = Mat::Random(200,50);
    Mat X = Mat::Zero(100,50);
    Problem<Mat> p(LS(A,B),NUCLEAR_NORM);
    p.mu = 0.01;
    Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1);
    opts.setClassical();
    ContOptions con_opts(opts, 10, 0.1, 100, 1e2, 1e2, 1e-1, 1e-1);

    Outputs out;
    
    X = cont_pgm(p, X, con_opts, out);
    // x = pgm(p, x, opts, out);
    //cout << x << endl;
    cout << "cputime=" << out.cputime << endl;
    cout << "f=" << out.F_cur << "; nrmG=" << out.nrmG << "; Flag=" << out.Flag << endl;
}
