#include <iostream>
#include <eigen3/Eigen/Dense>
#include "./src/problem.cpp"
using namespace std;

Real f(Vec x){return x(0);}
Vec gradf(Vec x){return x;}



int main(){
    Problem<Vec>  p(L1_NORM);
    p.f = f;
    p.gradf = gradf;
    p.mu = 1;
    Vec x(3);
    x << -2,0.5,4;
    cout << p.proxh(x,1)<<endl;
    cout << p.value(x) << endl;
}