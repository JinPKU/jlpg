#include <iostream>
#include <eigen3/Eigen/Dense>
#include "./src/problem.cpp"
using namespace std;

Real f(Vec x){return x(0);}
Vec gradf(Vec x){return x;}



int main(){
    Problem<Vec>  p;//(f,gradf,f,gradf);
    Vec x(3);
    p.f = f;
    p.h = f;
    p.mu = 1;
    x << 2,3,4;
    cout << p.value(x) << endl;
}