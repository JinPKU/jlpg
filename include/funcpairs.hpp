#ifndef __FUNCPAIRS_HPP__
#define __FUNCPAIRS_HPP__
#include <functional>
#include <iostream>
#include <eigen3/Eigen/Dense>
// func pairs provide an interface for problem, such as  LS, norm ball of norm 
using namespace std;
typedef double Real;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;

template <typename T>
class grad_pair{
    public:
    grad_pair(){}
    grad_pair(function<Real(T)> f, function<T(T)> gradf);

    function<Real(T)> f;
    function<T(T)> gradf;
};

template <typename T>
class prox_pair{
    public:
    prox_pair(){}
    prox_pair(function<Real(T)> h, function<T(T, Real)> proxh);

    function<Real(T)> h;
    function<T(T,Real)> proxh;
};


template<typename T>
grad_pair<T>:: grad_pair(function<Real(T)> f, function<T(T)> gradf){
    this->f = f; this->gradf = gradf;
}


template<typename T>
prox_pair<T>:: prox_pair(function<Real(T)> h, function<T(T, Real)> proxh){
    this->h = h; this->proxh = proxh;
}

grad_pair<Vec> LS(Mat, Vec);
extern prox_pair<Vec> L1_NORM;
#endif