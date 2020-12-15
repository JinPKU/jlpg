#include <functional>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "funcpairs.cpp"


template <typename T>
class Problem{
    public:
    Problem(){}
    Problem(function<Real(T)> f, function<T(T)> gradf, function<Real(T)> h, function<T(T,Real)> proxh);
    Problem(grad_pair<T> fpair, prox_pair<T> hpair);
    Problem(grad_pair<T> fpair);
    Problem(prox_pair<T> hpair);


    function<Real(T)> f;
    function<T(T)> gradf;
    function<Real(T)> h;
    function<T(T,Real)> proxh;
    Real mu;
    Real value(T x){
        return (this->f(x)) + (this->mu) * (this->h(x));
    }

};

template <typename T>
Problem<T>::Problem(function<Real(T)> f, function<T(T)> gradf, function<Real(T)> h, function<T(T,Real)> proxh){
    this->f = f; this->gradf = gradf; this->h = h; this->proxh = proxh;
}

template <typename T>
Problem<T>::Problem(grad_pair<T> fpair, prox_pair<T> hpair){
    this->f = fpair.f; this->gradf = fpair.gradf;
    this->h = hpair.h; this->proxh = hpair.proxh;
}

template <typename T>
Problem<T>::Problem(grad_pair<T> fpair){
    this->f = fpair.f; this->gradf = fpair.gradf;
}

template <typename T>
Problem<T>::Problem(prox_pair<T> hpair){
    this->h = hpair.h; this->proxh = hpair.proxh;
}