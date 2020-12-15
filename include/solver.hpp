#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#ifndef __PROBLEM_HPP__
#include "problem.hpp"
#endif

#include <cmath>


template <typename T>
T pgm(Problem<T> p, T x,int maxiter, Real tol){
    Real f_prev = 1e8;
    for(int iter = 0; iter < maxiter; iter ++){
        Real alpha = 0.01;
        x = p.proxh(x - alpha*p.gradf(x), alpha*p.mu);
        Real f_cur = p.f(x);
        if(abs(f_cur - f_prev) < tol){ break;} 
    }
    return x;
};

template <typename T>
T pgm_armijo(Problem<T>p, T x, int maxiter, Real tol, Real rho, int bt_time){
// todo

    return x;
}

template <typename T>
T pgm_nbb(Problem<T>p, T x, int maxiter, Real tol){

    return x;
}
#endif