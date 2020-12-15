#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#ifndef __PROBLEM_HPP__
#include "problem.hpp"
#endif

#include <cmath>


template <typename T>
T pgm(Problem<T> p, T x,int maxiter, Real tol, Real alpha){
#if VERBOSED
    cout << "Proximal Gradient Method with Constant Stepsize " <<  alpha << endl;
#endif
    Real f_prev = p.value(x);
    int iter;
    for( iter = 0; iter < maxiter; iter ++){
        x = p.proxh(x - alpha*p.gradf(x), alpha*p.mu);
        Real f_cur = p.value(x);
#if VERBOSED
    cout << "In iteration " << iter << ", objective function value = "<< f_cur << endl;
#endif
        if(abs(f_cur - f_prev) < tol){ break;} 
        f_prev = f_cur;
    }
#if VERBOSED
    cout << "Problem Solved within " << iter+1 << " Iteration(s)." << endl;
#endif
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