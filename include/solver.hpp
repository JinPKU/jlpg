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
T pgm_armijo(Problem<T>p, T x, int maxiter, Real tol, Real alpha){
#if VERBOSED
    cout << "Proximal Gradient Method by Armijo Backtracking Method" << endl;
#endif
	Real f_prev, f_cur, fs_cur, fs_prev, tmp;
	T g_prev, g_cur;
	T x_prev;

	fs_cur = p.f(x);
	f_cur = p.value(x);
	g_cur = p.gradf(x);

    int iter, nls;
    for( iter = 0; iter < maxiter; iter ++){
		f_prev = f_cur;
		fs_prev = fs_cur;
		g_prev = g_cur;
		x_prev = x;

        x = p.proxh(x - alpha * g_cur, alpha * p.mu);

		nls = 0;
		while (true) {
			tmp = p.value(x);
			if (tmp <= fs_prev + g_prev.dot(x - x_prev) + .5 / alpha * (x - x_prev).squaredNorm() || nls == 10) {
				break;
			}

			alpha = 0.5 * alpha;
			nls = nls + 1;
			x = p.proxh(x - alpha * g_cur, alpha * p.mu);
		}

		fs_cur = tmp;
		f_cur = p.value(x);
		g_cur = p.gradf(x);
	


#if VERBOSED
    cout << "In iteration " << iter << ", objective function value = "<< f_cur << endl;
#endif
        if(abs(f_cur - f_prev) < tol){ break;} 
    }
#if VERBOSED
    cout << "Problem Solved within " << iter+1 << " Iteration(s)." << endl;
#endif
    return x;
}

template <typename T>
T pgm_nbb(Problem<T>p, T x, int maxiter, Real tol){

    return x;
}
#endif
