#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#ifndef __PROBLEM_HPP__
#include "problem.hpp"
#endif

#include <cmath>

struct Options{
	bool ls, bb;
	int maxiter;
	Real ftol, gtol;
	Real alpha;

	Options(bool ls_, bool bb_, int maxiter_, Real ftol_, Real gtol_, Real alpha_): ls(ls_), bb(bb_), maxiter(maxiter_), ftol(ftol_), gtol(gtol_), alpha(alpha_) {};
};

template <typename T>
T pgm(Problem<T> p, T x, Options opts){
#if VERBOSED
    cout << "Proximal Gradient Method " << endl;
#endif
	Real F_prev, f_prev;	// F = f + h
	T x_prev, g_prev;
	Real F_cur = p.value(x);
    Real f_cur = p.f(x);
	T g_cur = p.gradf(x);
	Real alpha = opts.alpha;
	T dx, dg;
	Real dxg;

    int iter, nls;
	Real tmp, nrmG;

    for( iter = 0; iter < opts.maxiter; ++iter ){
		F_prev = F_cur;
		f_prev = f_cur;
		g_prev = g_cur;
		x_prev = x;

        x = p.proxh(x_prev - alpha * g_prev, alpha * p.mu);
		if (opts.ls) {
			nls = 0;
			while (true) {
				tmp = p.f(x);
				if (tmp <= f_prev + g_prev.dot(x - x_prev) + .5 / alpha * (x - x_prev).squaredNorm() || nls == 10) { break; }
				
				alpha = 0.5 * alpha; nls = nls + 1;
				x = p.proxh(x_prev - alpha * g_prev, alpha * p.mu);
			}

			f_cur = tmp;
			F_cur = p.value(x);
		}
		else {
			f_cur = p.f(x);
			F_cur = p.value(x);
		}
		
#if VERBOSED
    cout << "In iteration " << iter << ", objective function value = "<< F_cur << endl;
#endif
		
		g_cur = p.gradf(x);
		nrmG = (x - p.proxh(x - g_cur, p.mu)).norm();

		if (opts.bb && opts.ls) {
			dx = x - x_prev;
			dg = g_cur - g_prev;
			dxg = abs(dx.dot(dg));

			if (dxg > 0) {
				if (iter & 1) {
					alpha = dxg / dg.squaredNorm();
				}
				else {
					alpha = dx.squaredNorm() / dxg;
				} 
			}
			alpha = min(max(alpha, 1e-4), 1e4);
		} 
		else {
			alpha = opts.alpha;
		}

        if (abs(F_cur - F_prev) < opts.ftol || nrmG < opts.gtol){
			break;
		} 
    }
#if VERBOSED
    cout << "Problem Solved within " << iter+1 << " Iteration(s)." << endl;
#endif
    return x;
};


#endif
