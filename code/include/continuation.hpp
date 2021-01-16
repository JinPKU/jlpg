#ifndef __CONTINUATION_HPP__
#define __CONTINUATION_HPP__


#include "solver.hpp"
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

struct ContOptions{
	const Options &opts;
	int maxiter_out;
	Real factor, mu1, ftol_init, gtol_init, etaf, etag;

	ContOptions(const Options &opts, int maxiter_out, Real factor, Real mu1, Real ftol_init, Real gtol_init, Real etaf, Real etag):
		opts(opts), maxiter_out(maxiter_out), factor(factor), mu1(mu1), ftol_init(ftol_init), gtol_init(gtol_init), etaf(etaf), etag(etag) {};
};

template <typename T>
T cont_pgm(Problem<T> p, T x, ContOptions conopts, Outputs &conoutput) {
	int iter_out;
	Options opts_tmp = conopts.opts;
	opts_tmp.ftol = opts_tmp.ftol * conopts.ftol_init;
	opts_tmp.gtol = opts_tmp.gtol * conopts.gtol_init;
	Real mu_tmp = conopts.mu1;
	Real F_cur = p.value(x);
	Real F_prev;

	Problem<T> p_tmp = p;
	p_tmp.mu = mu_tmp;
	Outputs output;
	clock_t Tstart, Tend;
	Tstart = clock();

	conoutput.Flag = false;
	conoutput.iter = 0;

	for ( iter_out = 0; iter_out < conopts.maxiter_out; ++iter_out ) {
		F_prev = F_cur;

		x = pgm(p_tmp, x, opts_tmp, output);
		conoutput.iter += output.iter;
		F_cur = p.value(x);
		

		if (mu_tmp == p.mu && (abs(F_cur - F_prev) < conopts.opts.ftol || output.nrmG < conopts.opts.gtol)) {
			conoutput.Flag = true;
			break;
		}

		if (output.Flag) {
			mu_tmp = max(mu_tmp * conopts.factor, p.mu);
			p_tmp.mu = mu_tmp;
		}

		opts_tmp.ftol = max(opts_tmp.ftol * conopts.etaf, conopts.opts.ftol);
		opts_tmp.gtol = max(opts_tmp.gtol * conopts.etag, conopts.opts.gtol);
	}

	Tend = clock();
	Real during = (double)(Tend - Tstart)/CLOCKS_PER_SEC;

	conoutput.cputime = during;
	conoutput.F_cur = p.value(x);
	conoutput.f_cur = p.f(x);
	conoutput.nrmG = (x - p.proxh(x - p.gradf(x), p.mu)).norm();

	return x;
}


#endif
