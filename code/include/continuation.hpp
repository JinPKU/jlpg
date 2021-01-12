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
	Options opts;
	int maxiter_out;
	Real factor, mu1, ftol_init, gtol_init, etaf, etag;
};

template <typename T>
T cont_pgm(Problem<T> p, T x, ContOptions conopts) {
	int iter_out;
	Options opts_tmp = conopts.opts;
	opts_tmp.ftol = opts_tmp.ftol * conopts.ftol_init;
	opts_tmp.gtol = opts_tmp.gtol * conopts.gtol_init;
	Real mu_t = conopts.mu1;
	Real mu = p.mu;
	Real F_cur = p.value(x);
	Real F_prev;
	p.mu = mu_t;

	for ( iter_out = 0; iter_out < conopts.maxiter_out; ++iter_out ) {
		F_prev = F_cur;

		x = pgm(p, x, opts_tmp);

		if (mu_t == mu && (abs(F_cur - F_prev) < conopts.opts.ftol || out.nrmG < conopts.opts.gtol)) {
			break;
		}

		if (out.flag) {
			mu_t = max(mu_t * conopts.factor, mu);
		}

		opts_tmp.ftol = max(opts_tmp.ftol * conopts.etaf, conopts.opts.ftol);
		opts_tmp.gtol = max(opts_tmp.gtol * conopts.etag, conopts.opts.gtol);
	}
}


#endif
