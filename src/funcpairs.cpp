#include "funcpairs.hpp"

Real LS_F(Mat A, Vec b, Vec x){return .5*(A*x-b).squaredNorm();}
Vec LS_GRADF(Mat A, Vec b, Vec x){return A.transpose()*(A*x - b);}
grad_pair<Vec> LS(Mat A, Vec b){
    return grad_pair<Vec> (bind(LS_F,A,b,placeholders::_1),bind(LS_GRADF,A,b,placeholders::_1));
}
Real L1_NORM_H(Vec x){return x.lpNorm<1>();}
Vec L1_NORM_PROXH(Vec x, Real t){return (x.array().sign())*(x.array() - t).max(0);}
prox_pair<Vec> L1_NORM(L1_NORM_H,L1_NORM_PROXH);

