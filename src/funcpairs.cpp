
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

template<typename T>
grad_pair<T>:: grad_pair(function<Real(T)> f, function<T(T)> gradf){
    this->f = f; this->gradf = gradf;
}

template <typename T>
class prox_pair{
    public:
    prox_pair(){}
    prox_pair(function<Real(T)> h, function<T(T, Real)> proxh);

    function<Real(T)> h;
    function<T(T,Real)> proxh;
};

template<typename T>
prox_pair<T>:: prox_pair(function<Real(T)> h, function<T(T, Real)> proxh){
    this->h = h; this->proxh = proxh;
}



Real LS_F(Mat A, Vec b, Vec x){return .5*(A*x-b).squaredNorm();}
Vec LS_GRADF(Mat A, Vec b, Vec x){return A.transpose()*(A*x - b);}
grad_pair<Vec> LS(Mat A, Vec b){
    return grad_pair<Vec> (bind(LS_F,A,b,placeholders::_1),bind(LS_GRADF,A,b,placeholders::_1));
}
Real L1_NORM_H(Vec x){return x.lpNorm<1>();}
Vec L1_NORM_PROXH(Vec x, Real t){return (x.array().sign())*(x.array() - t).max(0);}


prox_pair<Vec> L1_NORM(L1_NORM_H, L1_NORM_PROXH);

