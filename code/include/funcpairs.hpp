// func pairs provide an interface for problem, such as  LS, norm ball of norm 
// with easy construction method, the classes provide a simple but general/useful abstractions for proximal gradient methods.


/* ------------------------------------
        FUNCPAIRS.HPP
-------------------------------------*/

#ifndef __FUNCPAIRS_HPP__
#define __FUNCPAIRS_HPP__
#include <functional>
#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace std;

/* type definition, you can change it if speicial requirement is needed. */
typedef double Real;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;
double inf = std::numeric_limits<double>::infinity();

/* ------------------------------------------------------------------------- 
        CLASS DECLARATIONS
------------------------------------------------------------------------- */

template <typename T>
class grad_pair{
    /*
    Class for smooth function. Containing 
        - function<Real(T)> f : the value of function.
        - function<T(T)> gradf: the gradient of function. 
    */
    public:

    /* construction functions */
    grad_pair(){}
    grad_pair(function<Real(T)> f, function<T(T)> gradf);

    /* members */ 

    function<Real(T)> f;
    function<T(T)> gradf;
};

template <typename T>
class prox_pair{
    /*
        Class for convex function. Containing 
            - function<Real(T)> h: the value of function. 
            - function<T(T,Real)> proxh: the proximal operator of function. proxh(x,t) returns the scaled proximal operator
    */

    public:

    /* construction functions */
    prox_pair(){}
    prox_pair(function<Real(T)> h, function<T(T, Real)> proxh);

    /* members */
    function<Real(T)> h;
    function<T(T,Real)> proxh;
};



/* -----------------------------
        CLASS DEFINITION
----------------------------  */

template<typename T>
grad_pair<T>:: grad_pair(function<Real(T)> f, function<T(T)> gradf){
    this->f = f; this->gradf = gradf;
}


template<typename T>
prox_pair<T>:: prox_pair(function<Real(T)> h, function<T(T, Real)> proxh){
    this->h = h; this->proxh = proxh;
}

/*--------------------------------
        BUILT-IN IMPLEMENTATIONS
------------------------------- */
// No grad pair


// No prox pair
Real NO_PROX_H(Mat x){return 0;}
Mat NO_PROX_PROXH(Mat x, Real t){return x;}
prox_pair<Mat> NO_PROX(NO_PROX_H,NO_PROX_PROXH);


// vec LS
Real LS_F(Mat A, Vec b, Vec x){return .5*(A*x-b).squaredNorm();}
Vec LS_GRADF(Mat A, Vec b, Vec x){return A.transpose()*(A*x - b);}
grad_pair<Vec> LS(Mat A, Vec b){
    return grad_pair<Vec> (bind(LS_F,A,b,placeholders::_1),bind(LS_GRADF,A,b,placeholders::_1));
}

// mat LS 
Real LS_FM(Mat A, Mat b, Mat x){return .5*(A*x-b).squaredNorm();}
Mat LS_GRADFM(Mat A, Mat b, Mat x){return A.transpose()*(A*x - b);}
grad_pair<Mat> LS(Mat A, Mat b){
    return grad_pair<Mat> (bind(LS_FM,A,b,placeholders::_1),bind(LS_GRADFM,A,b,placeholders::_1));
}


// logistic regression
Real LOGISTIC_F(Mat A, Vec b, Vec x){
return ((-b.array()*(A.transpose()*x).array()).exp()+1).log().mean();
}
Vec LOGISTIC_GRADF(Mat A, Vec b, Vec x){
    int m = A.cols();
    Vec xx = b.array()/((b.array()*(A.transpose()*x).array()).exp()+1);
    return -A*xx.matrix()/m;
    // return -A*(b.array()/((-b.array()*(A.transpose()*x).array()).exp()+1)/A.cols()).matrix();
}
grad_pair<Vec> LOGISTIC(Mat A, Vec b){
    return grad_pair<Vec> (bind(LOGISTIC_F,A,b,placeholders::_1),bind(LOGISTIC_GRADF,A,b,placeholders::_1));
}


// L0 form
Real L0_NORM_H(Vec x){return (Real)(x.nonZeros());}
Vec L0_NORM_PROXH(Vec x, Real t){
    Real s = sqrt(2*t);
    int sz = x.size();
    Vec u(x);
    for(int i=0;i<sz;i++){u[i] = abs(x[i])>s?x[i]:0;}
    return u;
}
prox_pair<Vec> L0_NORM(L0_NORM_H,L0_NORM_PROXH);



// L1 norm 
Real L1_NORM_H(Vec x){return x.lpNorm<1>();}
Vec L1_NORM_PROXH(Vec x, Real t){return (x.array().sign())*(x.array() - t).max(0);}
prox_pair<Vec> L1_NORM(L1_NORM_H,L1_NORM_PROXH);



// L2 Norm
Real L2_NORM_H(Vec x){return x.lpNorm<2>();}
Vec L2_NORM_PROXH(Vec x, Real t){
Real nrm = x.lpNorm<2>();
return max(nrm-t,0.0)/nrm*x;
}
prox_pair<Vec> L2_NORM(L2_NORM_H,L2_NORM_PROXH);



// L inf norm 




// L12 norm 
Real L12_NORM_H(Mat x){
    return x.rowwise().lpNorm<2>().sum();
}
Mat L12_NORM_PROXH(Mat x, Real t){
    Vec u = x.rowwise().lpNorm<2>();

    return ((x.array().colwise()) * (1-t/u.array()).max(0)).matrix();
}

prox_pair<Mat> L12_NORM(L12_NORM_H,L12_NORM_PROXH);


// L21 norm






// nuclear norm 
Real NUCLEAR_NORM_H(Mat x){
    Eigen::BDCSVD<Mat> svd(x, Eigen::ComputeThinU|Eigen::ComputeThinV); // this might lead to inefficiency.
    return svd.singularValues().sum();
}

Mat NUCLEAR_NORM_PROXH(Mat x, Real t){
        Eigen::BDCSVD<Mat> svd(x, Eigen::ComputeThinU|Eigen::ComputeThinV);
        Vec s = svd.singularValues();
        return svd.matrixU() * (s.array() - t).max(0).matrix().asDiagonal() * svd.matrixV().transpose();
}

prox_pair<Mat> NUCLEAR_NORM(NUCLEAR_NORM_H,NUCLEAR_NORM_PROXH);



// elastic net 
Real ELASTIC_NET_H(Vec x, Real lam){
    return x.lpNorm<1>() + lam/2*x.squaredNorm();
}
Vec ELASTIC_NET_PROXH(Vec x, Real lam, Real t){return (x.array().sign())*(x.array() - t).max(0)/(1+lam*t);}

prox_pair<Vec> ELASTIC_NET(Real lam){
return prox_pair<Vec> (bind(ELASTIC_NET_H,placeholders::_1, lam), bind(ELASTIC_NET_PROXH, placeholders::_1, lam, placeholders::_2));
}

// sum log
Real SUM_LOG_H(Vec x){
    return -x.array().log().sum();
}

Vec SUM_LOG_PROXH(Vec x, Real t){return (x.array()+(x.array()*x.array()+4*t).sqrt())/2;}

prox_pair<Vec> SUM_LOG(SUM_LOG_H, SUM_LOG_PROXH);

// L0 ball 
Real L0_BALL_H(Vec x, Real R){
return x.nonZeros()<R?0:inf;
}

Vec L0_BALL_PROXH(Vec x, Real R, Real t){
    Vec u = x.array().abs();
    int sz = x.size();
    std::sort(u.data(),u.data()+u.size(),std::greater<Real>());
    int R0 = min((int)floor(R),sz);
    double th = R0>0?u[R0-1]:0;
    for(int i = 0 ; i <sz; i++){u[i] = abs(x[i])>=th?abs(x[i]):0;}
    return u;
}

prox_pair<Vec> L0_BALL(Real R){
    return prox_pair<Vec>(bind(L0_BALL_H, placeholders::_1, R), bind(L0_BALL_PROXH, placeholders::_1, R, placeholders::_2));
}


// L1 ball 




// L2 ball 
Real L2_BALL_H(Vec x, Real R){
return x.lpNorm<2>()<R?0:inf;
}

Vec L2_BALL_PROXH(Vec x, Real R, Real t){
    return R/max(x.lpNorm<2>(),R)*x;
}

prox_pair<Vec> L2_BALL(Real R){
    return prox_pair<Vec>(bind(L2_BALL_H, placeholders::_1, R), bind(L2_BALL_PROXH, placeholders::_1, R, placeholders::_2));
}



// Linf ball 



// simple box
Real BOX_H(Vec x, Vec lx, Vec ux){
    return (x.array()>lx.array()).all() && (x.array()<ux.array()).all()?0:inf;
}
Vec BOX_PROXH(Vec x, Vec lx, Vec ux){
    return x.cwiseMax(lx).cwiseMin(ux);
}
prox_pair<Vec> BOX(Vec lx, Vec ux){
    return prox_pair<Vec> (bind(BOX_H,placeholders::_1,lx,ux),bind(BOX_PROXH,placeholders::_1,lx,ux));
}
#endif    /* __JLPG_HPP__ */

