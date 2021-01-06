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


Real LS_F(Mat A, Vec b, Vec x){return .5*(A*x-b).squaredNorm();}
Vec LS_GRADF(Mat A, Vec b, Vec x){return A.transpose()*(A*x - b);}
grad_pair<Vec> LS(Mat A, Vec b){
    return grad_pair<Vec> (bind(LS_F,A,b,placeholders::_1),bind(LS_GRADF,A,b,placeholders::_1));
}
Real L1_NORM_H(Vec x){return x.lpNorm<1>();}
Vec L1_NORM_PROXH(Vec x, Real t){return (x.array().sign())*(x.array() - t).max(0);}
prox_pair<Vec> L1_NORM(L1_NORM_H,L1_NORM_PROXH);




#endif    // FUNCPAIRS_HPP

