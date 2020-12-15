#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#ifndef __PROBLEM_HPP__
#include "problem.hpp"
#endif

template <typename T> 
class StepsizeScheduler{
    public:
    StepsizeScheduler(){}
};
template <typename T>
T pgm(Problem<T> p, T x,int maxiter, double tol){
    return x;
};

#endif