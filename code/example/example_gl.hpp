#ifndef __EXAMPLE_GL__
#define __EXAMPLE_GL__

#include <eigen3/Eigen/Dense>
#include "../include/jlpg.hpp"

void example_gl(Problem<Mat> &p) {
    Mat A(3,2);
    A<<1,2,3,4,5,6;
    
    Mat b(3,2); b<<1,1,1,1,1,1;
    p = Problem<Mat>(LS(A, b), L12_NORM);
    p.mu = 0.01;
}


#endif