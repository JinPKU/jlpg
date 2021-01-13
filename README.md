# A Software Package of Proximal Gradient Descent Method


Implemented by \
**Zeyu Jin**, PKU\
**Ting Lin**, PKU

## Structure of `/code`
`/include/` contains header file of our jlpg package, including
  - `jlpg.hpp` the wrapper of our all header file.
  - `funcpairs.hpp` define useful function pairs like `LS`, `L1_NORM`, `L2_NORM`
  - `problem.hpp` the construction and basic function of the objective.
  - `solver.hpp` the solver of our proximal gradient method
  - `continuation.hpp` provide an interface to use continuation method accelerating our program.

`/example` contains some examples we create in order to reveal the power of our package. 
 - `ex1.cpp` the lasso problem of size 1000.
 - `ex2.cpp` least square under a norm ball constraint. 

Inside `/code` dictionary, there is also a README file, introducing the basic framework of our package and some advanced parameters/settings. 


## Dependency
  1. Need C++11 support, **G++(>=5.0)** is recommended. It is welcome to inform us that the performance under other compile enviroments.
  2. **Eigen(>=3.0)** is required.  

## Basic assumption
We assume the following typedef
```C++
typedef Real double
typedef Vec Eigen::VectorXd;
typedef Mat Eigen::MatrixXd;
```
as our default setting. These lines locate in the beginning of `./code/include/funcpairs.hpp`. Feel free to change it!

## Mini example to illustrate how our package works:
1. Include the necessary files 
  ```C++
  #include <iostream>  
  // for IO
  #include <eigen3/Eigen/Dense> 
  // for using eigen
  #include "jlpg.hpp"
  // our package, change to the right file
  using namespace std;
  ```

2. Create Data for $A$ and $b$.
```C++
  Mat A(3,3); // Mat is MatrixXd in eigen.
  A << 1,2,3,4,5,6,7,8,9; // Assign value
  Vec b(3);  // Vec is VectorXd in eigen
  b << 1,4,9; // Assign Value
```

3. Set up the problem 
 $$\frac{1}{2}\|Ax-b\|_2^2+0.01\|x\|_1$$
```C++
Problem<Vec> p(LS(A,b),L1_NORM) // set up the problem with LeastSqaure and L1 norm.
p.mu = 0.01; // set up the coefficient of problem
```

4. Set up options for solver.
```C++
Options opts(10000, 1e-8, 1e-6, 1e-0, 5e-1); // Create the option
opts.setClassical();  // set strategy for the option: classical
```

5. Solve it happily
```C++
Outputs out; //output structure
Vec x(3); x << 0,0,0; //init value
x = pgm(p,x,opts,out);
cout << x << endl;
```

6. Suitable compile command
```
g++ -O3 -march=native -std=c++11 naivelasso -o naivelasso.out -I../include -DVERBOSED=1 
```
Here `-DVERBOSED=1` enable us to get the information at each iteration.



## What we support 
1. Support Least Square(both vector and matrix version) `LS(A,b)` and logistic regression `LOGISTIC(A,b)`. See doc for further reference.
2. Support the following proximal pair
  - `L1_NORM`, `L2_NORM`, `Linf_NORM` and `L0_NORM`
  - `L1_BALL(R)`, `L2_NORM(R)`, `Linf_NORM(R)` and `L0_NORM(R)`
  - Matrix norm in generalized LASSO problem, currently including `L12_NORM` only.
  - Spectral-relevant norm, currently including `NUCLEAR_NORM`  only.
  - `ELASTIC_NET(lam)` and `LOG_SUM`.
  - Naive gradient method, use `NO_PROX` or `NO_PROX_MAT`.

3. Support four kinds of optimization strategies to choose the step size:
  - `Constant`: constant step size.
  - `Armijo`: backtracking line search to achieve the Armijo condition.
  - `Nonmonotone`: backtracking line search to achieve a non-monotone condition using the BB step size.
  - `Classical`: a classical strategy to choose the step size in the proximal gradient method (See the document for details).


4. Easy construction of new `grad_pair` and `prox_pair`. 
```C++
Real h(Vec x){
  return L1_NORM.h(x.segment<3>(0))+L2_NORM.h(x.segment<3>(3));
}
Vec h(Vec x, Real t){
  Vec u = x;
  u.segment<3>(0) = L1_NORM.proxh(x.segment<3>(0),t);
  u.segment<3>(3) = L2_NORM.proxh(x.segment<3>(0),t);
  return u;
}
prox_pair<Vec> my_block_prox(h,proxh);

// new problem setting....
```
An automatic setting for block proximal pair is under constructed...

5. Support a continuation strategy.

