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

`/example` contains some examples we create in order to reveal the power of our package. 
 - `ex1.cpp` the lasso problem of size 1000.
 - `ex2.cpp` least square under a norm ball constraint. 



## Dependency
  1. Need C++11 support, GCC(>=5.0) is recommended. It is welcome to inform us that the performance under other compile enviroments.
  2. Eigen(>=3.0) is required.  


## mini example to illustrate how our package works:
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