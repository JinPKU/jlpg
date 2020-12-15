# A Software Package of Proximal Gradient Descent Method


Implemented by \
**Zeyu Jin**, PKU\
**Ting Lin**, PKU


`\include` contains header file of our jlpg package, including
  - `jlpg.hpp` the wrapper of our all header file.
  - `funcpairs.hpp` define useful function pairs like `LS`, `L1_NORM`, `L2_NORM`
  - `problem.hpp` the construction and basic function of the objective.
  - `solver.hpp` the solver of our proximal gradient method

Things todo:
1. Armijo and Nonmonotone BB search 
2. various gradpair and proxpair
3. examples
