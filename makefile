

% : %.cpp
	g++  -O3 -march=native -std=c++11  $< ./src/problem.cpp  ./src/funcpairs.cpp -o $@  -Iinclude