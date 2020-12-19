
main : main
	@g++ -O3 -march=native -std=c++11 main.cpp -o main.out -Iinclude -DVERBOSED=1

