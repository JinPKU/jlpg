
# todo : change the make file
% : %.cpp
	@g++  -O3 -march=native -std=c++11  $< -o $@  -Iinclude




.PHONY: clean
clean:
	@-rm ./test1
