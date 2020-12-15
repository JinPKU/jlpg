
# todo : change the make file
% : %.cpp
	@g++  -O3 -march=native -std=c++11  $< -o $@  -Iinclude -DVERBOSED=1




.PHONY: clean
clean:
	@-rm ./test1
