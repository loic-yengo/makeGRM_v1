makeGRM-GPD: 	makeGRM-GPD.cpp
	g++ -fopenmp -I"./eigen-3.4-rc1/" -Wall -std=c++11 -g -O2 makeGRM-GPD.cpp -o makeGRM-GPD -lm
clean:
	rm makeGRM-GPD

