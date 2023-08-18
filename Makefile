all: runCumulant

runCumulant: 
	g++ -std=c++11 -o runCumulant Cumulant.cpp ECorr.cpp `root-config --libs --cflags`

cbwc: 
	g++ -std=c++11 -o cbwc CBWC.cpp `root-config --libs --cflags`

duoCBWC: 
	g++ -std=c++11 -o duoCBWC duoCBWC.cpp `root-config --libs --cflags`
	

