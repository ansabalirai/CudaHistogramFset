all: program

program: cudacode.o
	g++ main.falcon.cpp kernel.falcon.o -o program -L/usr/local/cuda/lib64 -lcudart -lcuda


cudacode.o:
	nvcc -c kernel.falcon.cu 

clean: rm -rf *o program
