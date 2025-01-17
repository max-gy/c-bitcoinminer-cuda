# this is a btc miner for nvidia xavier

based on bitcoin miner from MattyAB https://github.com/MattyAB

but modified to run with CUDO on a nvidia xavier

for compiling I used nvcc like so

nvcc -o main.out main.cpp miner.cpp util.cpp mainCuda.cu
