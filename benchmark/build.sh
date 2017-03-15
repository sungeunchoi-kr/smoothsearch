#!/bin/bash
#g++ -O3 main.cpp -orun -L"../lib/ntl-10.3.0/src" -lntl -lgmp -lm --std=c++11
#g++ -O0 main.cpp -orun -L/usr/local/lib -lntl -lgmp -lm --std=c++11
g++ -O0 main.cpp -orun -L/usr/local/lib -lntl -lpari -lgmp -lm --std=c++11

