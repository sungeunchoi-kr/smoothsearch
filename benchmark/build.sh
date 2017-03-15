#!/bin/bash
g++ -O3 main.cpp -orun -L/usr/local/lib -lntl -lpari -lgmp -lm --std=c++11
