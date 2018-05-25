#!/bin/sh

gcc -c main.cpp particle.cpp base.cpp levin.cpp analysis.cpp ran2_lib.cpp mc.cpp ewald3D.cpp direct.cpp ewald2D.cpp -Xpreprocessor -fopenmp -O3
gcc -o main main.o particle.o base.o levin.o analysis.o ran2.o ran2_lib.o mc.o ewald3D.o direct.o ewald2D.o -lstdc++ -lboost_program_options -lomp