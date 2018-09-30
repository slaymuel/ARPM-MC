#!/bin/sh

time icpc -c main.cpp img_rep.cpp imagitron.cpp particle.cpp levin.cpp analysis.cpp ran2_lib.cpp mc.cpp ewald3D.cpp direct.cpp ewald2D.cpp valleau.cpp hard_sphere.cpp -O3 -march=native -pipe -std=c++11 -ipo
icpc -o main main.o img_rep.o imagitron.o particle.o levin.o analysis.o ran2.o ran2_lib.o mc.o ewald3D.o direct.o ewald2D.o valleau.o hard_sphere.o -O3 -lstdc++ -lboost_program_options