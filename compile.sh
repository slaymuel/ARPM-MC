#!/bin/sh

gcc -c main.cpp particle.cpp base.cpp analysis.cpp ran2_lib.cpp mc.cpp ewald3D.cpp direct.cpp
gcc -o main main.o particle.o base.o analysis.o ran2.o ran2_lib.o mc.o ewald3D.o direct.o -lstdc++