#ifndef BASE_H
#define BASE_H
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <float.h>
#include "constants.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>
//#include "Eigen/Dense"

class Base{
    public:
        static double xL;
        static double yL;
        static double zL;
        static double eCummulative;
        static int acceptedMoves;
        static int totalMoves;
        static double T;
        static double lB;
        static double wall;
        static void set_lB(){
            lB = EC*EC/(4 * PI * VP * 2 * 1e-10 * KB * T);
            //lB = 69.6094879234987;//EC*EC/(4 * PI * VP * 2 * 1e-10 * KB * T);
        }
};

#endif

