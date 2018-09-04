#ifndef BASE_H
#define BASE_H
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include <float.h>
#include "constants.h"
#include <Eigen/Dense>
#include <vector>
//#include <omp.h>

//#include "Eigen/Dense"

class Base{
    public:
        static double xL;
        static double yL;
        static double zL;
        static double eCummulative;
        static double T;
        static double lB;
        static double wall;
        static double beta;
        static double P;
        static double volume;
        static int acceptedMoves;
        static int totalMoves;
        static std::vector<double> volumes;
        static bool d2;

        static void set_beta(){
            beta = 1/(KB * T);
        }

        static void set_lB(){
            lB = EC * EC / (4.0 * PI * VP * 2.0 * 1e-10 * KB * T);
            //lB = 69.6094879234987;//EC*EC/(4 * PI * VP * 2.0 * 1e-10 * KB * T);
        }
};

#endif