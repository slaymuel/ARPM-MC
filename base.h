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
        static double xLHalf;
        static double yL;
        static double yLHalf;
        static double zL;
        static double zLBox;
        static double zLBoxHalf;
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
        static std::vector<double> box;

        static void set_beta(){
            beta = 1.0 / (KB * T);
        }

        static void set_lB(){
            lB = EC * EC / (4.0 * PI * VP * 78.3 * 1e-10 * KB * T); //7.16 // 7.1599645637319043;
            //lB = 69.6094879234987;//EC*EC/(4 * PI * VP * 2.0 * 1e-10 * KB * T);
        }
};

#endif
