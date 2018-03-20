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
        double lB;

};

#endif

