#ifndef MC_H
#define MC_H

#include "base.cpp"
#include "particle.h"
#include "math.h"
#include "ran2_lib.cpp"
#include "constants.h"

class MC: public Base{
    public:
        void equilibrate(Particle **particles);
        static int mcmove(Particle **particles, double dr);
        static double getParticleEnergy(int pInd, Particle *p, Particle **particles);
        double getEnergy(Particle **particles);
        void disperse(Particle **particles);
};

#endif