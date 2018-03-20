#ifndef MC_H
#define MC_H

#include "base.cpp"
#include "particle.h"
#include "ran2_lib.cpp"
#include "ewald3D.h"
#include "direct.h"

class MC: public Base{
    public:
        void equilibrate(Particle **particles);
        static int mcmove(Particle **particles, double dr);
        static double getParticleEnergy(int pInd, Particle *p, Particle **particles);
        double getEnergy(Particle **particles);
        void disperse(Particle **particles);
        static Ewald3D ewald;
        static Direct direct;
};

#endif