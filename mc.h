#ifndef MC_H
#define MC_H

#include "base.cpp"
#include "particle.h"
#include "ran2_lib.cpp"
#include "ewald3D.h"
#include "ewald2D.h"
#include "direct.h"
//#include "levin.cpp"

class MC: public Base{
    public:
        void equilibrate(Particle **particles);
        static int trans_move(Particle **particles, double dr);
        static int charge_rot_move(Particle **particles);
        static double get_particle_energy(int pInd, Particle *p, Particle **particles);
        double get_energy(Particle **particles);
        void disperse(Particle **particles);
        static Ewald3D ewald3D;
        static Ewald2D ewald2D;
        static Direct direct;
        //static Levin levin;
};

#endif