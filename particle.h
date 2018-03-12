#ifndef PARTICLE_H
#define PARTICLE_H

#include "base.cpp"

class Particle: public Base{
    public:
        double b;
        double d;
        int q;
        int index;
        double *pos;
        char name[3];

        static int numOfParticles;
    //    Particle();
    //    ~Particle();
        Particle();

        void pbc();
        void pbc_xy();
        void randomMove();
        double distance(Particle *p);
        double distance_xy(Particle *p);
        int hardSphere(Particle **particles);
};

#endif