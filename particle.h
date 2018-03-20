#ifndef PARTICLE_H
#define PARTICLE_H

#include "base.cpp"
#include <fstream>
#include <sstream>

class Particle: public Base{
    public:
        double b;
        double d;
        int q;
        int index;
        double *pos;
        char name[3];

        static int numOfParticles;
        Particle();

        void pbc();
        void pbc_xy();
        void randomMove(double stepSize);
        void randomMove_xy(double stepSize);
        double distance(Particle *p);
        double distance_xy(Particle *p);
        double distance_z(Particle *p);
        int hardSphere(Particle **particles);

        static int get_overlaps(Particle ** particles);
        static void place_particles(Particle **particles);
        static Particle** create_particles(int num);
        static void write_coordinates(char name[], Particle **particles);
        static Particle** read_coordinates(char *name);
};

#endif