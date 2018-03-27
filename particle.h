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
        double *pos;    //Coordinates of the charge
        char name[3];
        double chargeDisp[3];   //Charge displacement vector
        double com[3];  //Center of mass

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
        static Particle** create_dummies(Particle **particles);
        static void write_coordinates(char name[], Particle **particles);
        static Particle** read_coordinates(std::string name, bool relative);
};

#endif