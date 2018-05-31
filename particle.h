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
        //double *pos;    //Coordinates of the charge
        Eigen::Vector3d pos;
        Eigen::Vector3d chargeDisp;
        Eigen::Vector3d com;
        char name[3];
        //double chargeDisp[3];   //Charge displacement vector
        //double com[3];  //Center of mass
        static double **distances;
        static int numOfParticles;
        Particle(bool dummie);
        Particle();
        void pbc_xy();
        void random_move(double stepSize);
        void random_charge_rot();
        void randomMove_xy(double stepSize);
        double distance(Particle *p);
        double com_distance(Particle *p);
        double distance_xy(Particle *p);
        double distance_z(Particle *p);
        int hard_sphere(Particle **particles);

        static void update_distances(Particle **particles);
        static void update_distances(Particle **particles, Particle *p);
        static int get_overlaps(Particle ** particles);
        static void place_particles(Particle **particles);
        static Particle** create_particles(int num);
        static Particle** create_dummies(Particle **particles);
        static void write_coordinates(char name[], Particle **particles);
        static void write_charge_coordinates(char name[], Particle **particles);
        static Particle** read_jan(std::string pName, std::string nName);
        static Particle** read_arpm_jan(std::string fileName);
        static Particle** read_coordinates(std::string name, bool relative, bool nanometers);

    private:
        void pbc(Eigen::Vector3d& x);
        void pbc();
};

#endif