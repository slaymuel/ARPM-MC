#ifndef PARTICLE_H
#define PARTICLE_H

#include "base.h"
#include <fstream>
#include <sstream>

class Particle: public Base{
    public:
        Particle(bool dummie);
        Particle();

        static double **distances;  //Distance matrix
        static int numOfParticles;
        static double oldEnergy;
        static int numOfElectrons;
        double b;   //Length of charge displacement
        double d;   //Diameter of particle
        int q;  //Charge
        int index;
        int belongsTo; //Belongs to molecule
        Eigen::Vector3d pos;    //Coordinates of the charge
        Eigen::Vector3d oldPos;
        Eigen::Vector3d chargeDisp; //Charge vector
        Eigen::Vector3d com;    //Center of mass coordinates
        Eigen::Vector3d oldCom;
        char name[3];

        void random_move(double stepSize);
        void random_charge_rot();
        void randomMove_xy(double stepSize);
        int hard_sphere(Particle **particles);
        
        double com_distance(Particle *p);
        double com_distance_xy(Particle *p);
        double distance_xy(Particle *p);
        double distance_z(Particle *p);
        static void update_distances(Particle **particles);
        static void update_distances(Particle **particles, Particle *p);
        static int get_overlaps(Particle ** particles);
        static void place_particles(Particle **particles);
        static Particle** create_particles(int nNum, int pNum, int eNum);
        static void create_electrons(Particle** particles, int num);
        static Particle** create_dummies(Particle **particles);
        static void write_coordinates(char name[40], Particle **particles);
        static void write_charge_coordinates(char name[], Particle **particles);
        static Particle** read_jan(std::string pName, std::string nName);
        static Particle** read_arpm_jan(std::string fileName);
        static Particle** read_coordinates(std::string name, bool relative, bool nanometers);
        static Particle** read_coordinates_gro(std::string name);
        void pbc_pos();
        void pbc();
    private:

        void pbc(Eigen::Vector3d& x);
        void pbc_xy(Eigen::Vector3d& x);
        double distance(Particle *p);
};

#endif