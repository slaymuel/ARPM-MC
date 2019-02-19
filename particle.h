#ifndef PARTICLE_H
#define PARTICLE_H

#include "base.h"
#include <fstream>
#include <sstream>

class Particle: public Base{
    public:
        Particle(bool dummie);
        Particle();

        //static double **distances;  //Distance matrix
        static int numOfParticles;
        static double oldEnergy;
        static int numOfElectrons;
        double b;   //Length of charge displacement
        double d;   //Diameter of particle
        double q;  //Charge
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

        
        double com_distance(Particle &p);
        double com_distance_xy(Particle &p);
        double distance_xy(Particle &p);
        double distance_z(Particle *p);
        bool wall_2d();

        static void place_particles(Particle **particles);
        static void create_electrons(std::vector<Particle> &particles, int num);
        static Particle** create_dummies(Particle **particles);


        static Particle** read_arpm_jan(std::string fileName);
        static Particle** read_coordinates(std::string name, bool relative, bool nanometers);
        //static Particle** read_coordinates_gro(std::string name);

        void pbc_pos();
        void pbc();
        static void pbc_xy(Eigen::Vector3d& x);
        double distance(Particle &p);
        void pbc(Eigen::Vector3d& x);     
};

#endif
