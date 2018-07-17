#ifndef Levin_H
#define Levin_H

#include "base.h"
//#include <eigen3/Eigen/Dense>
#include "particle.h"
#include "direct.h"
#include "valleau.h"
#include "ewald3D.h"
//#include <vector>
//#include <complex>
//#include <eigen3/Eigen/Dense>

namespace energy{ namespace levin{

    double get_energy(Particle **particles);
    double get_particle_energy(Particle **particles, Particle *p);
    void initialize(Particle **particles);
    void update_f(Particle *_old, Particle *_new);
    double u_gamma(Particle **particles);
    double get_polarization();

    extern int kNumMax;
    extern double gamma;
    //std::vector< std::vector<double> > kVec;
    //std::complex<double> *rkVec;
    extern int kNum;
    extern int kMax;
    extern Eigen::MatrixXd kVectors;
    extern Eigen::ArrayXd kNorms;
    extern Eigen::ArrayXd f1;
    extern Eigen::ArrayXd f2;
    extern Eigen::ArrayXd f3;
    extern Eigen::ArrayXd f4;
    extern Eigen::ArrayXd eFactors;
    extern double uGamma;

    
} }

#endif