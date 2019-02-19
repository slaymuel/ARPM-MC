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

    double get_energy(Particles &particles);
    double get_particle_energy(Particles &particles, Particle &p);
    void initialize(Particles &particles);
    void update_f(Particle &_old, Particle &_new);
    double u_gamma(Particles &particles);
    double get_polarization();

    extern int kNumMax;
    extern double gamma;
    //std::vector< std::vector<double> > kVec;
    //std::complex<double> *rkVec;
    extern int kNum;
    extern int kMax;
    extern std::vector< std::vector<double> > kVectors;
    extern std::vector<double> kNorms;
    extern std::vector<double> f1;
    extern std::vector<double> f2;
    extern std::vector<double> f3;
    extern std::vector<double> f4;
    extern std::vector<double> eFactors;
    extern double uGamma;

    
} }

#endif
