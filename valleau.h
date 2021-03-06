#ifndef VALLEAU_H
#define VALLEAU_H

#include "base.h"
#include "particle.h"
#include "particles.h"
#include "direct.h"
#include "analysis.h"
#include "ewald3D.h"

namespace energy{ namespace valleau{
    extern Eigen::VectorXd chargeVector;
    extern Eigen::VectorXd ext;
    extern Eigen::VectorXd imgExt;
    extern Eigen::VectorXd pDensity;
    extern Eigen::VectorXd nDensity;
    extern int numOfSamples;
    extern double binWidth;
    void update_charge_vector(Particles &particles);
    void update_potential();
    double get_images(Particles &particles);
    double get_particle_images(Particles &particles, Particle &p);
    double get_energy(Particles &particles);
    double get_particle_energy(Particles &particles, Particle &p);
    double phiw(double z);
    void initialize();
} }

#endif
