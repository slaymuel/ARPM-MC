#ifndef VALLEAU_H
#define VALLEAU_H

#include "base.h"
#include "particle.h"
#include "direct.h"
#include "analysis.h"

namespace energy{ namespace valleau{
    extern Eigen::VectorXd chargeVector;
    extern Eigen::VectorXd ext;
    extern Eigen::VectorXd imgExt;
    extern Eigen::VectorXd pDensity;
    extern Eigen::VectorXd nDensity;
    extern int numOfSamples;
    extern double binWidth;
    void update_charge_vector(Particle **particles);
    void update_potential();
    double get_images(Particle **particles);
    double get_particle_images(Particle **particles, Particle *p);
    double get_energy(Particle **particles);
    double get_particle_energy(Particle **particles, Particle *p);
    double phiw(double z);
    void initialize();
} }

#endif