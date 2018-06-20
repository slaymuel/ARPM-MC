#ifndef VALLEAU_H
#define VALLEAU_H

#include "base.cpp"
#include "particle.h"
#include "direct.h"

namespace energy{ namespace valleau{
    extern Eigen::VectorXd chargeVector;
    extern Eigen::VectorXd ext;
    void update_charge_vector(Particle **particles);
    void update_potential();
    double get_energy(Particle **particles);
    double phiw(double z);
} }

#endif