#ifndef VALLEAU_H
#define VALLEAU_H

#include "base.cpp"
#include "particle.h"

namespace energy{ namespace valleau{
    extern Eigen::VectorXd chargeVector;
    extern Eigen::VectorXd ext;
    void update_charge_vector(Particle **particles);
    void calculate_potential();
    double phiw(double z);
} }

#endif