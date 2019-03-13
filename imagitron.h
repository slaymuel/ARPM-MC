#ifndef IMAGITRON_H
#define IMAGITRON_H

#include "base.h"
#include "particle.h"
#include "particles.h"

namespace energy { namespace imagitron{
    extern double wallCharge;

    void initialize(Particles &particles);
    double get_particle_energy(Particles &particles, Particle &p);
    double get_energy(Particles &particles);
    double wall_charge(double z);
} }

#endif