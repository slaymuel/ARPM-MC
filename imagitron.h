#include "base.h"
#include "particle.h"

namespace energy { namespace imagitron{
    extern double wallCharge;

    void initialize();
    double get_particle_energy(Particle **particles, Particle *p);
    double get_energy(Particle **particles);
    double wall_charge(double z);
} }
