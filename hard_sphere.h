#ifndef HS_H
#define HS_H
#include "particle.h"

namespace energy{ namespace hs{
    double get_energy(Particle **particles);
    double get_particle_energy(Particle **particles, Particle *p);
}}

#endif
