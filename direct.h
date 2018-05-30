#ifndef DIRECT_H
#define DIRECT_H

#include "base.cpp"
#include "particle.h"

class Direct{
    public:
        Direct();
        double get_energy(Particle **particles);
        double get_energy(Particle **particles, Particle *p);
        double get_central(Particle **particles);
        double get_central(Particle **particles, Particle *p);
        double get_replicates(Particle **particles);
        template<typename T>
        double norm(T vec);
};

#endif