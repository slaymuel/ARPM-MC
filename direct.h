#ifndef DIRECT_H
#define DIRECT_H

#include "base.cpp"
#include "particle.h"

class Direct{
    public:
        Direct();
        double get_energy(Particle *p, Particle **particles);
        template<typename T>
        double norm(T vec);
};

#endif