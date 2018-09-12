#ifndef sp3_H
#define sp3_H

#include "particle.h"

class sp3{
    private:
        double get_short(Particle **particles);
        double split(double distance, double cutoff);
        double get_self(Particle **particles, double cutoff);
    public:
        double get_energy(Particle **particles);
};

#endif
