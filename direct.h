// #ifndef DIRECT_H
// #define DIRECT_H

// #include "base.cpp"
// #include "particle.h"

// class Direct{
//     public:
//         Direct();
//         static double get_energy(Particle **particles);
//         double get_energy(Particle **particles, Particle *p);
//         static double get_central(Particle **particles);
//         double get_central(Particle **particles, Particle *p);
//         static double get_replicates(Particle **particles);
//         double get_replicates(Particle **particles, Particle *p);
//         template<typename T>
//         static double norm(T vec);
// };

// #endif

#ifndef DIRECT_H
#define DIRECT_H

#include "base.h"
#include "particle.h"
#include "particles.h"

namespace energy{ namespace direct{ 
    double get_energy(Particles &particles);
    double get_particle_energy(Particles &particles, Particle &p);
    double get_particle_pot(Particles &particles, Particle &p);
    //double get_energy(Particle **particles, Particle *p);
    double get_central(Particles &particles);
    double get_central(Particles &particles, Particle &p);
    double get_central_pot(Particles &particles, Particle &p);
    double get_replicates(Particles &particles);
    double get_replicates(Particles &particles, Particle &p);
    double get_replicates_pot(Particles &particles, Particle &p);
    double phiw(double z);
    template<typename T>
    double norm(T vec);
} }
#endif
