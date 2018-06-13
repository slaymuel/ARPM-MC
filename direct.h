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

#include "base.cpp"
#include "particle.h"

namespace energy{ namespace direct{ 
    double get_energy(Particle **particles);
    //double get_energy(Particle **particles, Particle *p);
    double get_central(Particle **particles);
    double get_central(Particle **particles, Particle *p);
    double get_replicates(Particle **particles);
    double get_replicates(Particle **particles, Particle *p);
    template<typename T>
    double norm(T vec);
} }
#endif