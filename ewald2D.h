#ifndef Ewald2D_H
#define Ewald2D_H

#include "base.h"
#include "particle.h"
#include <vector>
#include <complex>

namespace energy{ namespace ewald2D{

        extern std::vector< Eigen::Vector3d > kVec;
        extern double *kNorm;
        extern int kNum;
        extern double alpha;

        double get_energy(Particle **particles);
        double get_particle_energy(Particle **particles, Particle *p);
        void initialize();
        void set_alpha();

        template<typename T>
        T erfc_x( T x );

        template<typename T>
        T erf_x( T x );

        template<typename T>
        double norm(T x);
        double get_self_correction(Particle *p);
        double f(double norm, double zDist);
        double g(Particle *p1, Particle *p2);
        double get_reciprocal(Particle *p1, Particle *p2);
        double get_real(Particle *p1, Particle *p2);
        double dipole_correction(Particle *p);

} }

#endif
