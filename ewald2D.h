#ifndef Ewald2D_H
#define Ewald2D_H

#include "base.h"
#include "particle.h"
#include <vector>
#include <complex>

class Ewald2D{
    public:
        Ewald2D();
        double get_energy(Particle **particles);
        void initialize();
        void set_alpha();

    private:
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
        std::vector< std::vector<double> > kVec;
        double *kNorm;
        int kNum;
        double alpha;
};

#endif