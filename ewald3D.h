#ifndef Ewald3D_H
#define Ewald3D_H

#include "base.cpp"
#include "particle.h"
#include <vector>
#include <complex>
#include <eigen3/Eigen/Dense>

class Ewald3D{
    public:
        Ewald3D();
        double get_energy(Particle **particles);
        void initialize();

    private:
        template<typename T>
        T erfc_x( T x );

        template<typename T>
        T erf_x( T x );

        template<typename T>
        T Ewald3D_F(T x);

        template<typename T>
        double norm(T x);
        int kNumMax;
        double get_self_correction(Particle *p);
        double f(double norm, double zDist);
        double g(Particle *p1, Particle *p2);
        double get_reciprocal(Particle *p1, Particle *p2);
        double get_reciprocal2(Particle *p);
        double p(Particle **p);
        double get_real(Particle *p1, Particle *p2);
        std::vector< std::vector<double> > kVec;
        double *kNorm;
        double *resFac;
        double alpha;
        int kNum;
};

#endif