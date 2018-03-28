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
        void initialize(Particle **p);
        void update_reciprocal(Particle *_old, Particle *_new);
    private:
        template<typename T>
        T erfc_x( T x );

        template<typename T>
        T erf_x( T x );

        template<typename T>
        double norm(T x);
        
        int kNumMax;
        double get_self_correction(Particle *p);
        double get_reciprocal();
        double get_real(Particle *p1, Particle *p2);
        std::vector< std::vector<double> > kVec;
        std::complex<double> *rkVec;
        double *kNorm;
        double *resFac;
        double alpha;
        int kNum;
};

#endif