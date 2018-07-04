#ifndef Ewald3D_H
#define Ewald3D_H

#include "base.h"
#include "particle.h"
#include <vector>
#include <complex>

namespace energy { namespace ewald3D{
    //public:
    //Ewald3D();

    extern int kNumMax;
    extern double selfTerm;
    extern std::vector< std::vector<double> > kVec;
    extern std::complex<double> *rkVec;
    extern double *kNorm;
    extern double *resFac;
    extern double alpha;
    extern int kNum;

    template<typename T, typename G>
    double dot(T vec1, G vec2);

    template<typename T>
    static inline T erfc_x( T x );


    template<typename T>
    T erf_x( T x );

    template<typename T>
    double norm(T x);

    double get_reciprocal();
    double get_energy(Particle **particles);
    double get_particle_energy(Particle **particles, Particle* p);
    double get_self_correction(Particle *p);
    void initialize(Particle **p);
    void reset();
    void update_reciprocal(Particle *_old, Particle *_new);
    void set_alpha();
} }

/*
std::vector< std::vector<double> > Ewald3D::kVec;
std::complex<double> *Ewald3D::rkVec;
double *Ewald3D::kNorm;
double *Ewald3D::resFac;
double Ewald3D::alpha;
int Ewald3D::kNum;
int Ewald3D::kNumMax;
double Ewald3D::selfTerm;
*/
#endif