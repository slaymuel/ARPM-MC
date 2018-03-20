#include "base.cpp"
#include "Particle.h"
#include <vector>

class Ewald{
    public:
        Ewald();
        double get_energy(Particle *p, Particle **particles);

    private:
        template<typename T>
        T erfc_x( T x );

        template<typename T>
        T erf_x( T x );

        template<typename T>
        T ewald_F(T x);

        template<typename T>
        double norm(T x);

        void initialize();
        double get_self_correction(Particle *p);
        double f(double norm, double zDist);
        double g(Particle *p1, Particle *p2);
        double get_reciprocal(Particle *p1, Particle *p2);
        double get_real(Particle *p1, Particle *p2);
        std::vector< std::vector<double> > kVec;
        double *kNorm;
        static int kNum;
        double alpha;
};