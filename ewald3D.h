#ifndef Ewald3D_H
#define Ewald3D_H

#include "base.cpp"
#include "particle.h"
#include <vector>
#include <complex>

class Ewald3D{
    //public:
    //Ewald3D();
public:
    Ewald3D(){
        kNumMax = 1000000;
        kNum = 0;
    }
    int kNumMax;
    double selfTerm;
    std::vector< std::vector<double> > kVec;
    std::complex<double> *rkVec;
    double *kNorm;
    double *resFac;
    double alpha;
    int kNum;

    template<typename T, typename G>
    double dot(T vec1, G vec2){
        return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
    }

    template<typename T>
    static inline T erfc_x( T x )
    {
        //static_assert(std::is_floating_point<T>::value, "type must be floating point");
        if(x < 0){
            return ( 2.0 - erfc_x(-x) );
        }
        T t = 1.0 / (1.0 + 0.3275911 * x);
        const T a1 = 0.254829592;
        const T a2 = -0.284496736;
        const T a3 = 1.421413741;
        const T a4 = -1.453152027;
        const T a5 = 1.061405429;
        return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) *  exp(-x * x);
        //return 1;
    }


    template<typename T>
    T erf_x( T x ) {
        return (1 - erfc_x(x));
    }


    template<typename T>
    double norm(T x){
        double norm = 0;

        norm = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return sqrt(norm);
    }


    double get_reciprocal(){
        double energy = 0;

        for(int k = 0; k < kNum; k++){
            energy += std::norm(rkVec[k]) * resFac[k];
        }
        return energy;
    }

    inline double get_energy(Particle **particles){
            double real = 0;
            double self = 0;
            double reciprocal = 0;
            //double dipoleMoment[3] = {0, 0, 0};
            Eigen::Vector3d dipoleMoment;
            double corr = 0;
            double distance = 0;
            double energy = 0;

            reciprocal = get_reciprocal();
            //omp_set_num_threads(4);

            //double stime = omp_get_wtime();
            //#pragma omp parallel for reduction(+:real) schedule(dynamic) private(energy, distance, da, erfcRes)
            for(int i = 0; i < Particle::numOfParticles; i++){
                    for(int j = i + 1; j < Particle::numOfParticles; j++){
                        distance = Particle::distances[i][j];
                        energy = erfc_x(distance * alpha) / distance;
                        real += particles[i]->q * particles[j]->q * energy;
                    }
                    /*
                    dipoleMoment[0] += particles[i]->q * particles[i]->pos[0];
                    dipoleMoment[1] += particles[i]->q * particles[i]->pos[1];
                    dipoleMoment[2] += particles[i]->q * particles[i]->pos[2];
                     */
                    //reciprocal += get_reciprocal2(particles[i]);
                    //self += get_self_correction(particles[i]);
            }

            //corr = dipoleMoment.norm();
            //corr *= corr;
            //corr = 2 * PI * corr/(3 * Base::xL * Base::yL * Base::zL);
            reciprocal = 2 * PI/(Base::xL * Base::yL * Base::zL) * reciprocal;
            //self = alpha/sqrt(PI) * self;
            //printf("Dipole moment: %lf\n", corr);
            //printf("self term: %lf\n", selfTerm);
            //printf("Real: %lf, self: %lf, reciprocal: %lf\n", real, self, reciprocal);
            //return Base::lB * (real + reciprocal + corr) - selfTerm;
            return Base::lB * (real + reciprocal) - selfTerm;
    }

    inline double get_energy(Particle **particles, Particle* p){
        double real = 0;
        double self = 0;
        double reciprocal = 0;
        //double dipoleMoment[3] = {0, 0, 0};
        Eigen::Vector3d dipoleMoment;
        double corr = 0;
        double distance = 0;
        double energy = 0;

        //printf("alpha is: %lf\n", alpha);

        reciprocal = get_reciprocal();

        for(int i = p->index + 1; i < Particle::numOfParticles; i++){
            distance = Particle::distances[p->index][i];
            energy = erfc_x(distance * alpha) / distance;
            real += particles[i]->q * p->q * energy;
        }
        for(int i = 0; i < p->index; i++){
            distance = Particle::distances[i][p->index];
            energy = erfc_x(distance * alpha) / distance;
            real += particles[i]->q * p->q * energy;
        }
        reciprocal = 2 * PI/(Base::xL * Base::yL * Base::zL) * reciprocal;
        return Base::lB * (real + reciprocal) - selfTerm;
    }

    double get_self_correction(Particle *p){
        double self = 0;
        self = p->q * p->q;
        return self;
    }

    void initialize(Particle **p){
        int i = 0;
        double r = 0;
        double qq = 0;
        double kx2;
        double ky2;
        double kz2;
        double k2 = 0;
        resFac = (double*) malloc(kNumMax * sizeof(double));
        int kMax = 3;//8/Base::xL;
        //get k-vectors
        double factor = 1;
        std::vector<double> vec(3);
        printf("Calculating k-vectors");
        for(int kx = 0; kx <= kMax; kx++){
            for(int ky = -kMax; ky <= kMax; ky++){
                for(int kz = -kMax; kz <= kMax; kz++){
                    factor = 1.0;
                    if(kx > 0){
                        factor *= 2;
                    }
                    vec[0] = (2.0*PI*kx/Base::xL);
                    vec[1] = (2.0*PI*ky/Base::yL);
                    vec[2] = (2.0*PI*kz/Base::zL);
                    k2 = dot(vec, vec);

                    if(fabs(k2) > 1e-5){// && fabs(k2) < kMax) {
                        kVec.push_back(vec);
                        resFac[kNum] = factor * exp(-k2/(4.0 * alpha * alpha))/k2;
                        kNum++;
                    }
                }
            }
        }
        printf("\n");
        printf("3D: Found: %d k-vectors\n", kNum);
        //Calculate norms
        kNorm = (double*) malloc(kNum * sizeof(double));
        for(i = 0; i < kNum; i++){
            kNorm[i] = norm(kVec[i]);
        }

        rkVec = (std::complex<double>*) malloc(kNum * sizeof(std::complex<double>));
        std::complex<double> rho;
        std::complex<double> rk;
        std::complex<double> charge;
        for(int k = 0; k < kNum; k++){
            rho = 0;
            for(int i = 0; i < Particle::numOfParticles; i++){
                rk.imag(std::sin(dot(p[i]->pos, kVec[k])));
                rk.real(std::cos(dot(p[i]->pos, kVec[k])));
                charge = p[i]->q;
                rk = rk * charge;
                rho += rk;
            }
            rkVec[k] = rho;
        }
        selfTerm = 0;
        for(int i = 0; i < Particle::numOfParticles; i++){
            selfTerm += get_self_correction(p[i]);
        }
        selfTerm = alpha/sqrt(PI) * selfTerm * Base::lB;
    }


    void update_reciprocal(Particle *_old, Particle *_new){
        double energy = 0;
        std::complex<double> rk_new;
        std::complex<double> rk_old;
        std::complex<double> charge = _new->q;

        for(int k = 0; k < kNum; k++){

            rk_new.imag(std::sin(dot(_new->pos, kVec[k])));
            rk_new.real(std::cos(dot(_new->pos, kVec[k])));

            rk_old.imag(std::sin(dot(_old->pos, kVec[k])));
            rk_old.real(std::cos(dot(_old->pos, kVec[k])));

            rkVec[k] -= rk_old * charge;
            rkVec[k] += rk_new * charge;
        }
        //printf("rk: %f\n", rkVec[kNum - 1].real() + rkVec[kNum - 1].imag());
    }


    void set_alpha(){
        alpha = 5/Base::xL;
    }
};

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