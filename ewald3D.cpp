#include "ewald3D.h"
//#include "boost/math/special_functions/erf.hpp"
//#include "omp.h"
int energy::ewald3D::kNumMax;
double energy::ewald3D::selfTerm;
std::vector< std::vector<double> > energy::ewald3D::kVec;
std::complex<double> *energy::ewald3D::rkVec = NULL;
double *energy::ewald3D::kNorm = NULL;
double *energy::ewald3D::resFac = NULL;
double energy::ewald3D::alpha;
int energy::ewald3D::kNum;




template<typename T, typename G>
double energy::ewald3D::dot(T vec1, G vec2){
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}





template<typename T>
static inline T energy::ewald3D::erfc_x( T x )
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
T energy::ewald3D::erf_x( T x ) {
    return (1 - erfc_x(x));
}





template<typename T>
double energy::ewald3D::norm(T x){
    double norm = 0;

    norm = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    return sqrt(norm);
}





void energy::ewald3D::set_alpha(){
    alpha = 5.0 / Base::xL; //8.0 / alpha
}





void energy::ewald3D::reset(){
    kVec.clear();
    free(resFac);
    free(kNorm);
    free(rkVec);
}







void energy::ewald3D::initialize(Particles &particles){
    int i = 0;
    double r = 0;
    double qq = 0;
    double k2 = 0;
    kNumMax = 1000000;
    kNum = 0;
    resFac = (double*) malloc(kNumMax * sizeof(double));
    int kMax = 6;//8/Base::xL;
    int zMax = (int) (Base::zL / Base::xL * kMax);
    printf("Wavevectors in x: %i\n", kMax);
    printf("Wavevectors in z: %i\n", zMax);
    double recCut = 4;
    //get k-vectors
    double factor = 1;
    std::vector<double> vec(3);
    //printf("Calculating k-vectors");
    for(int kx = 0; kx <= kMax; kx++){
        for(int ky = -kMax; ky <= kMax; ky++){
            for(int kz = -zMax; kz <= zMax; kz++){
                factor = 1.0;
                if(kx > 0){
                    factor *= 2;
                }

                vec[0] = (2.0 * PI * kx / Base::xL);
                vec[1] = (2.0 * PI * ky / Base::yL);
                vec[2] = (2.0 * PI * kz / Base::zL);
                //printf("%i, %i, %i\n", kx, ky, kz);
                //if(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2] < recCut * recCut){ //Spherical reciprocal space cutoff
                    k2 = dot(vec, vec);

                    if(fabs(k2) > 1e-8){// && fabs(k2) < kMax) {
                        kVec.push_back(vec);
                        resFac[kNum] = factor * exp(-k2 / (4.0 * alpha * alpha)) / k2;
                        kNum++;
                    }
                //}
            }
        }
    }
    //printf("\n");
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
        for(int i = 0; i < particles.numOfParticles; i++){
            rk.imag(std::sin(dot(particles[i].pos, kVec[k])));
            rk.real(std::cos(dot(particles[i].pos, kVec[k])));
            charge = particles[i].q;
            rk = rk * charge;
            rho += rk;
        }
        rkVec[k] = rho;
    }
    selfTerm = 0;
    for(int i = 0; i < particles.numOfParticles; i++){
        selfTerm += get_self_correction(particles[i]);
    }
    selfTerm = alpha/sqrt(PI) * selfTerm * Base::lB;
}





void energy::ewald3D::update_reciprocal(Particle &_old, Particle &_new){
    double energy = 0;
    std::complex<double> rk_new;
    std::complex<double> rk_old;
    std::complex<double> charge = _new.q;

    for(int k = 0; k < kNum; k++){

        rk_new.imag(std::sin(dot(_new.pos, kVec[k])));
        rk_new.real(std::cos(dot(_new.pos, kVec[k])));

        rk_old.imag(std::sin(dot(_old.pos, kVec[k])));
        rk_old.real(std::cos(dot(_old.pos, kVec[k])));

        rkVec[k] -= rk_old * charge;
        rkVec[k] += rk_new * charge;
    }
    //printf("rk: %f\n", rkVec[kNum - 1].real() + rkVec[kNum - 1].imag());
}






double energy::ewald3D::get_reciprocal(){
    double energy = 0;

    for(int k = 0; k < kNum; k++){
        energy += std::norm(rkVec[k]) * resFac[k];
    }
    return energy;
}






double energy::ewald3D::get_self_correction(Particle &p){
    double self = p.q * p.q;
    return self;
}




double energy::ewald3D::get_energy(Particles &particles){
        double real = 0;
        //double self = 0;
        double reciprocal = 0;
        //double dipoleMoment[3] = {0, 0, 0};
        Eigen::Vector3d dipoleMoment;
        dipoleMoment.setZero();
        double corr = 0;
        double distance = 0;
        double energy = 0;

        reciprocal = get_reciprocal();

        for(int i = 0; i < particles.numOfParticles; i++){
                for(int j = i + 1; j < particles.numOfParticles; j++){
                    distance = particles.distances[i][j];

                    //if(distance <= 10){
                        energy = erfc_x(distance * alpha) / distance;
                        real += particles[i].q * particles[j].q * energy;
                    //}
                }
                
                dipoleMoment += particles[i].q * particles[i].pos;

                //reciprocal += get_reciprocal2(particles[i]);
                //self += get_self_correction(particles[i]);
        }

        corr = dipoleMoment[2];
        corr *= corr * 2.0 * PI / (Base::volume);
        reciprocal = 2.0 * PI / (Base::xL * Base::yL * Base::zL) * reciprocal;
        //self = alpha/sqrt(PI) * self;
        //printf("Dipole moment: %lf\n", corr);
        //printf("self term: %lf\n", selfTerm);
        printf("Real: %lf, self: %lf, reciprocal: %lf, correction: %lf, tot: %lf\n", real * Base::lB, selfTerm, reciprocal * Base::lB, corr * Base::lB, Base::lB * (real + reciprocal + corr) - selfTerm);
        printf("Tinfoil Energy: %.10lf\n", (real + reciprocal) - selfTerm/Base::lB);
        //printf("Vacuum Energy: %.10lf\n", (real + reciprocal + corr) - selfTerm/Base::lB);
        return Base::lB * (real + reciprocal + corr) - selfTerm;    //vacuum
        //return Base::lB * (real + reciprocal) - selfTerm;   //tinfoil
}






double energy::ewald3D::get_particle_energy(Particles &particles, Particle &p){
    double real = 0;
    double self = 0;
    double reciprocal = 0;
    //double dipoleMoment[3] = {0, 0, 0};
    Eigen::Vector3d dipoleMoment;
    dipoleMoment.setZero();
    double corr = 0;
    double distance = 0;
    double energy = 0;
    //printf("alpha is: %lf\n", alpha);

    reciprocal = get_reciprocal();

    for(int i = p.index + 1; i < particles.numOfParticles; i++){
        distance = particles.distances[p.index][i];

        //if(distance <= 10){
            energy = erfc_x(distance * alpha) / distance;
            real += particles[i].q * p.q * energy;
        //}
        dipoleMoment += particles[i].q * particles[i].pos;
    }
    for(int i = 0; i < p.index; i++){
        distance = particles.distances[i][p.index];
        
        //if(distance <= 10){
            energy = erfc_x(distance * alpha) / distance;
            real += particles[i].q * p.q * energy;
        //}
        dipoleMoment += particles[i].q * particles[i].pos;
    }

    dipoleMoment += p.q * p.pos;
    corr = dipoleMoment[2];
    corr *= corr;
    corr = 2.0 * PI * corr/(Base::volume);
    reciprocal = 2.0 * PI/(Base::xL * Base::yL * Base::zL) * reciprocal;

    return Base::lB * (real + reciprocal + corr) - selfTerm;    //vacuum
    //return Base::lB * (real + reciprocal) - selfTerm;   //tinfoil
}