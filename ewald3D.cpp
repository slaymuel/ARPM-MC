#include "ewald3D.h"

Ewald3D::Ewald3D(){
    alpha = 8/(Base::xL);
    kNumMax = 10000;
    kNum = 0;
}
/**
 * @brief Approximation of erfc-function
 * @param x Value for which erfc should be calculated 
 * @details Reference for this approximation is found in Abramowitz and Stegun, 
 *          Handbook of mathematical functions, eq. 7.1.26
 *
 * @f[
 *     \erf(x) = 1 - (a_1t + a_2t^2 + a_3t^3 + a_4t^4 + a_5t^5)e^{-x^2} + \epsilon(x)
 * @f]
 * @f[
 *     t = \frac{1}{1 + px}
 * @f]
 * @f[
 *     |\epsilon(x)| \le 1.5\times 10^{-7}
 * @f]
 * 
 * @warning Needs testing for x < 0.
 */
template<typename T>
T Ewald3D::erfc_x( T x )
{
    static_assert(std::is_floating_point<T>::value, "type must be floating point");
    if(x < 0)
    return ( 2.0 - erfc_x(-x) );
    T t = 1.0 / (1.0 + 0.3275911 * x);
    const T a1 = 0.254829592;
    const T a2 = -0.284496736;
    const T a3 = 1.421413741;
    const T a4 = -1.453152027;
    const T a5 = 1.061405429;
    return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) * exp(-x * x);
}

/**
 * @brief Approximate 1 - erfc_x
 * @param x Value for which erf should be calculated 
 */
template<typename T>
T Ewald3D::erf_x( T x ) { 
    return (1 - erfc_x(x)); 
}

template<typename T>
double Ewald3D::norm(T x){
    double norm = 0;

    norm = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    return sqrt(norm);
}

template<typename T, typename G>
double dot(T vec1, G vec2){
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

void Ewald3D::initialize(){
    int i = 0;
    double r = 0;
    double qq = 0;
    int kx = 0;
    double kx2;
    int ky = 0;
    double ky2;
    int kz = 0;
    double kz2;
    double k2 = 0;
    resFac = (double*) malloc(kNumMax * sizeof(double));
    int kMax = 11;//8/Base::xL;
    //get k-vectors
    double factor = 1;
    std::vector<double> vec(3);

    for(kx = 0; kx <= kMax; kx++){
        for(ky = -kMax; ky <= kMax; ky++){
            for(kz = -kMax; kz <= kMax; kz++){
                factor = 1.0;
                if(kx > 0){
                    factor *= 2;
                }
                vec[0] = (2.0*PI*kx/Base::xL);
                vec[1] = (2.0*PI*ky/Base::yL);
                vec[2] = (2.0*PI*kz/Base::zL);
                k2 = dot(vec, vec);

                if(k2 != 0){
                    kVec.push_back(vec);
                    resFac[kNum] = factor * exp(-k2/(4.0 * alpha * alpha))/k2;
                    kNum++;
                }
            }
        }
    }

    printf("Found: %d k-vectors\n", kNum);
    //Calculate norms
    kNorm = (double*) malloc(kNum * sizeof(double));
    for(i = 0; i < kNum; i++){
        kNorm[i] = norm(kVec[i]);
    }
}

double Ewald3D::get_self_correction(Particle *p){
    double self = 0;
    self = p->q * p->q;
    return self;
}

double Ewald3D::get_reciprocal(Particle **p){
    double energy = 0;
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
        energy += std::norm(rho) * resFac[k];
        //printf("resfac: %lf norm: %lf kvec: %lf %lf %lf\n", resFac[k], std::norm(rho) * std::norm(rho), kVec[k][0], kVec[k][1], kVec[k][2]);
    }
    return 2 * PI * energy;    
}

double Ewald3D::get_real(Particle *p1, Particle *p2){
    int i = 0;
    double energy = 0;
    double distance = sqrt(p1->distance(p2));

    energy = erfc_x(distance * alpha)/
                distance;
    return p1->q * p2->q * energy;
}

double Ewald3D::get_energy(Particle **particles){
    double real;
    double self;
    double reciprocal;
    int j = 0;
    reciprocal = get_reciprocal(particles);

    for(int i = 0; i < Particle::numOfParticles; i++){
        j = i + 1;
        while(j < Particle::numOfParticles){
            if(j != i){
                real += get_real(particles[i], particles[j]);
                //reciprocal += get_reciprocal(particles[i], particles[j]);
            }
            j++;
        }
        //reciprocal += get_reciprocal2(particles[i]);
        self += get_self_correction(particles[i]);
    }
    real = real;
    reciprocal = 1.0/(Base::xL * Base::yL * Base::zL) * reciprocal;
    self = alpha/sqrt(PI) * self;

    printf("Real: %lf, self: %lf, reciprocal: %lf\n", real, self, reciprocal);
    return (real + reciprocal - self);
}