#include "ewald.h"

Ewald::Ewald(){
    double *kNorm;
    static int kNum = 0;
    double alpha = 5/Base::xL;
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
T Ewald::erfc_x( T x )
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
T Ewald::erf_x( T x ) { 
    return (1 - erfc_x(x)); 
}

template<typename T>
T Ewald::ewald_F(T x){

}

template<typename T>
double Ewald::norm(T x){
    double norm = 0;

    norm = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    return sqrt(norm);
}

void Ewald::initialize(){
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
    int kNumMax = 400;
    double kMax = 8*PI/Base::xL;
    
    //get k-vectors
    for(kx = 0; kx < kMax; kx++){
        for(int ky = -kMax + 1; ky < kNumMax; ky++){
            kx2 = kx*kx;
            ky2 = ky*ky;
            k2 = kx2 + ky2;
            if(k2 != 0 && k2 < kMax){
                std::vector<double> vec;
                vec.push_back(2*PI*kx/Base::xL);
                vec.push_back(2*PI*ky/Base::yL);
                vec.push_back(0);
                kVec.push_back(vec);
                //kVec[kNum] = exp(-k2/());
                kNum++;
            }
        }
    }

    //Calculate norms
    kNorm = (double*) malloc(kNum * sizeof(double));
    for(i = 0; i < kNum; i++){
        kNorm[i] = norm(kVec[i]);
    }
}

double Ewald::get_self_correction(Particle *p){
    double self = 0;
    self = alpha/sqrt(PI) * p->q * p->q;
    return self;
}

double Ewald::f(double norm, double zDist){
    double f = 0;
    f = exp(norm * zDist) * erfc_x(alpha * zDist + norm/(2 * alpha)) +
        exp(-norm * zDist) * erfc(-alpha * zDist + norm/(2 * alpha));

    return f;
}

double Ewald::g(Particle *p1, Particle *p2){
    double zDist = p1->distance_z(p2);
    double g = 0;
    g = zDist * erfc_x(alpha * zDist) * exp(-(zDist*alpha*zDist*alpha))/(alpha*sqrt(PI));
}

double Ewald::get_reciprocal(Particle *p1, Particle *p2){
    int i = 0;
    double energy = 0;
    for(i = 0; i < kNum; i++){
        energy += cos(kNorm[i] * p1->distance_xy(p2)) * f(kNorm[i], p1->distance_z(p2));
    }
    return energy * PI/(Base::xL * Base::xL);
}

double Ewald::get_real(Particle *p1, Particle *p2){
    int i = 0;
    double energy = 0;
    energy = erfc_x(alpha * p1->distance_xy(p2))/
                p1->distance_xy(p2);
    return energy;
}

double Ewald::get_energy(Particle *p, Particle **particles){
    int i = 0;
    double energy = 0;

    for(i = 0; i < Particle::numOfParticles; i++){
        energy += get_real(p, particles[i]) + get_reciprocal(p, particles[i]) - g(p, particles[i]) - get_self_correction(p);
    }
    return 1/2 * energy;
}