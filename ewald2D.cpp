#include "ewald2D.h"

Ewald2D::Ewald2D(){
    kNum = 0;
    alpha = 8/Base::xL;
}

template<typename T>
T Ewald2D::erfc_x( T x ){
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

template<typename T>
T Ewald2D::erf_x( T x ) { 
    return (1 - erfc_x(x)); 
}

template<typename T, typename G>
double dot(T vec1, G vec2){
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

template<typename T>
double Ewald2D::norm(T x){
    double norm = 0;

    norm = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    return sqrt(norm);
}

void Ewald2D::initialize(){
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
    int kMax = 11;//8/Base::xL;
    //get k-vectors
    std::vector<double> vec(3);

    for(kx = 0; kx <= kMax; kx++){
        for(ky = -kMax; ky <= kMax; ky++){
            vec[0] = (2.0 * PI * kx/Base::xL);
            vec[1] = (2.0 * PI * ky/Base::yL);
            vec[2] = 0;
            k2 = dot(vec, vec);
            if(k2 > 0){
                kVec.push_back(vec);
                kNum++;
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

double Ewald2D::dipole_correction(Particle *p){
    return p->q * p->pos[2];
}

double Ewald2D::get_self_correction(Particle *p){
    double self = 0;
    self = p->q * p->q;
    return self;
}

double Ewald2D::f(double norm, double zDist){
    double f = 0;
    f = exp(norm * zDist) * erfc_x(alpha * zDist + norm/(2 * alpha)) +
        exp(-norm * zDist) * erfc_x(-alpha * zDist + norm/(2 * alpha));

    return f/(2 * norm);
}

double Ewald2D::g(Particle *p1, Particle *p2){
    double zDist = sqrt(p1->distance_z(p2));
    double g = 0;
    g = zDist * erf_x(alpha * zDist) + exp(-(zDist*zDist*alpha*alpha))/(alpha * sqrt(PI));

    return g;
}

double Ewald2D::get_reciprocal(Particle *p1, Particle *p2){
    double energy = 0;
    //double distance = p1->distance(p2);
    std::vector<double> dispVec(3);
    std::complex<double> rk;
    dispVec[0] = p1->pos[0] - p2->pos[0];
    dispVec[1] = p1->pos[1] - p2->pos[1];
    dispVec[2] = p1->pos[2] - p2->pos[2];
    
    // printf("p1: %lf %lf %lf\n", p1->pos[0], p1->pos[1], p1->pos[2]);
    // printf("p2: %lf %lf %lf\n", p2->pos[0], p2->pos[1], p2->pos[2]);
    // printf("dispVec: %lf %lf %lf\n\n", dispVec[0], dispVec[1], dispVec[2]);
    for(int i = 0; i < kNum; i++){
        rk.imag(std::sin(dot(dispVec, kVec[i])));
        rk.real(std::cos(dot(dispVec, kVec[i])));
        // printf("kVec: %lf %lf %lf\n", kVec[i][0], kVec[i][1], kVec[i][2]);
        // printf("dot: %lf\n", dot(dispVec, kVec[i]));
        // printf("real: %lf\n", rk.real());
        rk = rk/kNorm[i];
        energy += rk.real() * f(kNorm[i], sqrt(p1->distance_z(p2)));
        //printf("real after: %lf\n\n", rk.real());
        //energy += -1 * g(p1, p2);
    }
    return energy * PI/(Base::xL * Base::xL);
}

double Ewald2D::get_real(Particle *p1, Particle *p2){
    double energy = 0;
    double distance = sqrt(p1->distance(p2));
    energy = erfc_x(alpha * distance)/distance;
    return energy;
}

double Ewald2D::get_energy(Particle **particles){
    double energy = 0;
    double real = 0;
    double self = 0;
    double reciprocal = 0;
    double dipCorr = 0;

    for(int i = 0; i < Particle::numOfParticles; i++){
        for(int j = 0; j < Particle::numOfParticles; j++){
            if(i != j){
                real += particles[i]->q * particles[j]->q * get_real(particles[i], particles[j]);
            }
            reciprocal += particles[i]->q * particles[j]->q * get_reciprocal(particles[i], particles[j]);
            reciprocal += -1 * particles[i]->q * particles[j]->q * g(particles[i], particles[j]);
            printf("%lf\n", -1 * particles[i]->q * particles[j]->q * g(particles[i], particles[j]));
        }
        dipCorr += dipole_correction(particles[i]);
        self += get_self_correction(particles[i]);
    }

    dipCorr = -2 * PI/(Base::xL * Base::yL * Base::zL) * dipCorr * dipCorr;
    self = alpha/sqrt(PI) * self;
    real = 1.0/2.0 * real;
    reciprocal = 1.0/2.0 * reciprocal;
    energy = real + reciprocal - self + dipCorr;
    printf("Real: %lf, self: %lf, reciprocal: %lf dipCorr: %lf\n", real, self, reciprocal, dipCorr);
    return energy;
}