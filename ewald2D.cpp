#include "ewald2D.h"
#include "omp.h"
std::vector< Eigen::Vector3d > energy::ewald2D::kVec;
double *energy::ewald2D::kNorm;
int energy::ewald2D::kNum;
double energy::ewald2D::alpha;

void energy::ewald2D::set_alpha(){
    alpha = 5.0 / Base::xL;
}

template<typename T>
T energy::ewald2D::erfc_x( T x ){
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
T energy::ewald2D::erf_x( T x ) { 
    return (1 - erfc_x(x)); 
}

template<typename T, typename G>
double dot(T vec1, G vec2){
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

template<typename T>
double energy::ewald2D::norm(T x){
    double norm = 0;

    norm = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    return sqrt(norm);
}

void energy::ewald2D::initialize(){
    kNum = 0;
    int i = 0;
    int kx = 0;
    int ky = 0;
    int kz = 0;
    double k2 = 0;
    int kMax = 5;//8/Base::xL;

    //get k-vectors
    Eigen::Vector3d vec;

    for(kx = -kMax; kx <= kMax; kx++){
        for(ky = -kMax; ky <= kMax; ky++){
            vec[0] = (2.0 * PI * kx/Base::xL);
            vec[1] = (2.0 * PI * ky/Base::yL);
            vec[2] = 0;
            k2 = vec.dot(vec);
            if(fabs(k2) > 1e-5){// && fabs(k2) < kMax){
                kVec.push_back(vec);
                kNum++;
            }
        }
    }

    printf("2D: Found: %d k-vectors\n", kNum);
    //Calculate norms
    kNorm = (double*) malloc(kNum * sizeof(double));
    for(i = 0; i < kNum; i++){
        kNorm[i] = norm(kVec[i]);
    }
}

double energy::ewald2D::dipole_correction(Particle &p){
    return p.q * p.pos[2];
}

double energy::ewald2D::get_self_correction(Particle &p){
    double self = 0;
    self = p.q * p.q;
    return self;
}

double energy::ewald2D::f(double norm, double zDist){
    double f = 0;
    f = exp(norm * zDist) * erfc_x(alpha * zDist + norm/(2 * alpha)) +
        exp(-norm * zDist) * erfc_x(-alpha * zDist + norm/(2 * alpha));

    return f / (2.0 * norm);
}

double energy::ewald2D::g(Particle &p1, Particle &p2){
    //double zDist = sqrt(p1.distance_z(p2));
    double zDist = p1.pos[2] - p2.pos[2];
    double g = 0;
    g = zDist * erf_x(alpha * zDist) + exp(-(zDist * zDist * alpha * alpha)) / (alpha * sqrt(PI));

    return g;
}

double energy::ewald2D::get_reciprocal(Particle &p1, Particle &p2){
    double energy = 0;
    //double distance = p1.distance(p2);
    //std::vector<double> dispVec(3);
    std::complex<double> rk;
    double zDist = p1.pos[2] - p2.pos[2];
    Eigen::Vector3d dispVec = p1.pos - p2.pos;
    //dispVec[0] = p1.pos[0] - p2.pos[0];
    //dispVec[1] = p1.pos[1] - p2.pos[1];
    //dispVec[2] = p1.pos[2] - p2.pos[2];
    
 
    for(int i = 0; i < kNum; i++){
        rk.imag(std::sin(dispVec.dot(kVec[i])));
        rk.real(std::cos(dispVec.dot(kVec[i])));
    
        energy += rk.real() * (exp(kNorm[i] * zDist) * erfc_x(alpha * zDist + kNorm[i]/(2 * alpha)) +
            exp(-kNorm[i] * zDist) * erfc_x(-alpha * zDist + kNorm[i]/(2 * alpha))) / (2.0 * kNorm[i]);

        //energy += rk.real() * f(kNorm[i], p1.pos[2] - p2.pos[2]);
    }
    return energy;
}

double energy::ewald2D::get_real(Particle &p1, Particle &p2){
    double energy = 0;
    double distance = sqrt(p1.distance_xy(p2));
    energy = erfc_x(alpha * distance)/distance;
    if(distance >= Base::xL/2 && fabs(energy) >= 1e-4){
        printf("Energy is %lf at the boundary maybe you should increase alpha?\n", energy);
        //exit(1);
        //printf("Distance: %lf\n", distance);
    }
    return energy;
}

double energy::ewald2D::get_particle_energy(std::vector<Particle> &particles, Particle &p){
    double energy = 0;
    double distance = 0;
    double real = 0;
    double self = 0;
    double reciprocal = 0;
    double reci = 0;
    double dipCorr = 0;
    double done = 0;
    double gE = 0;
    int k = 0;

    double stime = omp_get_wtime();
    for(int i = p.index + 1; i < Particle::numOfParticles; i++){
        distance = Particle::distances[p.index][i];
        real += particles[i].q * p.q * erfc_x(alpha * distance) / distance;

    }
    
    for(int i = 0; i < p.index; i++){
        distance = Particle::distances[i][p.index];
        real += particles[i].q * p.q * erfc_x(alpha * distance) / distance;
    }

    for(int i = 0; i < Particle::numOfParticles; i++){

        for(int j = 0; j < Particle::numOfParticles; j++){
            double zDist = particles[i].pos[2] - particles[j].pos[2];
            reciprocal += particles[i].q * particles[j].q * get_reciprocal(particles[i], particles[j]) * 1.0 / 2.0;

            //reci += particles[i].q * particles[j].q * get_reciprocal(particles[i], particles[j]);
            gE += particles[i].q * particles[j].q * (zDist * erf_x(alpha * zDist) + exp(-(zDist * zDist * alpha * alpha)) / (alpha * sqrt(PI)));
            //gE += particles[i].q * particles[j].q * g(particles[i], particles[j]);
        }
        //dipCorr += dipole_correction(particles[i]);
        self += get_self_correction(particles[i]);
        //printf("Rec: %lf g: %lf: reci: %lf Done: %lf\r", reciprocal, gE, reci,(double)i/Particle::numOfParticles * 100);
        //fflush(stdout);
    }
    reciprocal = (reciprocal - gE) * PI/(Base::xL * Base::yL);
    
    //dipCorr = -2 * PI/(Base::xL * Base::yL * Base::zL) * dipCorr * dipCorr;
    self = alpha / sqrt(PI) * self;
    energy = (real + reciprocal - self);// + dipCorr);
    //printf("Real: %lf, self: %lf, reciprocal: %lf dipCorr: %lf\n", real * Base::lB, self, reciprocal, dipCorr);
    return energy * Base::lB;
}

double energy::ewald2D::get_energy(std::vector<Particle> &particles){
    double energy = 0;
    double distance = 0;
    double real = 0;
    double self = 0;
    double reciprocal = 0;
    double reci = 0;
    double dipCorr = 0;
    double done = 0;
    double gE = 0;
    int k = 0;

    for(int i = 0; i < Particle::numOfParticles; i++){

        k = i + 1;
        while(k < Particle::numOfParticles){
            distance = Particle::distances[i][k];
            if(k != i){
                //real += particles[i].q * particles[k].q * get_real(particles[i], particles[k]);
                real += particles[i].q * particles[k].q * erfc_x(alpha * distance)/distance;
            }
            k++;
        }

        for(int j = 0; j < Particle::numOfParticles; j++){
            double zDist = particles[i].pos[2] - particles[j].pos[2];
            reciprocal += particles[i].q * particles[j].q * get_reciprocal(particles[i], particles[j]) * 1.0 / 2.0;
            //reci += particles[i].q * particles[j].q * get_reciprocal(particles[i], particles[j]);
            //gE += particles[i].q * particles[j].q * g(particles[i], particles[j]);
            gE += particles[i].q * particles[j].q * (zDist * erf_x(alpha * zDist) + exp(-(zDist * zDist * alpha * alpha)) / (alpha * sqrt(PI)));
        }
        //dipCorr += dipole_correction(particles[i]);
        self += get_self_correction(particles[i]);
        //printf("Rec: %lf g: %lf: reci: %lf Done: %lf\r", reciprocal, gE, reci,(double)i/Particle::numOfParticles * 100);
        //fflush(stdout);
    }
    reciprocal = (reciprocal - gE) * PI/(Base::xL * Base::yL);
    dipCorr = -2 * PI/(Base::xL * Base::yL * Base::zL) * dipCorr * dipCorr;
    self = alpha/sqrt(PI) * self;
    energy = (real + reciprocal - self);// + dipCorr);
    //printf("Real: %lf, self: %lf, reciprocal: %lf dipCorr: %lf\n", real * Base::lB, self, reciprocal, dipCorr);
    return energy * Base::lB;
}
