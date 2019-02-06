#include "levin.h"

int energy::levin::kNumMax;
int energy::levin::kNum;
int energy::levin::kMax;
//(eIn - eOut)/(eIn + eOut)
double energy::levin::gamma;//0.95;
Eigen::MatrixXd energy::levin::kVectors;
Eigen::ArrayXd energy::levin::kNorms;
Eigen::ArrayXd energy::levin::f1;
Eigen::ArrayXd energy::levin::f2;
Eigen::ArrayXd energy::levin::f3;
Eigen::ArrayXd energy::levin::f4;
Eigen::ArrayXd energy::levin::eFactors;
double energy::levin::uGamma;

void energy::levin::initialize(Particles &particles){
    kNumMax = 1000000;
    kNum = 0;
    kMax = 12;
    //(eIn - eOut)/(eIn + eOut)
    gamma = -1;//0.95;
    kVectors.resize((2 * kMax + 1) * (2 * kMax + 1) - 1, 2);
    kNorms.resize((2 * kMax + 1) * (2 * kMax + 1) - 1);
    uGamma = 0;
    gamma = -1;
    Eigen::Vector2d tempVec;

    for(int x = -kMax; x <= kMax; x++){
        for(int y = -kMax; y <= kMax; y++){
            if(x != 0 && y != 0){
                tempVec << x, y;
                kVectors.row(kNum) = tempVec;
                kNorms(kNum) = 2.0 * PI * sqrt((double)(x * x)/(Base::xL * Base::xL) + (double)(y * y)/(Base::yL * Base::yL));
                kNum++;
            }
        }
    }
    printf("Found %d k-vectors.\n", kNum);

    //Calculate f-functions
    f1 = Eigen::ArrayXd::Zero(kNum);
    f2 = Eigen::ArrayXd::Zero(kNum);
    f3 = Eigen::ArrayXd::Zero(kNum);
    f4 = Eigen::ArrayXd::Zero(kNum);
    double factor = 0;
    for(int i = 0; i < kNum; i++){
        for(int j = 0; j < particles.numOfParticles; j++){
            factor = 2.0 * PI/Base::xL * (kVectors(i, 0) * particles[j].pos[0] + kVectors(i, 1) * particles[j].pos[1]);
            f1(i) += particles[j].q * std::cos(factor) * std::exp(-kNorms(i) * (particles[j].pos[2] - Base::wall));
            f2(i) += particles[j].q * std::sin(factor) * std::exp(-kNorms(i) * (particles[j].pos[2] - Base::wall));
            f3(i) += particles[j].q * std::cos(factor) * std::exp(kNorms(i) * (particles[j].pos[2] - Base::wall));
            f4(i) += particles[j].q * std::sin(factor) * std::exp(kNorms(i) * (particles[j].pos[2] - Base::wall));
        }  
    }
    printf("Calculated f-functions\n");
    eFactors = Eigen::ArrayXd::Zero(kNum);
    for(int i = 0; i < kNum; i++){
        eFactors(i) = exp(-2.0 * kNorms(i) * (Base::zL - 2 * Base::wall));
    }

    double dipol = 0;
    for(int i = 0; i < particles.numOfParticles; i++){
        dipol += particles[i].q * (particles[i].pos[2] - Base::wall);
    }
    uGamma = -2.0 * PI/(Base::xL * Base::xL) * (dipol * dipol / (Base::zL - 2.0 * Base::wall));
    //for(int i = 0; i < kNum; i++){
        //printf("Vec %d: %lf %lf\n", i, kVectors(i,0), kVectors(i,1));
    //}
    //std::cout << kVectors << std::endl;
    //printf("Element %lf\n", kVectors(0,1));

}


void energy::levin::update_f(Particle &_old, Particle &_new){
    for(int i = 0; i < kNum; i++){
        double oldFactor = 2.0 * PI/Base::xL * (kVectors(i, 0) * _old.pos[0] + kVectors(i, 1) * _old.pos[1]);
        double newFactor = 2.0 * PI/Base::xL * (kVectors(i, 0) * _new.pos[0] + kVectors(i, 1) * _new.pos[1]);
        f1(i) -= _old.q * std::cos(oldFactor) * std::exp(-kNorms(i) * (_old.pos[2] - Base::wall));
        f1(i) += _new.q * std::cos(newFactor) * std::exp(-kNorms(i) * (_new.pos[2] - Base::wall));

        f2(i) -= _old.q * std::sin(oldFactor) * std::exp(-kNorms(i) * (_old.pos[2] - Base::wall));
        f2(i) += _new.q * std::sin(newFactor) * std::exp(-kNorms(i) * (_new.pos[2] - Base::wall));

        f3(i) -= _old.q * std::cos(oldFactor) * std::exp(kNorms(i) * (_old.pos[2] - Base::wall));
        f3(i) += _new.q * std::cos(newFactor) * std::exp(kNorms(i) * (_new.pos[2] - Base::wall));

        f4(i) -= _old.q * std::sin(oldFactor) * std::exp(kNorms(i) * (_old.pos[2] - Base::wall));
        f4(i) += _new.q * std::sin(newFactor) * std::exp(kNorms(i) * (_new.pos[2] - Base::wall));
    }
}


double energy::levin::u_gamma(Particles &particles){
    double dipol = 0;
    double chargeProd = 0;

    for(int i = 0; i < particles.numOfParticles; i++){
        dipol += particles[i].q * (particles[i].pos[2] - Base::wall);
        chargeProd += particles[i].q;
    }
    //printf("dipol: %lf, chargeProd: %lf, u_gamma: %lf\n", dipol, chargeProd, -2 * PI/(Base::xL * Base::xL) * (dipol * dipol / Base::zL - chargeProd * dipol));
    return -2.0 * PI/(Base::xL * Base::xL) * (dipol * dipol / (Base::zL - 2.0 * Base::wall) - chargeProd * dipol);
}


double energy::levin::get_polarization(){
    double ePol = 0;
    
    //ePol += uGamma;
    
    for(int i = 0; i < kNum; i++){
        //printf("eFactor: %lf    f1: %lf     f2: %lf     f3: %lf     f4: %lf     kNorm: %lf\n", eFactors(i), f1(i), f2(i), f3(i), f4(i), kNorms(i));
        
        ePol += -1/(kNorms(i) * (1 - eFactors(i))) * 
                (f1(i) * f1(i) + f2(i) * f2(i) + eFactors(i) * (f3(i) * f3(i) + f4(i) * f4(i)) - 
                2 * eFactors(i) * (f3(i) * f1(i) + f2(i) * f4(i)));
                
        //printf("%.15lf     %.15lf\n", gamma/(kNorms(i) * (1 - gamma * gamma * eFactors(i))) * (f1(i) * f1(i) + f2(i) * f2(i) + eFactors(i) * (f3(i) * f3(i) + f4(i) * f4(i))), 
        //        gamma/(kNorms(i) * (1 - gamma * gamma * eFactors(i))) * 2 * gamma * eFactors(i) * (f3(i) * f1(i) + f2(i) * f4(i)));
        //printf("ePol %lf\n", ePol);
    }
    //printf("ePol %lf\n", ePol);
    return ePol;
}


double energy::levin::get_energy(Particles &particles){
    //double eDir = energy::ewald3D::get_energy(particles);
    //printf("u_gamma: %.15lf\n", u_gamma(particles) * Base::lB);
    //printf("normalized pol: %lf\n", PI/(Base::xL * Base::xL) * get_polarization() * Base::lB);

    double ePol = (u_gamma(particles) + PI/(Base::xL * Base::xL) * get_polarization()) * Base::lB;

    //printf("direct: %lf    polarization: %lf\n", eDir, ePol);
    //printf("Levin direct: %lf    polarization: %lf\n", eDir, ePol);
    return 0;//eDir + ePol;
}


double energy::levin::get_particle_energy(Particles &particles, Particle &p){
    //double eDir = energy::ewald3D::get_particle_energy(particles, p);
    //printf("u_gamma: %lf    polarization: %lf\n", u_gamma(particles), get_polarization());

    double ePol = (u_gamma(particles) + PI/(Base::xL * Base::xL) * get_polarization()) * Base::lB;
    
    //printf("direct: %lf    polarization: %lf\n", eDir, ePol);
    return 0;//eDir + ePol;
}
