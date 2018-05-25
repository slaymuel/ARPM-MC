#include "levin.h"

Levin::Levin(){
    kNumMax = 1000000;
    kNum = 0;
    kMax = 10;
    //(eIn - eOut)/(eIn + eOut)
    gamma = -1;//0.95;
    kVectors = Eigen::MatrixXd::Zero((2 * kMax + 1) * (2 * kMax + 1), 2);
    kNorms = Eigen::ArrayXd::Zero((2 * kMax + 1) * (2 * kMax + 1));
    uGamma = 0;
}

void Levin::initialize(Particle **particles){
    Eigen::Vector2d tempVec;
    for(int x = -kMax; x <= kMax; x++){
        for(int y = -kMax; y <= kMax; y++){
            tempVec << x, y;
            //std::cout << tempVec << std::endl;
            //printf("Nr: %d\n", x + y + 2 * kMax + (x + kMax) * kMax);
            kVectors.row(kNum) = tempVec;
            kNorms(kNum) = 2 * PI * sqrt(x * x/(Base::xL * Base::xL) + y * y/(Base::yL * Base::yL));
            kNum++;
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
        for(int j = 0; j < Particle::numOfParticles; j++){
            factor = 2 * PI/Base::zL * (kVectors(i, 0) * particles[j]->pos[0] + kVectors(i, 1) * particles[j]->pos[1]);
            f1(i) += particles[j]->q * cos(factor) * exp(-kNorms(i) * particles[j]->pos[2]);
            f2(i) += particles[j]->q * sin(factor) * exp(-kNorms(i) * particles[j]->pos[2]);
            f3(i) += particles[j]->q * cos(factor) * exp(kNorms(i) * particles[j]->pos[2]);
            f4(i) += particles[j]->q * sin(factor) * exp(kNorms(i) * particles[j]->pos[2]);
        }  
    }
    printf("Calculated f-functions\n");
    eFactors = Eigen::ArrayXd::Zero(kNum);
    for(int i = 0; i < kNum; i++){
        eFactors(i) = exp(2 * kNorms(i) * Base::zL);
    }

    double dipol = 0;
    for(int i = 0; i < Particle::numOfParticles; i++){
        dipol += particles[i]->q * particles[i]->pos[2];
    }
    uGamma = -2 * PI/(Base::zL * Base::zL) * (dipol * dipol / (Base::zL));
    //for(int i = 0; i < kNum; i++){
        //printf("Vec %d: %lf %lf\n", i, kVectors(i,0), kVectors(i,1));
    //}
    //std::cout << kVectors << std::endl;
    //printf("Element %lf\n", kVectors(0,1));

}

void Levin::update_f(Particle *_old, Particle *_new){
    for(int i = 0; i < kNum; i++){
        double oldFactor = 2 * PI/Base::zL * (kVectors(i, 0) * _old->pos[0] + kVectors(i, 1) * _old->pos[1]);
        double newFactor = 2 * PI/Base::zL * (kVectors(i, 0) * _new->pos[0] + kVectors(i, 1) * _new->pos[1]);
        f1(i) -= _old->q * cos(oldFactor) * exp(-kNorms(i) * _old->pos[2]);
        f1(i) += _new->q * cos(newFactor) * exp(-kNorms(i) * _new->pos[2]);

        f2(i) -= _old->q * sin(oldFactor) * exp(-kNorms(i) * _old->pos[2]);
        f2(i) += _new->q * sin(newFactor) * exp(-kNorms(i) * _new->pos[2]);

        f3(i) -= _old->q * cos(oldFactor) * exp(kNorms(i) * _old->pos[2]);
        f3(i) += _new->q * cos(newFactor) * exp(kNorms(i) * _new->pos[2]);

        f4(i) -= _old->q * sin(oldFactor) * exp(kNorms(i) * _old->pos[2]);
        f4(i) += _new->q * sin(newFactor) * exp(kNorms(i) * _new->pos[2]);
    }
}

double Levin::get_polarization(){
    double ePol = 0;
    Eigen::Vector2d v(1.0, 2.0);
    
    ePol += uGamma;

    for(int i = 0; i < kNum; i++){
        ePol += gamma/(kNorms(i) * (1 - gamma * gamma * eFactors(i))) * 
                f1(i) * f1(i) + f2(i) * f2(i) + eFactors(i) * (f3(i) * f3(i) + f4(i) * f4(i));
        ePol += 2 * gamma * eFactors(i) * (f3(i) * f1(i) + f2(i) * f4(i));
    }

    return ePol;
}

double Levin::get_energy(){
    double eDir = 0;
    double ePol = 0;

    return eDir + ePol;
}