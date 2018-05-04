#include "levin.h"

Levin::Levin(){
    kNumMax = 1000000;
    kNum = 0;
    kMax = 100;
    //(eIn - eOut)/(eIn + eOut)
    gamma = 0.95;
    kVectors = Eigen::MatrixXd::Zero((2 * kMax + 1) * (2 * kMax + 1), 2);
    kNorms = Eigen::ArrayXd::Zero((2 * kMax + 1) * (2 * kMax + 1));
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
    //for(int i = 0; i < kNum; i++){
        //printf("Vec %d: %lf %lf\n", i, kVectors(i,0), kVectors(i,1));
    //}
    //std::cout << kVectors << std::endl;
    //printf("Element %lf\n", kVectors(0,1));

}

double Levin::get_polarization(){
    double ePol = 0;
    Eigen::Vector2d v(1.0, 2.0);

    return ePol;
}

double Levin::get_energy(){
    double eDir = 0;
    double ePol = 0;

    return eDir + ePol;
}