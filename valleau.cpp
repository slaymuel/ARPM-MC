#include "valleau.h"


void energy::valleau::update_charge_vector(Particle **particles){
    double binWidth = 0.1; //Angstrom
    int numOfBins = Base::zL/binWidth; //number of bins
    Eigen::VectorXd pDensity(numOfBins);
    Eigen::VectorXd nDensity(numOfBins);
    for(int i = 0; i < numOfBins; i++){
        nDensity[i] = 0;
        pDensity[i] = 0;
        
    }

    for(int i = 0; i < Particle::numOfParticles; i++){
        if(particles[i]->q < 0){
            pDensity[(int)(particles[i]->pos[2]/binWidth)]++;
        }
        else{
            nDensity[(int)(particles[i]->pos[2]/binWidth)]++;
        }
    }

    chargeVector = pDensity - nDensity; //charge vector

    //Symmetrize charge vector
    int j = chargeVector.size() - 1;
    double avg;
    for(int i = 0; i < chargeVector.size()/2; i++){
        avg = (chargeVector[i] + chargeVector[j])/2;
        chargeVector[i] = avg;
        chargeVector[j] = avg;
        j--;
    }
    //std::cout << chargeVector << std::endl;

    double dZ = 0.1;
    for(int i = 0; i < chargeVector.size() / dZ; i++){
        //-2 * PI * i * dZ * chargeVector[i/dZ];
    }
}

double energy::valleau::phiw(double z){
    return z;
}

void energy::valleau::update_potential(){
    double dz = 0.5;
    double diffz = 0;
    ext.resize(chargeVector.size());
    for(int i = 0; i < Base::zL; i++){
        for(int j = 0; j < Base::zL; j++){
            diffz = fabs(j * dz - i * dz);
            ext[i] += chargeVector[j] * (-2 * PI * diffz - phiw(diffz));
        }
        ext[i] *= Base::lB * dz;
    }
}