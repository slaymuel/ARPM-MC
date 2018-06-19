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

    double dZ = 0.1;
    for(int i = 0; i < chargeVector.size() / dZ; i++){
        //-2 * PI * i * dZ * chargeVector[i/dZ];
    }
}

double energy::valleau::phiw(double z){
    double a = Base::xL/2;
    double asq = a * a;
    double b = Base::zL/2;
    double s = 1;
    double zsbsq = (z + s * b) * (z + s * b);
    double self = 8 * a * log((sqrt(2 * asq+ zsbsq) + a) / sqrt(asq + zsbsq)) - 
        2 * fabs(z + s * b) * (asin((asq * asq - zsbsq * zsbsq - s * asq * zsbsq) / (asq * asq + 2 * asq * zsbsq  + zsbsq * zsbsq))) +
        PI/2;
    s = -1;
    self = self + 
        8 * a * log((sqrt(2 * asq+ zsbsq) + a) / sqrt(asq + zsbsq)) - 
        2 * fabs(z + s * b) * (asin((asq * asq - zsbsq * zsbsq - s * asq * zsbsq) / (asq * asq + 2 * asq * zsbsq  + zsbsq * zsbsq)) +
        PI/2);
    
    return self;
}

void energy::valleau::update_potential(){
    double dz = 0.1;
    double diffz = 0;
    int numOfBins = chargeVector.size();
    ext.resize(numOfBins);
    for(int i = 0; i < numOfBins; i++){
        for(int j = 0; j < numOfBins; j++){
            diffz = fabs(j * dz - i * dz);
            ext[i] += chargeVector[j] * (-2 * PI * diffz - phiw(diffz));
        }
        ext[i] *= Base::lB * dz;
    }
    std::cout << ext << std::endl;
}

double energy::valleau::get_energy(Particle *particle){
    double energy = ext[(int)particle->pos[2] / 0.1];
    return energy;
}