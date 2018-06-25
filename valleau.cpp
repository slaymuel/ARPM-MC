#include "valleau.h"
double binWidth = 0.05; //Angstrom
int numOfBins = Base::zL/binWidth; //number of bins
Eigen::VectorXd energy::valleau::pDensity(numOfBins + 1);
Eigen::VectorXd energy::valleau::nDensity(numOfBins + 1);
int energy::valleau::numOfSamples = 0;

void energy::valleau::update_charge_vector(Particle **particles){
    //Eigen::VectorXd pDensity(numOfBins);
    //Eigen::VectorXd nDensity(numOfBins);

    for(int i = 0; i < Particle::numOfParticles; i++){
        if(particles[i]->q > 0){
            pDensity[(int)(particles[i]->pos[2]/binWidth)]++;
        }
        else{
            nDensity[(int)(particles[i]->pos[2]/binWidth)]++;
        }
    }
    numOfSamples++;
}


double energy::valleau::phiw(double z){
    double a = Base::xL/2;
    double asq = a * a;

    double zsbsq = z * z;//(z + b) * (z + b);
    double self = 8 * a * log((sqrt(2 * asq + zsbsq) + a) / sqrt(asq + zsbsq)) - 
                  2 * z * (asin((asq * asq - zsbsq * zsbsq - 2 * asq * zsbsq) / (asq * asq + 2 * asq * zsbsq  + zsbsq * zsbsq)) + PI/2);
        /*
    s = -1;
    self = self + 
        8 * a * log((sqrt(2 * asq+ zsbsq) + a) / sqrt(asq + zsbsq)) - 
        2 * fabs(z + s * b) * (asin((asq * asq - zsbsq * zsbsq - s * asq * zsbsq) / (asq * asq + 2 * asq * zsbsq  + zsbsq * zsbsq)) +
        PI/2);
    */
    return self;
}


void energy::valleau::update_potential(){
    chargeVector = pDensity - nDensity; //charge vector
    chargeVector = chargeVector/(Base::xL * Base::zL * binWidth * numOfSamples);

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

    //Intergrate:
    double dz = 0.05;
    double diffz = 0;
    double iIt = 0;
    double jIt = 0;
    int numOfBins = chargeVector.size();

    ext.resize(numOfBins);
    ext.setZero();
    for(int i = 0; i < numOfBins; i++){
        iIt = 0.5 * dz + i * dz;
        for(int j = 0; j < numOfBins; j++){
            jIt = 0.5 * dz + j * dz;
            //diffz = fabs(j * dz - i * dz);
            diffz = fabs(jIt - iIt);
            ext[i] += chargeVector[j] * (-2.0 * PI * diffz - phiw(diffz));
        }
        ext[i] *= Base::lB * dz;
    }
    printf("External potential:\n");
    std::cout << ext << std::endl;
}

double energy::valleau::get_images(Particle **particles){

    return 0.0;
}

double energy::valleau::get_energy(Particle **particles){
    double energy = 0;
    double y0, y1, x0, x1;
    double dz = 0.05;
    for(int i = 0; i < Particle::numOfParticles; i++){
        //Linear interpolation
        x0 = (int)particles[i]->pos[2] / dz;
        x1 = (int)particles[i]->pos[2] / dz + 1;
        y0 = ext[(int)particles[i]->pos[2] / dz];
        y1 = ext[(int)particles[i]->pos[2] / dz + 1];
        if(particles[i]->q > 0){
            energy += y0 + (particles[i]->pos[2] / dz - x0) * (y1 - y0) / (x1 - x0);
            //energy += ext[(int)particles[i]->pos[2] / 0.1];
        }
        else{
            energy -= y0 + (particles[i]->pos[2] / dz - x0) * (y1 - y0) / (x1 - x0);
            //energy -= ext[(int)particles[i]->pos[2] / 0.1];
        }
    }
    //energy *= Base::lB;
    energy += energy::direct::get_energy(particles);
    return energy;
}