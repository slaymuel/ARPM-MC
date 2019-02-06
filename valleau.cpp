#include "valleau.h"
double energy::valleau::binWidth; //Angstrom
int numOfBins; //number of bins
Eigen::VectorXd energy::valleau::pDensity;
Eigen::VectorXd energy::valleau::nDensity;
Eigen::VectorXd energy::valleau::chargeVector;
Eigen::VectorXd energy::valleau::ext;
Eigen::VectorXd energy::valleau::imgExt;
int energy::valleau::numOfSamples;


void energy::valleau::initialize(){
    energy::valleau::binWidth = 0.05;
    numOfBins = (Base::zL - 2 * Base::wall)/energy::valleau::binWidth;
    energy::valleau::pDensity.resize(numOfBins);
    energy::valleau::nDensity.resize(numOfBins);
    energy::valleau::ext.resize(numOfBins);
    energy::valleau::imgExt.resize(numOfBins);
    pDensity.setZero();
    nDensity.setZero();
    ext.setZero();
    imgExt.setZero();
    energy::valleau::numOfSamples = 0;   
}


void energy::valleau::update_charge_vector(Particle **particles){
    //Eigen::VectorXd pDensity(numOfBins);
    //Eigen::VectorXd nDensity(numOfBins);
    
    for(int i = 0; i < Particle::numOfParticles; i++){
        if(particles[i]->q > 0){
            pDensity[(int)((particles[i]->pos[2] - Base::wall)/binWidth)]++;
        }
        else{
            nDensity[(int)((particles[i]->pos[2] - Base::wall)/binWidth)]++;
        }
    }
    numOfSamples++;
}


double energy::valleau::phiw(double z){
    double a = Base::xL/2.0;
    double asq = a * a;
    double zsq = z * z;//(z + b) * (z + b);

    double self = 8.0 * a * std::log((std::sqrt(2.0 * asq + zsq) + a) / std::sqrt(asq + zsq)) - 
                  2.0 * z * (std::asin((asq * asq - zsq * zsq - 2.0 * asq * zsq) / std::pow(asq + zsq, 2.0)) + PI/2.0);
    return self;
}


void energy::valleau::update_potential(){
    chargeVector = pDensity - nDensity;
    chargeVector = chargeVector/(Base::xL * Base::yL * binWidth * numOfSamples);

    //Symmetrize charge vector
    int j = chargeVector.size() - 1;
    double avg;
    for(int i = 0; i < (int)(chargeVector.size()/2); i++){
        avg = (chargeVector[i] + chargeVector[j]) / 2.0;
        chargeVector[i] = avg;
        chargeVector[j] = avg;
        j--;
    }

    //Intergrate:
    double dz = 0.05;
    double diffz = 0;
    double iIt = 0;
    double jIt = 0;
    double dhs = 2.5; //hard-sphere radius
    int numOfBins = chargeVector.size();

    ext.setZero();
    imgExt.setZero();
    iIt = -0.5 * dz;
    for(int i = 0; i < numOfBins; i++){
        iIt += dz;
        jIt = -0.5 * dz;

        for(int j = 0; j < numOfBins; j++){
            jIt += dz;

            //xy replicates of box
            //diffz = std::fabs(jIt - iIt);
            //ext[i] += chargeVector[j] * (-2.0 * PI * diffz - phiw(diffz));

            //First reflection
            diffz = jIt + iIt - dhs;
            imgExt[i] += chargeVector[j] * (-2.0 * PI * diffz - phiw(diffz));
            diffz = (Base::zL - 2 * Base::wall) - jIt + ((Base::zL - 2 * Base::wall) - iIt) - dhs;
            imgExt[i] += chargeVector[j] * (-2.0 * PI * diffz - phiw(diffz));
            //Second reflection
            diffz = 2 * (Base::zL - 2 * Base::wall) + jIt - iIt - 2 * dhs;
            imgExt[i] -= chargeVector[j] * (-2.0 * PI * diffz - phiw(diffz));
            diffz = 2 * (Base::zL - 2 * Base::wall) - jIt + iIt - 2 * dhs;
            imgExt[i] -= chargeVector[j] * (-2.0 * PI * diffz - phiw(diffz));
            
        }
    }
    ext = dz * Base::lB * (ext - imgExt);

    pDensity.setZero();
    nDensity.setZero();
    numOfSamples = 0;

    char histo_name[64] = "ext.txt";
    FILE *f = fopen(histo_name, "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    for(int i = 0; i < numOfBins; i++){
        fprintf(f, "%lf     %.12lf\n", (double)i * binWidth, ext[i]);
    }
    fclose(f);
}


double energy::valleau::get_images(Particle **particles){
    int numOfReflections = 2;
    double energy = 0;
    double distance = 0;
    Eigen::Vector3d disp;
    double tempX = 0;
    double tempY = 0;

    for(int k = 1; k <= numOfReflections; k += 2){  //Uneven reflections of opposite sign
        for(int i = 0; i < Particle::numOfParticles; i++){
            for(int j = i; j < Particle::numOfParticles; j++){
                tempX = particles[j]->pos[0];
                tempY = particles[j]->pos[1];

                if(particles[i]->pos[0] - particles[j]->pos[0] < -1.0 * Base::xL / 2.0){
                    tempX = particles[j]->pos[0] - Base::xL;
                }
                if(particles[i]->pos[0] - particles[j]->pos[0] > Base::xL / 2.0){
                    tempX = particles[j]->pos[0] + Base::xL;
                }
                if(particles[i]->pos[1] - particles[j]->pos[1] < -1 * Base::yL / 2.0){
                    tempY = particles[j]->pos[1] - Base::yL;
                }
                if(particles[i]->pos[1] - particles[j]->pos[1] > Base::yL / 2.0){
                    tempY = particles[j]->pos[1] + Base::yL;
                }

                //Left reflections
                disp << particles[i]->pos[0] - tempX,
                        particles[i]->pos[1] - tempY,
                        (particles[i]->pos[2] - Base::wall) + (k - 1) * (Base::zL - 2.0 * Base::wall) + (particles[j]->pos[2] - Base::wall);
                distance = disp.norm();
                if(i == j){
                    energy -= particles[i]->q * particles[j]->q * 1.0 /(distance * 2.0);
                }
                else{
                    energy -= particles[i]->q * particles[j]->q * 1.0 /distance;
                }

                //Right reflections
                disp << particles[i]->pos[0] - tempX,
                        particles[i]->pos[1] - tempY,
                        (k + 1) * (Base::zL - 2 * Base::wall) - (particles[i]->pos[2] - Base::wall) - (particles[j]->pos[2] - Base::wall);
                distance = disp.norm();
                if(i == j){
                    energy -= particles[i]->q * particles[j]->q * 1.0 / (distance * 2.0);
                }
                else{
                    energy -= particles[i]->q * particles[j]->q * 1.0 / distance;
                }
            }
        }
    }

    for(int m = 2; m <= numOfReflections; m += 2){  //Even reflections of same sign
        for(int i = 0; i < Particle::numOfParticles; i++){
            for(int j = i; j < Particle::numOfParticles; j++){
                tempX = particles[j]->pos[0];
                tempY = particles[j]->pos[1];

                if(particles[i]->pos[0] - particles[j]->pos[0] < -1.0 * Base::xL/2.0){
                    tempX = particles[j]->pos[0] - Base::xL;
                }
                if(particles[i]->pos[0] - particles[j]->pos[0] > Base::xL / 2.0){
                    tempX = particles[j]->pos[0] + Base::xL;
                }
                if(particles[i]->pos[1] - particles[j]->pos[1] < -1.0 * Base::yL / 2.0){
                    tempY = particles[j]->pos[1] - Base::yL;
                }
                if(particles[i]->pos[1] - particles[j]->pos[1] > Base::yL / 2.0){
                    tempY = particles[j]->pos[1] + Base::yL;
                }
                //Left reflections
                disp << particles[i]->pos[0] - tempX,
                        particles[i]->pos[1] - tempY,
                        m * (Base::zL - 2 * Base::wall) + (particles[i]->pos[2] - Base::wall) - (particles[j]->pos[2] - Base::wall);
                distance = disp.norm();
                if(i == j){
                    energy += particles[i]->q * particles[j]->q * 1.0 / (distance * 2.0);
                }
                else{
                    energy += particles[i]->q * particles[j]->q * 1.0 / distance;
                }

                //Right reflections
                disp << particles[i]->pos[0] - tempX,
                        particles[i]->pos[1] - tempY,
                        m * (Base::zL - 2 * Base::wall) + (particles[j]->pos[2] - Base::wall) - (particles[i]->pos[2] - Base::wall);
                distance = disp.norm();
                if(i == j){
                    energy += particles[i]->q * particles[j]->q * 1.0 / (distance * 2.0);
                }
                else{
                    energy += particles[i]->q * particles[j]->q * 1.0 / distance;
                }
            }
        }
    }
    energy *= Base::lB;
    return energy;
}


double energy::valleau::get_particle_images(Particle **particles, Particle *p){
    int numOfReflections = 2;
    double energy = 0;
    double distance = 0;
    Eigen::Vector3d disp;
    double tempX = 0;
    double tempY = 0;
    for(int k = 1; k <= numOfReflections; k += 2){  //Uneven reflections of opposite sign
        for(int j = 0; j < Particle::numOfParticles; j++){
            tempX = particles[j]->pos[0];
            tempY = particles[j]->pos[1];

            if(p->pos[0] - particles[j]->pos[0] < -1.0 * Base::xL/2.0){
                tempX = particles[j]->pos[0] - Base::xL;
            }
            if(p->pos[0] - particles[j]->pos[0] > Base::xL/2){
                tempX = particles[j]->pos[0] + Base::xL;
            }
            if(p->pos[1] - particles[j]->pos[1] < -1 * Base::yL/2){
                tempY = particles[j]->pos[1] - Base::yL;
            }
            if(p->pos[1] - particles[j]->pos[1] > Base::yL/2){
                tempY = particles[j]->pos[1] + Base::yL;
            }

            //Left reflections
            disp << p->pos[0] - tempX,
                    p->pos[1] - tempY,
                    (p->pos[2] - Base::wall) + (k - 1) * (Base::zL - 2 * Base::wall) + (particles[j]->pos[2] - Base::wall);
            distance = disp.norm();
            if(p->index == j){
                energy -= p->q * particles[j]->q * 1/(distance * 2);
            }
            else{
                energy -= p->q * particles[j]->q * 1/distance;
            }

            //Right reflections
            disp << p->pos[0] - tempX,
                    p->pos[1] - tempY,
                    (k + 1) * (Base::zL - 2 * Base::wall) - (p->pos[2] - Base::wall) - (particles[j]->pos[2] - Base::wall);
            distance = disp.norm();
            if(p->index == j){
                energy -= p->q * particles[j]->q * 1.0/(distance * 2.0);
            }
            else{
                energy -= p->q * particles[j]->q * 1.0/distance;
            }
        }
    }

    for(int m = 2; m <= numOfReflections; m += 2){  //Even reflections of same sign
        for(int j = 0; j < Particle::numOfParticles; j++){
            tempX = particles[j]->pos[0];
            tempY = particles[j]->pos[1];

            if(p->pos[0] - particles[j]->pos[0] < -1.0 * Base::xL/2.0){
                tempX = particles[j]->pos[0] - Base::xL;
            }
            if(p->pos[0] - particles[j]->pos[0] > Base::xL/2){
                tempX = particles[j]->pos[0] + Base::xL;
            }
            if(p->pos[1] - particles[j]->pos[1] < -1 * Base::yL/2){
                tempY = particles[j]->pos[1] - Base::yL;
            }
            if(p->pos[1] - particles[j]->pos[1] > Base::yL/2){
                tempY = particles[j]->pos[1] + Base::yL;
            }
            //Left reflections
            disp << p->pos[0] - tempX,
                    p->pos[1] - tempY,
                    m * (Base::zL - 2 * Base::wall) + (p->pos[2] - Base::wall) - (particles[j]->pos[2] - Base::wall);
            distance = disp.norm();
            if(p->index == j){
                energy += p->q * particles[j]->q * 1.0/(distance * 2.0);
            }
            else{
                energy += p->q * particles[j]->q * 1.0/distance;
            }

            //Right reflections
            disp << p->pos[0] - tempX,
                    p->pos[1] - tempY,
                    m * (Base::zL - 2 * Base::wall) + (particles[j]->pos[2] - Base::wall) - (p->pos[2] - Base::wall);
            distance = disp.norm();
            if(p->index == j){
                energy += p->q * particles[j]->q * 1.0/(distance * 2.0);
            }
            else{
                energy += p->q * particles[j]->q * 1.0/distance;
            }
        }
    }
    energy *= Base::lB;
    return energy;
}


double energy::valleau::get_energy(Particle **particles){
    double energy = 0;
    double y0, y1, x0, x1;
    double dz = 0.05;
    
    for(int i = 0; i < Particle::numOfParticles; i++){
        //Linear interpolation
        x0 = (int)(particles[i]->pos[2] - Base::wall) / dz;
        x1 = (int)(particles[i]->pos[2] - Base::wall) / dz + 1;
        y0 = ext[(int)(particles[i]->pos[2] - Base::wall) / dz];
        y1 = ext[(int)(particles[i]->pos[2] - Base::wall) / dz + 1];
        if(particles[i]->q > 0){
            energy += y0 + ((particles[i]->pos[2] - Base::wall) / dz - x0) * (y1 - y0) / (x1 - x0);
        }
        else{
            energy -= y0 + ((particles[i]->pos[2] - Base::wall) / dz - x0) * (y1 - y0) / (x1 - x0);
        }
    }
    //printf("Valleau dir: %lf    pol: %lf\n", energy::direct::get_energy(particles), get_images(particles));
    //energy += energy::direct::get_energy(particles);
    energy += energy::ewald3D::get_energy(particles);
    energy += get_images(particles);
    return energy;
}


double energy::valleau::get_particle_energy(Particle **particles, Particle *p){
    double energy = 0;
    double y0, y1, x0, x1;
    double dz = 0.05;
    
    //Linear interpolation
    x0 = (int)(p->pos[2] - Base::wall) / dz;
    x1 = (int)(p->pos[2] - Base::wall) / dz + 1;
    y0 = ext[(int)(p->pos[2] - Base::wall) / dz];
    y1 = ext[(int)(p->pos[2] - Base::wall) / dz + 1];
    if(p->q > 0){
        energy += y0 + ((p->pos[2] - Base::wall) / dz - x0) * (y1 - y0) / (x1 - x0);
    }
    else{
        energy -= y0 + ((p->pos[2] - Base::wall) / dz - x0) * (y1 - y0) / (x1 - x0);
    }
    //printf("dir: %lf    pol: %lf\n", energy::direct::get_particle_energy(particles, p), get_particle_images(particles, p));
    //energy += energy::direct::get_particle_energy(particles, p);
    energy += energy::ewald3D::get_particle_energy(particles, p);
    energy += get_particle_images(particles, p);
    
    return energy;
}
