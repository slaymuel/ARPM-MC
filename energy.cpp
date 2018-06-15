#include "particle.h"
#include "constants.h"

/*
double getEnergy(Particle **particles){
    int i = 0;
    int j = 0;
    double dist = 0;
    double energy = 0;
    double r6 = 0;
    double sigma = 1;
    double epsilon = EPSILON;

    for(i = 0; i < Particle::numOfParticles; i++){
        j = i + 1;
        while(j < Particle::numOfParticles){
            //LJ
            dist = distance(particles[j], particles[i]);
            r6 = pow(dist, 3);
            energy += 4 * (1/(r6*r6) - 1/r6);
            //Coloumb
            dist = sqrt(dist);
            energy += 1/(4 * PI * VP * 2 * T) * (particles[i]->q * particles[j]->q)/dist;
            j++;
        }
    }
    return energy;
}

double getParticleEnergy(int pInd, Particle *p, Particle **particles){
    int i = 0;
    double energy = 0;
    double dist = 0;
    double r6 = 0;
    double sigma = 1;
    double epsilon = EPSILON;

    for(i = 0; i < Particle::numOfParticles; i++){
        if(i != pInd){
            //LJ
            dist = distance(p, particles[i]);
            r6 = pow(dist, 3);
            energy += 4 * (1/(r6*r6) - 1/r6);
            //printf("Lennard: %lf\n", energy);
            //Coloumb
            dist = sqrt(dist);
            energy += 1/(4 * PI * VP * 2 * T) * (p->q * particles[i]->q)/dist;
            //printf("Coloumb: %lf\n", energy);
        }
    }
    return energy;
}*/
