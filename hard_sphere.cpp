#include "hard_sphere.h"

double energy::hs::get_energy(Particle **particles){
    for(int i = 0; i < Particle::numOfParticles; i++){
        for(int j = i + 1; j < Particle::numOfParticles; j++){
            if(Particle::distances[i][j] < particles[i]->d/2 + particles[j]->d/2)
            {
                return 10000000000000.0;
            }
        }
    }
    return 0;
}

double energy::hs::get_particle_energy(Particle **particles, Particle *p){
    for(int i = 0; i < p->index; i++){
        if(Particle::distances[i][p->index] < particles[i]->d/2 + p->d/2)
        {
            return 100000000000000.0;
        }
    }
    for(int i = p->index + 1; i < Particle::numOfParticles; i++){
        if(Particle::distances[p->index][i] < particles[i]->d/2 + p->d/2)
        {
            return 10000000000000.0;
        }
    }
    return 0;
}