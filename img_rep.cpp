#include "img_rep.h"

void energy::imgrep::set_positions(Particle **particles){
    int numOfParticles = Particle::numOfParticles;
    for(int i = 0; i < numOfParticles; i++){
        //Second left reflection
        particles[i + numOfParticles] = new Particle();
        particles[i + numOfParticles]->pos = particles[i]->pos;
        particles[i + numOfParticles]->com = particles[i]->com;
        particles[i + numOfParticles]->chargeDisp[2] = particles[i]->chargeDisp[2];
        particles[i + numOfParticles]->q = particles[i]->q;
        strcpy(particles[i]->name, "DUM\0");
        particles[i + numOfParticles]->pos[2] = -2 * (Base::zL - 2 * Base::wall) + particles[i]->pos[2];
        //First left reflection
        particles[i + 2 * numOfParticles] = new Particle();
        particles[i + 2 * numOfParticles]->pos = particles[i]->pos;
        particles[i + 2 * numOfParticles]->com = particles[i]->com;
        particles[i + 2 * numOfParticles] = particles[i];
        particles[i + 2 * numOfParticles]->chargeDisp[2] = -1 * particles[i]->chargeDisp[2];
        particles[i + 2 * numOfParticles]->q = -1 * particles[i]->q;
        strcpy(particles[i + 2 * numOfParticles]->name, "DUM\0");
        particles[i + 2 * numOfParticles]->pos[2] = -1 * particles[i]->pos[2];

        //First right reflection
        particles[i + 3 * numOfParticles] = new Particle();
        particles[i + 3 * numOfParticles]->pos = particles[i]->pos;
        particles[i + 3 * numOfParticles]->com = particles[i]->com;
        particles[i + 3 * numOfParticles]->chargeDisp[2] = -1 * particles[i]->chargeDisp[2];
        particles[i + 3 * numOfParticles]->q = -1 * particles[i]->q;
        strcpy(particles[i + 3 * numOfParticles]->name, "DUM\0");
        particles[i + 3 * numOfParticles]->pos[2] = -1 * particles[i]->pos[2] + 2 * (Base::zL - 2 * Base::wall);
        //Second right reflection
        particles[i + 4 * numOfParticles] = new Particle();
        particles[i + 4 * numOfParticles]->pos = particles[i]->pos;
        particles[i + 4 * numOfParticles]->com = particles[i]->com;
        particles[i + 4 * numOfParticles]->chargeDisp[2] = particles[i]->chargeDisp[2];
        particles[i + 4 * numOfParticles]->q = particles[i]->q;
        strcpy(particles[i + 4 * numOfParticles]->name, "DUM\0");
        particles[i + 4 * numOfParticles]->pos[2] = particles[i]->pos[2] + 2 * (Base::zL - 2 * Base::wall);

        printf("%d  %d  %d  %d\n", i + numOfParticles, i + 2 * numOfParticles, i + 3 * numOfParticles, i + 4 * numOfParticles);
    }
}

void energy::imgrep::update_position(Particle **particles, Particle *p){
    int numOfParticles = Particle::numOfParticles / 5;
    //Second left reflection
    particles[p->index + numOfParticles]->pos[2] = -2 * (Base::zL - 2 * Base::wall) + p->pos[2];
    //First left reflection
    particles[p->index + 2 * numOfParticles]->pos[2] = -1 * p->pos[2];
    //First right reflection
    particles[p->index + 3 * numOfParticles]->pos[2] = -1 * p->pos[2] + 2 * (Base::zL - 2 * Base::wall);
    //Second right reflection
    particles[p->index + 4 * numOfParticles]->pos[2] = p->pos[2] + 2 * (Base::zL - 2 * Base::wall);
}

double energy::imgrep::get_particle_energy(Particle **particles, Particle *p){
    return 5;
}

double energy::imgrep::get_energy(Particle **particles){
    return 5;
}