#include "particle.h"
#include <iostream>
#include <math.h>

//Constructor
Particle::Particle(){
    //Keep track of the number of particles
    numOfParticles++;
}

void Particle::pbc(){
    //Translate particles according to periodic boundary conditions
    if(this->pos[0] > xL){
        this->pos[0] = this->pos[0] - xL;
    }
    if(this->pos[0] < 0){
        this->pos[0] = this->pos[0] + xL;
    }
    if(this->pos[1] > yL){
        this->pos[1] = this->pos[1] - yL;
    }
    if(this->pos[1] < 0){
        this->pos[1] = this->pos[1] + yL;
    }
    if(this->pos[2] > zL){
        this->pos[2] = this->pos[2] - zL;
    }
    if(this->pos[2] < 0){
        this->pos[2] = this->pos[2] + zL;
    }
}

void Particle::pbc_xy(){
    //Translate particles according to periodic boundary conditions in the xy-directions
    if(this->pos[0] > xL){
        this->pos[0] = this->pos[0] - xL;
    }
    if(this->pos[0] < 0){
        this->pos[0] = this->pos[0] + xL;
    }
    if(this->pos[1] > yL){
        this->pos[1] = this->pos[1] - yL;
    }
    if(this->pos[1] < 0){
        this->pos[1] = this->pos[1] + yL;
    }
}

// void Particle::randomMove(){
//     this->pos[0] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pos[1] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pos[2] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pbc();
// }

void Particle::randomMove_xy(double stepSize){
    this->pos[0] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;
    this->pos[1] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;
    //printf("%lf\n", ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize);
    if(this->pos[2] < zL - (this->d/2 + stepSize) && this->pos[2] > this->d/2 + stepSize){
        this->pos[2] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;
    }
    this->pbc_xy();
}

void Particle::randomMove(double stepSize){
    this->pos[0] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;
    this->pos[1] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;
    this->pos[2] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;

    this->pbc();
}

double Particle::distance(Particle *p){
    //Calculate distance between particles
    double xP1 = p->pos[0];
    double yP1 = p->pos[1];
    double zP1 = p->pos[2];
    double xP2 = this->pos[0];
    double yP2 = this->pos[1];
    double zP2 = this->pos[2];

    if(xP1 - xP2 < -1 * xL/2){
        xP2 = xP2 - xL;
    }
    if(xP1 - xP2 > xL/2){
        xP2 = xP2 + xL;
    }
    if(yP1 - yP2 < -1 * yL/2){
        yP2 = yP2 - yL;
    }
    if(yP1 - yP2 > yL/2){
        yP2 = yP2 + yL;
    }
    if(zP1 - zP2 < -1 * zL/2){
        zP2 = zP2 - zL;
    }
    if(zP1 - zP2 > zL/2){
        zP2 = zP2 + zL;
    }
    return (pow((xP1 - xP2), 2) + pow((yP1 - yP2), 2) + pow((zP1 - zP2), 2));
}

double Particle::distance_xy(Particle *p){
    //Calculate distance between particles
    double xP1 = p->pos[0];
    double yP1 = p->pos[1];
    double zP1 = p->pos[2];
    double xP2 = this->pos[0];
    double yP2 = this->pos[1];
    double zP2 = this->pos[2];

    if(xP1 - xP2 < -1 * xL/2){
        xP2 = xP2 - xL;
    }
    if(xP1 - xP2 > xL/2){
        xP2 = xP2 + xL;
    }
    if(yP1 - yP2 < -1 * yL/2){
        yP2 = yP2 - yL;
    }
    if(yP1 - yP2 > yL/2){
        yP2 = yP2 + yL;
    }
    return (pow((xP1 - xP2), 2) + pow((yP1 - yP2), 2) + pow((zP1 - zP2), 2));
}

int Particle::hardSphere(Particle **particles){
    int i = 0;
    int j = 0;

    for(i = 0; i < Particle::numOfParticles; i++){
        if(i != this->index){
            if(this->distance(particles[i]) < (this->d+particles[i]->d)/2){
                return 0;
            }
        }
    }
    return 1;
}

