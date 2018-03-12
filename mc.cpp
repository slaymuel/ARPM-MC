#include "particle.h"
#include "base.cpp"

Particle* pbc(Particle *p){
    //Translate particles according to periodic boundary conditions
    if(p->pos[0] > Base::box[0]){
        p->pos[0] = p->pos[0] - Base::box[0];
    }
    if(p->pos[0] < 0){
        p->pos[0] = p->pos[0] + Base::box[0];
    }
    if(p->pos[1] > Base::box[1]){
        p->pos[1] = p->pos[1] - Base::box[1];
    }
    if(p->pos[1] < 0){
        p->pos[1] = p->pos[1] + Base::box[1];
    }
    if(p->pos[2] < 0){
        p->pos[2] = p->pos[2] + Base::box[2];
    }
    if(p->pos[2] > Base::box[2]){
        p->pos[2] = p->pos[2] - Base::box[2];
    }

    return p;
}