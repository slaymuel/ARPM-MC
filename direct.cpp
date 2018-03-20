#include "direct.h"

Direct::Direct(){

}

template<typename T>
double Direct::norm(T vec){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

double Direct::get_energy(Particle *p, Particle **particles){
    double energy;
    int mc = 40;

    for(int i = -mc; i <= mc; i++){
        //printf("i: %d\n", i);
        for(int j = -mc; j <= mc; j++){
            for(int k = -mc; k <= mc; k++){
                for(int l = 0; l < Particle::numOfParticles; l++){
                    if(i == 0 && j == 0 && k == 0 && p->index == l){
                        continue;
                    }
                    double den[] = {p->pos[0] - particles[l]->pos[0] + i * Base::xL, 
                                    p->pos[1] - particles[l]->pos[1] + j * Base::yL, 
                                    p->pos[2] - particles[l]->pos[2] + k * Base::zL};

                    energy += p->q * particles[l]->q * 1/norm(den);
                    //printf("norm: %lf\n", norm(den));
                    //energy += pow(EC, 2)/(1e-10) * (p->q * particles[l]->q)/norm(den);
                }
            }
        }
    }
    return energy;
}