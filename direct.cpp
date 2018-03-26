#include "direct.h"

Direct::Direct(){

}

template<typename T>
double Direct::norm(T vec){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

double Direct::get_energy(Particle **particles){
    double energy;
    int mc = 0;

    for(int i = -mc; i <= mc; i++){
        for(int j = -mc; j <= mc; j++){
            for(int k = -mc; k <= mc; k++){
                for(int l = 0; l < Particle::numOfParticles; l++){
                    for(int m = 0; m < Particle::numOfParticles; m++){
                        if(i == 0 && j == 0 && k == 0 && m == l){}
                        else{
                            
                            double den[] = {particles[m]->pos[0] - particles[l]->pos[0] + i * Base::xL, 
                                            particles[m]->pos[1] - particles[l]->pos[1] + j * Base::yL, 
                                            particles[m]->pos[2] - particles[l]->pos[2] + k * Base::zL};

                            energy += particles[m]->q * particles[l]->q * 1/norm(den);
                        }
                    }
                        
                }
            }
        }
    }
    return 1.0/2.0 * energy;
}