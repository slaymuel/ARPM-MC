#include "direct.h"

Direct::Direct(){

}

template<typename T>
double Direct::norm(T vec){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

double Direct::get_energy(Particle **particles){
    double energy;
    double central = get_central(particles);
    double replicates = get_replicates(particles);

    printf("Central box: %lf Replicates: %lf\n", central, replicates);
    return replicates + central;
}

double Direct::get_replicates(Particle **particles){
    double energy = 0;
    double dist = 0;
    int rep = 8;
    int mx = rep;
    int my = rep;
    int mz = rep;
    int count = 0;
    int numOfRep = 0;

    printf("Calculating energy for %d replicas.\n", (2 * rep+1) * (2 * rep+1) * (2 * rep+1) - 1);

    for(int i = -mx; i <= mx; i++){
        for(int j = -my; j <= my; j++){
            for(int k = -mz; k <= mz; k++){
                if(sqrt(i * i + j * j + k * k) <= rep){
                    numOfRep++;
                    for(int l = 0; l < Particle::numOfParticles; l++){
                        for(int m = 0; m < Particle::numOfParticles; m++){
                            if(i == 0 && j == 0 && k == 0){}
                            else{
                                double den[] = {particles[l]->pos[0] - particles[m]->pos[0] + i * Base::xL, 
                                                particles[l]->pos[1] - particles[m]->pos[1] + j * Base::yL, 
                                                particles[l]->pos[2] - particles[m]->pos[2] + k * Base::zL};
                                dist = norm(den);
                                //if(dist <= (rep + 1) * Base::xL){//(dist >= Base::xL && dist <= rep * Base::xL){
                                energy += particles[m]->q * particles[l]->q * 1/dist;
                                //}
                            }
                        }
                    }
                }
            }
        }
        count++;
        printf("Done: %lf\r", (double)count/(2 * mx + 1) * 100.0);
        fflush(stdout);
    }
    printf("\n");
    printf("Number of replicas: %d\n", numOfRep);
    return energy;
}

double Direct::get_central(Particle **particles){
    int k = 0;
    double energy = 0;
    double dist = 0;

    for(int i = 0; i < Particle::numOfParticles; i++){
        k = i + 1;
        while(k < Particle::numOfParticles){
            if(i != k){
                double den[] = {particles[i]->pos[0] - particles[k]->pos[0], 
                                particles[i]->pos[1] - particles[k]->pos[1], 
                                particles[i]->pos[2] - particles[k]->pos[2]};
                dist = norm(den);
                //dist = sqrt(particles[i]->distance(particles[k]));
                energy += particles[i]->q * particles[k]->q * 1/dist;
            }
            k++;
        }  
    }

    return energy;
}   
