#include "direct.h"
//#include "particle.h"
// Direct::Direct(){

// }

    // double get_replicates(Particle **particles, Particle *p);
    // double get_replicates(Particle **particles);
    // double get_central(Particle **particles);
    // double get_central(Particle **particles, Particle *p);

template<typename T>
double energy::direct::norm(T vec){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

double energy::direct::get_energy(Particle **particles){
    double energy;
    double central = get_central(particles);
    double replicates = 1.0/2.0 * get_replicates(particles);

    //printf("Central box: %lf Replicates: %lf\n", central, replicates);
    return Base::lB * (replicates + central);
}

// double energy::direct::get_energy(Particle **particles, Particle *p){
//     double energy;
//     double central = get_central(particles, p);
//     double replicates = get_replicates(particles, p);

//     //printf("Central box: %lf Replicates: %lf\n", central, replicates);
//     return Base::lB * (replicates + central);
// }

double energy::direct::get_replicates(Particle **particles){
    double energy = 0;
    double dist = 0;
    int rep = 0;
    int mx = rep;
    int my = rep;
    int mz = rep;
    int count = 0;
    int numOfRep = 0;

    //printf("Calculating energy for %d replicas.\n", (2 * rep+1) * (2 * rep+1) * (2 * rep+1) - 1);
    //#pragma omp parallel for if(rep > 10) reduction(+:energy) private(dist)
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
        //printf("Done: %lf\r", (double)count/(2 * mx + 1) * 100.0);
        //fflush(stdout);
    }
    //printf("\n");
    //printf("Number of replicas: %d\n", numOfRep - 1);
    return energy;
}

double energy::direct::get_central(Particle **particles){
    int k = 0;
    double energy = 0;
    double dist = 0;

    for(int i = 0; i < Particle::numOfParticles; i++){
        k = i + 1;
        while(k < Particle::numOfParticles){
            // double den[] = {particles[i]->pos[0] - particles[k]->pos[0], 
            //                 particles[i]->pos[1] - particles[k]->pos[1], 
            //                 particles[i]->pos[2] - particles[k]->pos[2]};
            // dist = norm(den);
            dist = Particle::distances[i][k];
            //dist = sqrt(particles[i]->distance(particles[k]));
            energy += particles[i]->q * particles[k]->q * 1/dist;
            k++;
        }  
    }

    return energy;
}

double energy::direct::get_central(Particle **particles, Particle *p){
    int k = 0;
    double energy = 0;
    //double lenergy = 0;
    double dist = 0;
    //#pragma omp parallel for reduction(+:energy) private(lenergy, dist)
    //{
    for(int i = 0; i < p->index; i++){
        dist = Particle::distances[i][p->index];
        energy += particles[i]->q * p->q * 1/dist;
    }
    //}
    //#pragma omp parallel for reduction(+:energy) private(lenergy, dist)
    //{
    for(int i = p->index + 1; i < Particle::numOfParticles; i++){
        dist = Particle::distances[p->index][i];
        energy += particles[i]->q * p->q * 1/dist;
    }
    //}
    return energy;
}   

double energy::direct::get_replicates(Particle **particles, Particle *p){
    double energy = 0;
    //double lenergy = 0;
    double dist = 0;
    //Eigen::Vector3d disp;
    int rep = 0;
    int mx = rep;
    int my = rep;
    int mz = rep;

    //printf("Calculating energy for %d replicas.\n", (2 * rep+1) * (2 * rep+1) * (2 * rep+1) - 1);
    //#pragma omp parallel for if(rep > 10) reduction(+:energy) private(dist)
    for(int i = -mx; i <= mx; i++){
        for(int j = -my; j <= my; j++){
            for(int k = -mz; k <= mz; k++){
                if(sqrt(i * i + j * j + k * k) <= rep){
                    //numOfRep++;
                    for(int l = 0; l < Particle::numOfParticles; l++){
                        if(i == 0 && j == 0 && k == 0){}
                        else{
                            double den[] = {p->pos[0] - particles[l]->pos[0] + i * Base::xL, 
                                            p->pos[1] - particles[l]->pos[1] + j * Base::yL , 
                                            p->pos[2] - particles[l]->pos[2] + k * Base::zL};
                            dist = norm(den);
                            //if(dist <= (rep + 1) * Base::xL){//(dist >= Base::xL && dist <= rep * Base::xL){
                            energy += p->q * particles[l]->q * 1/dist;
                            //}
                        }
                    }
                }
            }
        }
        //printf("Done: %lf\r", (double)count/(2 * mx + 1) * 100.0);
        //fflush(stdout);
    }
    //energy += lenergy;
    //printf("Energy: %lf\n", energy);
    //printf("\n");
    //printf("Number of replicas: %d\n", numOfRep - 1);
    return energy;
}

