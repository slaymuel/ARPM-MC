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

double energy::direct::get_energy(Particles &particles){
    double central = get_central(particles);
    double replicates = 1.0/2.0 * get_replicates(particles);

    //printf("Central box: %lf Replicates: %lf\n", central, replicates);
    //printf("Energy: %lf\n", central + replicates);
    return Base::lB * (central);
}


double energy::direct::get_particle_energy(Particles &particles, Particle &p){
    double central = get_central(particles, p);
    double replicates = 1.0/2.0 * get_replicates(particles);

    //printf("Central box: %lf Replicates: %lf\n", central, replicates);
    return Base::lB * (central);
}

double energy::direct::get_particle_pot(Particles &particles, Particle &p){
    double central = get_central_pot(particles, p);
    //printf("Central box: %lf Replicates: %lf\n", central, replicates);
    return Base::lB * (central);
}

// double energy::direct::get_energy(Particle **particles, Particle *p){
//     double energy;
//     double central = get_central(particles, p);
//     double replicates = get_replicates(particles, p);

//     //printf("Central box: %lf Replicates: %lf\n", central, replicates);
//     return Base::lB * (replicates + central);
// }

double energy::direct::get_replicates(Particles &particles){
    double energy = 0;
    double dist = 0;
    int rep = 0;
    int mx = rep;
    int my = rep;
    int mz = 0;
    int count = 0;

    //printf("Calculating energy for %d replicas.\n", (2 * rep+1) * (2 * rep+1) * (2 * rep+1) - 1);
    //#pragma omp parallel for if(rep > 10) reduction(+:energy) private(dist)
    for(int i = -mx; i <= mx; i++){
        for(int j = -my; j <= my; j++){
            for(int k = -mz; k <= mz; k++){
                if(sqrt(i * i + j * j + k * k) <= rep){
                    for(int l = 0; l < particles.numOfParticles; l++){
                        for(int m = 0; m < particles.numOfParticles; m++){
                            if(i == 0 && j == 0 && k == 0){}
                            else{
                                double den[] = {particles[l].pos[0] - particles[m].pos[0] + i * Base::xL, 
                                                particles[l].pos[1] - particles[m].pos[1] + j * Base::yL, 
                                                particles[l].pos[2] - particles[m].pos[2] + k * Base::zLBox};
                                dist = norm(den);
                                //if(dist <= (rep + 1) * Base::xL){//(dist >= Base::xL && dist <= rep * Base::xL){
                                energy += particles[m].q * particles[l].q * 1.0 / dist;
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

double energy::direct::get_central(Particles &particles){
    double energy = 0;
    double dist = 0;
    Eigen::Vector3d disp;
    for(int i = 0; i < particles.numOfParticles; i++){
        for(int k = i + 1; k < particles.numOfParticles; k++){
            //disp = particles[i]->pos - particles[k]->pos;
            //dist = disp.norm();
            //dist = particles.distances[i][k];

            dist = particles[i].distance_xy(particles[k]);
            energy += particles[i].q * particles[k].q / dist;
        }  
        //std::cout << particles[i].com << std::endl;
        //std::cout << particles[i].pos << std::endl;
        if(i != particles[i].index){
            printf("Index for particle %i is wrong, is has %i!!\n", i, particles[i].index);
        }
        if(particles[i].pos != particles[i].com){
            printf("com != pos for %i!!\n", i);
            std::cout << particles[i].com << std::endl;
            std::cout << particles[i].pos << std::endl;
        }
        if(particles[i].pos[2] < -Base::zLBox / 2.0 - 2.0 + particles[i].d / 2.0 || particles[i].pos[1] < -Base::yL / 2.0 || particles[i].pos[0] < -Base::xL / 2.0 ||
                particles[i].pos[2] > Base::zLBox / 2.0 + 2.0 - particles[i].d / 2.0 || particles[i].pos[1] > Base::yL / 2.0 || particles[i].pos[0] > Base::xL / 2.0){
            printf("Error calculating energy, particle %i was found outside the box..\n", i);
            std::cout << particles[i].com << std::endl;
            std::cout << particles[i].pos << std::endl;
            exit(1);
        }
        //if(particles[i]->com != particles[i]->pos){
        //    printf("pos and com are not equal!\n");
        //    exit(1);
        //}
        //energy += phiw(particles[i].pos[2]);
    }
    return energy;
}

double energy::direct::get_central(Particles &particles, Particle &p){
    int k = 0;
    double energy = 0;
    //double lenergy = 0;
    double dist = 0;
    for(int i = 0; i < particles.numOfParticles; i++){
        if(p.index == i) continue;
        dist = p.distance_xy(particles[i]);
        energy += particles[i].q * p.q * 1.0 / dist;
    }

    //energy += phiw(p.pos[2]);
    //}
    return energy;
}   

double energy::direct::get_central_pot(Particles &particles, Particle &p){
    int k = 0;
    double energy = 0;
    //double lenergy = 0;
    double dist = 0;
    //#pragma omp parallel for reduction(+:energy) private(lenergy, dist)
    //{
    for(int i = 0; i < particles.numOfParticles; i++){
        if(p.index == i) continue;

        dist = p.distance_xy(particles[i]);
        if(dist <= 1.0){
            dist = 1.0;
        }
        energy += particles[i].q * p.q * 1.0 / dist;
    }

    //energy += phiw(p.pos[2]);
    //}
    return energy;
}  



double energy::direct::phiw(double z){
    double a = Base::xL/2.0;
    double asq = a * a;
    double zsq = z * z;//(z + b) * (z + b);

    double self = 8.0 * a * std::log((std::sqrt(2.0 * asq + zsq) + a) / std::sqrt(asq + zsq));// - 
                  //2.0 * z * (std::asin((asq * asq - zsq * zsq - 2.0 * asq * zsq) / std::pow(asq + zsq, 2.0)) + PI/2.0);
    return self;
}




double energy::direct::get_replicates(Particles &particles, Particle &p){
    double energy = 0;
    //double lenergy = 0;
    double dist = 0;
    //Eigen::Vector3d disp;
    int rep = 0;
    int mx = rep;
    int my = rep;
    int mz = 0;

    //printf("Calculating energy for %d replicas.\n", (2 * rep+1) * (2 * rep+1) * (2 * rep+1) - 1);
    //#pragma omp parallel for if(rep > 10) reduction(+:energy) private(dist)
    for(int i = -mx; i <= mx; i++){
        for(int j = -my; j <= my; j++){
            for(int k = -mz; k <= mz; k++){
                if(sqrt(i * i + j * j + k * k) <= rep){
                    //numOfRep++;
                    for(int l = 0; l < particles.numOfParticles; l++){
                        if(i == 0 && j == 0 && k == 0){}
                        else{
                            double den[] = {p.pos[0] - particles[l].pos[0] + i * Base::xL, 
                                            p.pos[1] - particles[l].pos[1] + j * Base::yL , 
                                            p.pos[2] - particles[l].pos[2] + k * Base::zLBox};
                            dist = norm(den);
                            //if(dist <= (rep + 1) * Base::xL){//(dist >= Base::xL && dist <= rep * Base::xL){
                            energy += p.q * particles[l].q * 1.0 / dist;
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

double energy::direct::get_replicates_pot(Particles &particles, Particle &p){
    double energy = 0;
    //double lenergy = 0;
    double dist = 0;
    //Eigen::Vector3d disp;
    int rep = 0;
    int mx = rep;
    int my = rep;
    int mz = 0;

    //printf("Calculating energy for %d replicas.\n", (2 * rep+1) * (2 * rep+1) * (2 * rep+1) - 1);
    //#pragma omp parallel for if(rep > 10) reduction(+:energy) private(dist)
    for(int i = -mx; i <= mx; i++){
        for(int j = -my; j <= my; j++){
            for(int k = -mz; k <= mz; k++){
                if(sqrt(i * i + j * j + k * k) <= rep){
                    //numOfRep++;
                    for(int l = 0; l < particles.numOfParticles; l++){
                        if((i == 0 && j == 0 && k == 0) || p.index == l){}
                        else{
                            double den[] = {p.pos[0] - particles[l].pos[0] + i * Base::xL, 
                                            p.pos[1] - particles[l].pos[1] + j * Base::yL , 
                                            p.pos[2] - particles[l].pos[2] + k * Base::zLBox};
                            dist = norm(den);
                            //if(dist <= (rep + 1) * Base::xL){//(dist >= Base::xL && dist <= rep * Base::xL){
                            energy += p.q * particles[l].q * 1.0 / dist;
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
