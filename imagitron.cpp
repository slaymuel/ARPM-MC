#include "imagitron.h"

double energy::imagitron::get_particle_energy(Particle **particles, Particle *p){
    int numOfReflections = 1;
    double energy = 0;
    Eigen::Vector3d disp;
    double distance = 0;

    for(int i = 0; i < p->index; i++){
        energy += particles[i]->q * p->q * 1.0/Particle::distances[i][p->index];
    }

    for(int i = p->index + 1; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
        energy += particles[i]->q * p->q * 1.0/Particle::distances[p->index][i];
    }
    if(p->index < Particle::numOfParticles){
        energy += 2 * PI * 1.0/50.0 * std::fabs(p->pos[2]) * p->q;
    }

    for(int k = 1; k <= numOfReflections; k += 2){

        //Reflections created by right hand side electrons
        for(int i = Particle::numOfParticles; i < Particle::numOfParticles + Particle::numOfElectrons / 2; i++){
            //Left images, uneven, opposite sign

            disp << p->pos[0] - particles[i]->pos[0],
                    p->pos[1] - particles[i]->pos[1],
                    p->pos[2] + particles[i]->pos[2];//p->pos[2] - ((1 - k) * Base::zL - particles[i]->pos[2]);
            //disp = p->pos - particles[i]->pos;
            if(p->index == i){
                energy -= p->q * particles[i]->q * 1.0 / (disp.norm() * 2);
            }
            else{
                energy -= p->q * particles[i]->q * 1.0 / disp.norm();
            }

            
            //Right images, even, same sign
            /*disp << p->pos[0] - particles[i]->pos[0],
                    p->pos[1] - particles[i]->pos[1],
                    p->pos[2] - (k * Base::zL + particles[i]->pos[2]);
            if(p->index == i){
                energy += p->q * particles[i]->q * 1.0 / (disp.norm() * 2);
            }
            else{
                energy += p->q * particles[i]->q * 1.0 / disp.norm();
            }*/
        }
/*
        //Reflections created by left hand side electrons    
        for(int i = Particle::numOfParticles + Particle::numOfElectrons / 2; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
            //Left images, even, same sign
            disp << p->pos[0] - particles[i]->pos[0],
                    p->pos[1] - particles[i]->pos[1],
                    p->pos[2] - (-k * Base::zL + particles[i]->pos[2]);
            if(p->index == i){
                energy += p->q * particles[i]->q * 1.0 / (disp.norm() * 2);
            }
            else{
                energy += p->q * particles[i]->q * 1.0 / disp.norm();
            }

            //Right images, uneven, opposite sign
            disp << p->pos[0] - particles[i]->pos[0],
                    p->pos[1] - particles[i]->pos[1],
                    p->pos[2] - (k * Base::zL - particles[i]->pos[2]);
            if(p->index == i){
                energy -= p->q * particles[i]->q * 1.0 / (disp.norm() * 2);
            }
            else{
                energy -= p->q * particles[i]->q * 1.0 / disp.norm();
            }
        }*/
    }
   
    return energy * Base::lB;
}

double energy::imagitron::get_energy(Particle **particles){
    int numOfReflections = 1;
    double energy = 0;
    Eigen::Vector3d disp;
    double distance = 0;

    for(int i = 0; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
        if(i < Particle::numOfParticles){
            energy += 2 * PI * 1.0/50.0 * std::fabs(particles[i]->pos[2]) * particles[i]->q;
        }
        for(int j = i + 1; j < Particle::numOfParticles + Particle::numOfElectrons; j++){
            energy += particles[i]->q * particles[j]->q * 1.0/Particle::distances[i][j];
        }
    }

    for(int k = 1; k <= numOfReflections; k += 2){
        for(int j = 0; j < Particle::numOfParticles + Particle::numOfElectrons; j++){

            //Reflections created by right hand side electrons
            for(int i = Particle::numOfParticles; i < Particle::numOfParticles + Particle::numOfElectrons / 2; i++){
                //Left images, uneven, opposite sign
                disp << particles[j]->pos[0] - particles[i]->pos[0],
                        particles[j]->pos[1] - particles[i]->pos[1],
                        particles[j]->pos[2] + particles[i]->pos[2];//particles[j]->pos[2] - ((1 - k) * Base::zL - particles[i]->pos[2]);
                //printf("%lf\n", particles[i]->pos[2]);
                //disp = particles[j]->pos - particles[i]->pos;
                
                if(j == i){
                    energy -= particles[j]->q * particles[i]->q * 1.0 / (disp.norm() * 2);
                }
                else{
                    energy -= particles[j]->q * particles[i]->q * 1.0 / disp.norm();
                }

                /*
                //Right images, even, same sign
                disp << particles[j]->pos[0] - particles[i]->pos[0],
                        particles[j]->pos[1] - particles[i]->pos[1],
                        particles[j]->pos[2] - (k * Base::zL + particles[i]->pos[2]);
                if(j == i){
                    energy += particles[j]->q * particles[i]->q * 1.0 / (disp.norm() * 2);
                }
                else{
                    energy += particles[j]->q * particles[i]->q * 1.0 / disp.norm();
                }*/

            }
/*
            //Reflections created by left hand side electrons    
            for(int i = Particle::numOfParticles + Particle::numOfElectrons / 2; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
                //Left images, even, same sign
                disp << particles[j]->pos[0] - particles[i]->pos[0],
                        particles[j]->pos[1] - particles[i]->pos[1],
                        particles[j]->pos[2] - (-k * Base::zL + particles[i]->pos[2]);
                if(j == i){
                    energy += particles[j]->q * particles[i]->q * 1.0 / (disp.norm() * 2);
                }
                else{
                    energy += particles[j]->q * particles[i]->q * 1.0 / disp.norm();
                }

                //Right images, uneven, opposite sign
                disp << particles[j]->pos[0] - particles[i]->pos[0],
                        particles[j]->pos[1] - particles[i]->pos[1],
                        particles[j]->pos[2] - (k * Base::zL - particles[i]->pos[2]);
                if(j == i){
                    energy -= particles[j]->q * particles[i]->q * 1.0 / (disp.norm() * 2);
                }
                else{
                    energy -= particles[j]->q * particles[i]->q * 1.0 / disp.norm();
                }

            }*/

        }
    }

    return energy * Base::lB;
}