#include "mc.h"

double MC::getEnergy(Particle **particles){
    int i = 0;
    int j = 0;
    double dist = 0;
    double energy = 0;
    double r6 = 0;
    double sigma = 1;
    double epsilon = 1;

    for(i = 0; i < Particle::numOfParticles; i++){
        j = i + 1;
        while(j < Particle::numOfParticles){
            //LJ
            dist = particles[j]->distance_xy(particles[i]);
            //r6 = pow(dist, 3);
            //energy += 4 * (1/(r6*r6) - 1/r6);
            //Coloumb
            dist = sqrt(dist);
            energy += pow(EC, 2)/(4 * VP * PI * 2 * KB * Base::T * 1e-10) * (particles[i]->q * particles[j]->q)/dist;
            j++;
        }
    }
    return energy;
}

double MC::getParticleEnergy(int pInd, Particle *p, Particle **particles){
    int i = 0;
    double energy = 0;
    double dist = 0;
    double r6 = 0;
    double sigma = 1;
    double epsilon = 1;

    for(i = 0; i < Particle::numOfParticles; i++){
        if(i != pInd){
            //LJ
            dist = p->distance_xy(particles[i]);
            r6 = pow(dist, 3);
            //energy += 1/kb*(4 * (1/(r6*r6) - 1/r6));
            //printf("Lennard%lf\n", 1/kb*(4 * (1/(r6*r6) - 1/r6)));
            //printf("Lennard: %lf\n", energy);
            //Coloumb
            dist = sqrt(dist);
            //printf("Coloumb %lf\n", pow(EC, 2)/(4 * VP * PI * 2 * kb * T  * 1e-10) * (p->q * particles[i]->q)/dist);
            energy += pow(EC, 2)/(4 * VP * PI * 2 * KB * Base::T * 1e-10) * (p->q * particles[i]->q)/dist;
            //printf("Coloumb: %lf\n", energy);
        }
    }
    return energy;
}

int MC::mcmove(Particle **particles, double dr){
    int i = 0;
    double eOld = 0;
    double eNew = 0;
    double dist = 0;
    double acceptProp = 0;
    double random = ran2::get_random();
    double dE = 0;
    int accepted= 0;

    int p =  random * Particle::numOfParticles;

    //Calculate old energy
    eOld = MC::getParticleEnergy(p, particles[p], particles);

    //Generate new trial coordinates
    double oldPos[3] = {particles[p]->pos[0], particles[p]->pos[1], particles[p]->pos[2]};
    particles[p]->pos[0] = particles[p]->pos[0] + (ran2::get_random()*2 - 1) * dr;
    particles[p]->pos[1] = particles[p]->pos[1] + (ran2::get_random()*2 - 1) * dr;
    particles[p]->pos[2] = particles[p]->pos[2] + (ran2::get_random()*2 - 1) * dr;

    //Appy PBC
    particles[p]->pbc_xy();

    //If there is no overlap in new position and it's inside the box
    if(particles[p]->hardSphere(particles) && particles[p]->pos[2] > particles[p]->d/2 && particles[p]->pos[2] < Base::zL - particles[p]->d/2 ){
        //Get new energy
        eNew = MC::getParticleEnergy(p, particles[p], particles);

        //Accept move?
        dE = eNew - eOld;
        acceptProp = exp(-1*dE);
        if(acceptProp > 1 || eNew < eOld){
            acceptProp = 1;
        }
        double rand = ran2::get_random();
        if(rand < acceptProp){  //Accept move
            Base::eCummulative += dE; //Cummulative energy
            accepted = 1;
            Base::acceptedMoves++;
        }
        else{   //Reject move
            particles[p]->pos[0] = oldPos[0];
            particles[p]->pos[1] = oldPos[1];
            particles[p]->pos[2] = oldPos[2];
        }
    }

    else{   //Reject move
        particles[p]->pos[0] = oldPos[0];
        particles[p]->pos[1] = oldPos[1];
        particles[p]->pos[2] = oldPos[2];
    }
    return accepted;
}

void MC::equilibrate(Particle **particles){
    int overlaps = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    double diameter2 = pow(particles[0]->d, 2);
    double *oldPos = (double*) malloc(3*sizeof(double));
    float random;
    int p;

    //Move particles to prevent overlap
    while(k < 100000){
        overlaps = 0;
        for(i = 0; i < Particle::numOfParticles; i++){
            j = i + 1;
            while(j < Particle::numOfParticles){
                if(i != j){
                    if(particles[i]->distance_xy(particles[j]) < diameter2){
                        particles[i]->randomMove_xy(1);
                        // if(particles[i]->pos[2] < particles[i]->d/2 || particles[i]->pos[2] > Base::zL - particles[i]->d/2){
                        //     printf("Particle found in forbidden zone!\n");
                        //     exit(1);
                        // }

                        overlaps += 1;
                    }
                }
                j++;
            }
        }
        if(overlaps < 1){
            break;
        }
        printf("Overlaps: %d\n", overlaps);
        k++;
    }
}

void MC::disperse(Particle **particles){
    int i = 0;
    int j = 0;
    int k = 0;
    double *oldPos = (double*)malloc(3*sizeof(double));

    for(j = 0; j < 100; j++){
        printf("Dispersion iteration: %d of 100, acepted dispersion moves %d\r", j, k);
        fflush(stdout);
        for(i = 0; i < Particle::numOfParticles; i++){
            oldPos[0] = particles[i]->pos[0];
            oldPos[1] = particles[i]->pos[1];
            oldPos[2] = particles[i]->pos[2];

            if(i % 100 == 0){
                particles[i]->randomMove_xy(Base::xL);
            }
            else{
                particles[i]->randomMove_xy(5);
            }

            if(!particles[i]->hardSphere(particles) || particles[i]->pos[2] < particles[i]->d/2 || particles[i]->pos[2] > Base::zL - particles[i]->d/2){
                particles[i]->pos[0] = oldPos[0];
                particles[i]->pos[1] = oldPos[1];
                particles[i]->pos[2] = oldPos[2];
            }
            else{
                k++;
            }
        }
    }
    printf("\n");
}