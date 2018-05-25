#include "mc.h"
#include <vector>

double MC::get_energy(Particle **particles){
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
            dist = particles[j]->distance(particles[i]);
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

double MC::get_particle_energy(int pInd, Particle *p, Particle **particles){
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
    double eOld = 0;
    double eNew = 0;
    double dist = 0;
    double acceptProp = 0;

    double random = ran2::get_random();

    double dE = 0;
    int accepted= 0;
    double ewald3DEnergy = 0;
    double ewald2DEnergy = 0;
    double directEnergy = 0;
    double pbc = 0;
    Particle *_old = new Particle(true);

    int p =  random * Particle::numOfParticles;

    //ewald3DEnergy = MC::ewald3D.get_energy(particles);
    //printf("Ewald3D: %lf\n", ewald3DEnergy);
    //directEnergy = MC::direct.get_energy(particles);
    //printf("Direct: %lf\n", directEnergy);
    //ewald2DEnergy = MC::ewald2D.get_energy(particles);
    //printf("Ewald2D: %lf\n", ewald2DEnergy);

    //exit(1);

    //Calculate old energy
    //eOld = MC::get_particle_energy(p, particles[p], particles);
    eOld = MC::ewald3D.get_energy(particles);
    //Save old particle state
    
    _old->pos = (double*)malloc(3 * sizeof(double));
    _old->pos[0] = particles[p]->pos[0];
    _old->pos[1] = particles[p]->pos[1];
    _old->pos[2] = particles[p]->pos[2];
    _old->q = particles[p]->q;
    _old->index = particles[p]->index;
    
    //printf("oldpos: %lf %lf %lf\n", particles[p]->pos[0], particles[p]->pos[1], particles[p]->pos[2]);
    //Generate new trial coordinates
    particles[p]->pos[0] = particles[p]->pos[0] + (ran2::get_random()*2 - 1) * dr;
    particles[p]->pos[1] = particles[p]->pos[1] + (ran2::get_random()*2 - 1) * dr;
    particles[p]->pos[2] = particles[p]->pos[2] + (ran2::get_random()*2 - 1) * dr;

    //printf("newpos: %lf %lf %lf\n", particles[p]->pos[0], particles[p]->pos[1], particles[p]->pos[2]);
    //Appy PBC
    particles[p]->pbc();
    Particle::update_distances(particles[p], particles);
    //If there is no overlap in new position and it's inside the box
    if(particles[p]->hardSphere(particles) && particles[p]->pos[2] > particles[p]->d/2 + Base::wall && 
                                              particles[p]->pos[2] < Base::zL - Base::wall - particles[p]->d/2 ){
        //Get new energy
        //eNew = MC::get_particle_energy(p, particles[p], particles);
        MC::ewald3D.update_reciprocal(_old, particles[p]);
        eNew = MC::ewald3D.get_energy(particles);

        //Accept move?
        dE = eNew - eOld;
        acceptProp = exp(-1*dE);
        if(acceptProp > 1 || eNew < eOld){
            acceptProp = 1;
        }

        double rand = ran2::get_random();

        if(rand <= acceptProp){  //Accept move
            Base::eCummulative += dE; //Update cummulative energy
            accepted = 1;
            Base::acceptedMoves++;
            //printf("Accept\n");
        }
        else{   //Reject move
            MC::ewald3D.update_reciprocal(particles[p], _old);
            particles[p]->pos[0] = _old->pos[0];
            particles[p]->pos[1] = _old->pos[1];
            particles[p]->pos[2] = _old->pos[2];
            Particle::update_distances(particles[p], particles);
            //if(fabs(Base::lB * MC::ewald3D.get_energy(particles) - eOld) > 1e-4){
            //    printf("oldpos: %lf %lf %lf\n", particles[p]->pos[0], particles[p]->pos[1], particles[p]->pos[2]);
            //    printf("New energy: %lf\n", Base::lB * MC::ewald3D.get_energy(particles));
            //    printf("Old energy: %lf\n", eOld);
            //    printf("Energy failed..\n");
            //    exit(1);
            //}
        }
    }

    else{   //Reject move
        particles[p]->pos[0] = _old->pos[0];
        particles[p]->pos[1] = _old->pos[1];
        particles[p]->pos[2] = _old->pos[2];
        Particle::update_distances(particles[p], particles);
    }

    free(_old->pos);
    delete _old;
    return accepted;
}

void MC::equilibrate(Particle **particles){
    int overlaps = 1;
    int prevOverlaps = 3000;
    int i = 0;
    int j = 0;
    double diameter2 = pow(particles[0]->d, 2);
    double *oldPos = (double*) malloc(3*sizeof(double));
    double random;
    double stepSize = 5;
    int p;

    //Move particles to prevent overlap
    while(overlaps > 0){
        random = ran2::get_random();
        p =  random * Particle::numOfParticles;
        oldPos[0] = particles[p]->pos[0];
        oldPos[1] = particles[p]->pos[1];
        oldPos[2] = particles[p]->pos[2];
        particles[p]->randomMove(5);
        Particle::update_distances(particles[p], particles);
        if(!particles[p]->hardSphere(particles) || particles[p]->pos[2] < particles[p]->d/2 + Base::wall || particles[p]->pos[2] > Base::zL - Base::wall - particles[p]->d/2){
            particles[p]->pos[0] = oldPos[0];
            particles[p]->pos[1] = oldPos[1];
            particles[p]->pos[2] = oldPos[2];
            Particle::update_distances(particles[p], particles);
        }

        if(i % 50000 == 0){
            overlaps = Particle::get_overlaps(particles);
            //stepSize = log(abs(prevOverlaps - overlaps) + 3.0);
            //prevOverlaps = overlaps;
            printf("Overlaps: %d\n", overlaps);
            //fflush(stdout);
        }
        i++;
    }
    printf("\n");
}

void MC::disperse(Particle **particles){
    int i = 0;
    int j = 0;
    int k = 0;
    int index = 0;
    int overlaps = Particle::get_overlaps(particles);
    double random;
    int p;
    double *oldPos = (double*)malloc(3*sizeof(double));
    std::vector<int> indices(Particle::numOfParticles);

    for(i = 0; i < indices.size(); i++){
        indices[i] = i;
    }

    while(overlaps > 0){
        for(i = 0; i < indices.size(); i++){
            if(k % 100 == 0){
                random = ran2::get_random();
                p =  random * Particle::numOfParticles;
                oldPos[0] = particles[p]->pos[0];
                oldPos[1] = particles[p]->pos[1];
                oldPos[2] = particles[p]->pos[2];
                particles[p]->randomMove(Base::xL);
                if(!particles[p]->hardSphere(particles) || particles[p]->pos[2] < particles[p]->d/2 || 
                                                                        particles[p]->pos[2] > Base::zL - particles[p]->d/2){
                        particles[p]->pos[0] = oldPos[0];
                        particles[p]->pos[1] = oldPos[1];
                        particles[p]->pos[2] = oldPos[2];
                    }
            }
            else{
                oldPos[0] = particles[indices[i]]->pos[0];
                oldPos[1] = particles[indices[i]]->pos[1];
                oldPos[2] = particles[indices[i]]->pos[2];
                particles[indices[i]]->randomMove(5);
                if(!particles[indices[i]]->hardSphere(particles) || particles[indices[i]]->pos[2] < particles[indices[i]]->d/2 || 
                                                                        particles[indices[i]]->pos[2] > Base::zL - particles[indices[i]]->d/2){
                        particles[indices[i]]->pos[0] = oldPos[0];
                        particles[indices[i]]->pos[1] = oldPos[1];
                        particles[indices[i]]->pos[2] = oldPos[2];
                    }
                else{
                        indices.erase(indices.begin() + i);
                    }
            }
            k++;
        }

        j++;
        if(j % 10 == 0){
            overlaps = Particle::get_overlaps(particles);
            printf("Overlapping elements: %lu, Overlaps: %d\r", indices.size(), overlaps);
            fflush(stdout);
        }
    }


    // for(j = 0; j < 10; j++){
    //     printf("Iteration: %d of 10, accepted dispersion moves %d\r", j, k);
    //     fflush(stdout);
    //     for(i = 0; i < Particle::numOfParticles; i++){
    //         oldPos[0] = particles[i]->pos[0];
    //         oldPos[1] = particles[i]->pos[1];
    //         oldPos[2] = particles[i]->pos[2];

    //         if(i % 9 == 0){
    //             particles[i]->randomMove(Base::xL);
    //         }
    //         else{
    //             particles[i]->randomMove(6);
    //         }

    //         if(!particles[i]->hardSphere(particles) || particles[i]->pos[2] < particles[i]->d/2 || particles[i]->pos[2] > Base::zL - particles[i]->d/2){
    //             particles[i]->pos[0] = oldPos[0];
    //             particles[i]->pos[1] = oldPos[1];
    //             particles[i]->pos[2] = oldPos[2];
    //         }
    //         else{
    //             k++;
    //         }
    //     }
    // }
    printf("\n");
}