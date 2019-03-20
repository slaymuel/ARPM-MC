#include "mc.h"
#include <vector>

// double MC::get_energy(Particle **particles){
//     int i = 0;
//     int j = 0;
//     double dist = 0;
//     double energy = 0;
//     double r6 = 0;
//     double sigma = 1;
//     double epsilon = 1;

//     for(i = 0; i < Particle::numOfParticles; i++){
//         j = i + 1;
//         while(j < Particle::numOfParticles){
//             //LJ
//             dist = particles[j]->distance(particles[i]);
//             //r6 = pow(dist, 3);
//             //energy += 4 * (1/(r6*r6) - 1/r6);
//             //Coloumb
//             dist = sqrt(dist);
//             energy += pow(EC, 2)/(4 * VP * PI * 2 * KB * Base::T * 1e-10) * (particles[i]->q * particles[j]->q)/dist;
//             j++;
//         }
//     }
//     return energy;
// }

// double MC::get_particle_energy(int pInd, Particle *p, Particle **particles){
//     int i = 0;
//     double energy = 0;
//     double dist = 0;
//     double r6 = 0;
//     double sigma = 1;
//     double epsilon = 1;

//     for(i = 0; i < Particle::numOfParticles; i++){
//         if(i != pInd){
//             //LJ
//             dist = p->distance_xy(particles[i]);
//             r6 = pow(dist, 3);
//             //energy += 1/kb*(4 * (1/(r6*r6) - 1/r6));
//             //printf("Lennard%lf\n", 1/kb*(4 * (1/(r6*r6) - 1/r6)));
//             //printf("Lennard: %lf\n", energy);
//             //Coloumb
//             dist = sqrt(dist);
//             //printf("Coloumb %lf\n", pow(EC, 2)/(4 * VP * PI * 2 * kb * T  * 1e-10) * (p->q * particles[i]->q)/dist);
//             energy += pow(EC, 2)/(4 * VP * PI * 2 * KB * Base::T * 1e-10) * (p->q * particles[i]->q)/dist;
//             //printf("Coloumb: %lf\n", energy);
//         }
//     }
//     return energy;
// }
/*
int MC::charge_rot_move(Particle **particles){
    double eOld;
    double eNew;
    double acceptProb;
    int r = ran2::get_random() * Particle::numOfParticles;
    double dE;
    int accepted = 0;
    Particle *_old = new Particle(true);

    _old->pos = particles[r]->pos;
    _old->com = particles[r]->com;
    _old->chargeDisp = particles[r]->chargeDisp;
    _old->q = particles[r]->q;
    _old->index = particles[r]->index;
    
    //eOld = MC::direct.get_energy(particles);
    eOld = MC::ewald3D.get_energy(particles);

    particles[r]->random_charge_rot();
    Particle::update_distances(particles, particles[r]);
    MC::ewald3D.update_reciprocal(_old, particles[r]);
    //eNew = MC::direct.get_energy(particles);
    eNew = MC::ewald3D.get_energy(particles);
    dE = eNew - eOld;
    acceptProb = exp(-1 * dE);
    double p = ran2::get_random();

    //Accept move
    if(p <= acceptProb){
        accepted = 1;
        Base::eCummulative += dE;
        Base::acceptedMoves++;
    }
    //Reject
    else{
        MC::ewald3D.update_reciprocal(particles[r], _old);
        particles[r]->chargeDisp = _old->chargeDisp;
        particles[r]->pos = _old->pos;
        Particle::update_distances(particles, particles[r]);
    }
    delete _old;
    return accepted;
}
 */


/*
int MC::trans_move(Particle **particles, double dr){
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
    //eOld = MC::ewald3D.get_energy(particles);
    eOld = MC::ewald3D.get_energy(particles, particles[p]);
    //eOld = energy::valleau::get_energy(particles);
    
    _old->pos = particles[p]->pos;
    _old->com = particles[p]->com;
    _old->q = particles[p]->q;
    _old->index = particles[p]->index;
    
    //Generate new trial coordinates
    particles[p]->random_move(dr);
    Particle::update_distances(particles, particles[p]);

    //If there is no overlap in new position and it's inside the box
    if(particles[p]->hard_sphere(particles) && particles[p]->pos[2] > particles[p]->d/2 + Base::wall &&
                                              particles[p]->pos[2] < Base::zL - Base::wall - particles[p]->d/2 ){

        //Get new energy
        MC::ewald3D.update_reciprocal(_old, particles[p]);
        //eNew = MC::ewald3D.get_energy(particles);
        eNew = MC::ewald3D.get_energy(particles, particles[p]);
        //eNew = energy::valleau::get_energy(particles);

        //printf("new %lf\n", eNew);
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
            //printf("Bad move\n");
            MC::ewald3D.update_reciprocal(particles[p], _old);
            particles[p]->pos = _old->pos;
            particles[p]->com = _old->com;
            Particle::update_distances(particles, particles[p]);
        }
    }

    else{   //Reject move
        //printf("Outside box\n");
        particles[p]->pos = _old->pos;
        particles[p]->com = _old->com;
        Particle::update_distances(particles, particles[p]);
    }

    //free(_old->pos);
    delete _old;
    return accepted;
}
*/

//void MC::equilibrate(Particle **particles){
void MC::equilibrate(){
    int overlaps = 1;
    int prevOverlaps = 500000;
    int i = 0;
    int j = 0;
    double diameter2 = pow(particles[0].d, 2);
    Eigen::Vector3d oldPos;
    Eigen::Vector3d oldCom;
    double random;
    double stepSize = 5;
    int p;

    //Move particles to prevent overlap
    while(overlaps > 0){
        random = ran2::get_random();
        p =  random * particles.numOfParticles;
        oldCom = particles[p].com;
        oldPos = particles[p].pos;
        // oldPos[0] = particles[p]->com[0];
        // oldPos[1] = particles[p]->com[1];
        // oldPos[2] = particles[p]->com[2];
        particles[p].random_move(5);
        //particles.update_distances(particles[p]);
        if(!particles.hard_sphere(particles[p]) || particles[p].com[2] < -Base::zLBox / 2.0 + particles[p].d / 2.0 ||
                                                    particles[p].com[2] > Base::zLBox / 2.0 - particles[p].d / 2.0){
            particles[p].com = oldCom;
            particles[p].pos = oldPos;
            //particles.update_distances(particles[p]);
        }

        else if(particles[p].com[2] < -Base::zLBox / 2.0 + particles[p].d / 2.0 || particles[p].com[1] < -Base::yL / 2.0 || particles[p].com[0] < -Base::xL / 2.0){
            printf("Failed to equilibrate system, a particle was found outside the box, %lf, %lf, %lf...\n", particles[p].com[0], particles[p].com[1], particles[p].com[2]);
            exit(1);
        }

        if(i % 5000 == 0){
            overlaps = particles.get_overlaps();
            //stepSize = log(abs(prevOverlaps - overlaps) + 3.0);
            //prevOverlaps = overlaps;
            printf("Overlaps: %d, iteration: %d\r", overlaps, i);
            fflush(stdout);
        }
        i++;
        if(i > 1E9){
            i = 0;
        }
    }
    printf("\nEquilibration done\n\n");
}
/*
void MC::disperse(Particle **particles) {
    int i = 0;
    int j = 0;
    int k = 0;
    int index = 0;
    int overlaps = Particle::get_overlaps(particles);
    double random;
    int p;
    double *oldPos = (double *) malloc(3 * sizeof(double));
    std::vector<int> indices(Particle::numOfParticles);

    for (i = 0; i < indices.size(); i++) {
        indices[i] = i;
    }

    while (overlaps > 0) {
        for (i = 0; i < indices.size(); i++) {
            if (k % 100 == 0) {
                random = ran2::get_random();
                p = random * Particle::numOfParticles;
                oldPos[0] = particles[p]->pos[0];
                oldPos[1] = particles[p]->pos[1];
                oldPos[2] = particles[p]->pos[2];
                particles[p]->random_move(Base::xL);
                if (!particles[p]->hard_sphere(particles) || particles[p]->pos[2] < particles[p]->d / 2 ||
                    particles[p]->pos[2] > Base::zL - particles[p]->d / 2) {
                    particles[p]->pos[0] = oldPos[0];
                    particles[p]->pos[1] = oldPos[1];
                    particles[p]->pos[2] = oldPos[2];
                }
            } else {
                oldPos[0] = particles[indices[i]]->pos[0];
                oldPos[1] = particles[indices[i]]->pos[1];
                oldPos[2] = particles[indices[i]]->pos[2];
                particles[indices[i]]->random_move(5);
                if (!particles[indices[i]]->hard_sphere(particles) ||
                    particles[indices[i]]->pos[2] < particles[indices[i]]->d / 2 ||
                    particles[indices[i]]->pos[2] > Base::zL - particles[indices[i]]->d / 2) {
                    particles[indices[i]]->pos[0] = oldPos[0];
                    particles[indices[i]]->pos[1] = oldPos[1];
                    particles[indices[i]]->pos[2] = oldPos[2];
                } else {
                    indices.erase(indices.begin() + i);
                }
            }
            k++;
        }

        j++;
        if (j % 10 == 0) {
            overlaps = Particle::get_overlaps(particles);
            printf("Overlapping elements: %lu, Overlaps: %d\r", indices.size(), overlaps);
            fflush(stdout);
        }
    }
}
*/


/*
template<typename F>
void MC::run(F&& energy_function, Particle** particles, int iter){
    double energy;
    int prevAccepted = 0;
    Base::eCummulative = energy_function(particles);
    for(int i = 0; i < iter; i++){
        //trans_move(particles, 0.1, energy_function);
        //energy = energy_function(particles);
        //printf("rvalue reference energy function called: %lf\n", energy);
        Base::totalMoves++;

        if(i % 1 == 0 && i != 0){
            energy = energy_function(particles);
            //Particle::write_coordinates(outName , particles);
            printf("Iteration: %d\n", i);
            printf("Energy: %lf\n", energy);
            printf("Acceptance ratio: %lf\n", (double)Base::acceptedMoves/Base::totalMoves);
            printf("Acceptance ratio for the last 10000 steps: %lf\n\n", (double)prevAccepted/1.0);
            if(fabs(energy - Base::eCummulative)/fabs(energy) > pow(10, -12)){
                printf("Error is too large!\n");
                printf("Error: %lf\n", fabs(energy - Base::eCummulative)/fabs(energy));
                exit(1);
            }
            prevAccepted = 0;
        }
    }
}
*/
    // for(j = 0; j < 10; j++){
    //     printf("Iteration: %d of 10, accepted dispersion moves %d\r", j, k);
    //     fflush(stdout);
    //     for(i = 0; i < Particle::numOfParticles; i++){
    //         oldPos[0] = particles[i]->pos[0];
    //         oldPos[1] = particles[i]->pos[1];
    //         oldPos[2] = particles[i]->pos[2];

    //         if(i % 9 == 0){
    //             particles[i]->random_move(Base::xL);
    //         }
    //         else{
    //             particles[i]->random_move(6);
    //         }

    //         if(!particles[i]->hard_sphere(particles) || particles[i]->pos[2] < particles[i]->d/2 || particles[i]->pos[2] > Base::zL - particles[i]->d/2){
    //             particles[i]->pos[0] = oldPos[0];
    //             particles[i]->pos[1] = oldPos[1];
    //             particles[i]->pos[2] = oldPos[2];
    //         }
    //         else{
    //             k++;
    //         }
    //     }
    // }
    //printf("\n");
//}
