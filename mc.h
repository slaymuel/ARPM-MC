#ifndef MC_H
#define MC_H

#include "base.cpp"
#include "particle.h"
#include "ran2_lib.cpp"
#include "ewald3D.h"
#include "ewald2D.h"
//#include "valleau.h"
//#include "direct.h"
//#include "levin.cpp"
class MC{
    public:
        void equilibrate(Particle **particles);
        //static int trans_move(Particle **particles, double dr);
        static int charge_rot_move(Particle **particles);
        static double get_particle_energy(int pInd, Particle *p, Particle **particles);
        double get_energy(Particle **particles);
        void disperse(Particle **particles);
        /*template<typename F>
        static void run(F&& energy_function, Particle** particles, int iter);*/
        static Ewald3D ewald3D;
        static Ewald2D ewald2D;
        //static Direct direct;

        template<typename E>
        static int trans_move(Particle **particles, double dr, E energy_function){
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

            //Calculate old energy
            eOld = energy_function(particles);
            _old->pos = particles[p]->pos;
            _old->com = particles[p]->com;
            _old->q = particles[p]->q;
            _old->index = particles[p]->index;

            //Generate new trial coordinates
            particles[p]->random_move(dr);
            Particle::update_distances(particles, particles[p]);
            //If there is no overlap in new position and it's inside the box
            if(particles[p]->hard_sphere(particles)){// && particles[p]->pos[2] > particles[p]->d/2 + Base::wall &&
                //particles[p]->pos[2] < Base::zL - Base::wall - particles[p]->d/2 ){
                //Get new energy
                //MC::ewald3D.update_reciprocal(_old, particles[p]);
                eNew = energy_function(particles);

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
                }
                else{   //Reject move
                    //MC::ewald3D.update_reciprocal(particles[p], _old);
                    particles[p]->pos = _old->pos;
                    particles[p]->com = _old->com;
                    Particle::update_distances(particles, particles[p]);
                }
            }

            else{   //Reject move
                particles[p]->pos = _old->pos;
                particles[p]->com = _old->com;
                Particle::update_distances(particles, particles[p]);
            }

            delete _old;
            return accepted;
        }

        template<typename F>
        static void run(F&& energy_function, Particle** particles, double dr, int iter){
            double energy;
            int prevAccepted = 0;
            int outFreq = 1000;

            Base::eCummulative = energy_function(particles);
            for(int i = 0; i < iter; i++){
                if(trans_move(particles, dr, energy_function)){
                    prevAccepted++; 
                }    

                Base::totalMoves++;

                if(i % outFreq == 0 && i != 0){
                    energy = energy_function(particles);
                    //Particle::write_coordinates(outName , particles);
                    printf("Iteration: %d\n", i);
                    printf("Energy: %lf\n", energy);
                    printf("Acceptance ratio: %lf\n", (double)Base::acceptedMoves/Base::totalMoves);
                    printf("Acceptance ratio for the last 10000 steps: %lf\n\n", (double)prevAccepted/outFreq);
                    if(fabs(energy - Base::eCummulative)/fabs(energy) > pow(10, -12)){
                        printf("Error is too large!\n");
                        printf("Error: %lf\n", fabs(energy - Base::eCummulative)/fabs(energy));
                        exit(1);
                    }
                    prevAccepted = 0;
                }
            }
        }

};

#endif