#ifndef MC_H
#define MC_H

#include "base.h"
#include "particle.h"
#include "ran2_lib.cpp"
#include "ewald3D.h"
#include "ewald2D.h"
#include "analysis.h"
#include "valleau.h"
//#include "direct.h"
#include "levin.h"
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

        template <typename E>
        static int vol_move(Particle **particles, E energy_function){
            double vMax = 2.0;
            double oldEnergy = energy_function(particles);
            double lnNewLength = std::log(Base::xL) + (ran2::get_random() - 0.5) * vMax;
            double newVolume = std::exp(lnNewLength);
            double newLength = cbrt(newVolume);

            for(int i = 0; i < Particle::numOfParticles; i++){
                particles[i]->pos *= newLength / Base::xL;
            }
            double newEnergy = energy_function(particles);

            double prob = exp(-Base::beta * ((newEnergy - oldEnergy) + Base::P * (newVolume - Base::volume) - 
                                             (Particle::numOfParticles + 1) * std::log(newVolume / Base::volume)/Base::beta));

            if(ran2::get_random() > prob){  //Reject
                for(int i = 0; i < Particle::numOfParticles; i++){
                    particles[i]->pos *= Base::xL / newLength;
                }
            }
            else{
                Base::xL = newLength;
                Base::yL = newLength;
                Base::zL = newLength;
                Base::volume = newVolume;
                return 1;
            }
            return 0;
        }

        template<typename E>
        static int trans_move(Particle **particles, double dr, E energy_function){
            double eOld = 0.0;
            double eNew = 0.0;
            double dist = 0.0;
            double acceptProp = 0.0;

            double random = ran2::get_random();

            double dE = 0.0;
            int accepted= 0.0;
            double ewald3DEnergy = 0.0;
            double ewald2DEnergy = 0.0;
            double directEnergy = 0.0;
            Particle *_old = new Particle(true);

            int p =  random * Particle::numOfParticles;

            //Calculate old energy
            eOld = energy_function(particles, particles[p]);
            _old->pos = particles[p]->pos;
            _old->com = particles[p]->com;
            _old->q = particles[p]->q;
            _old->index = particles[p]->index;

            //Generate new trial coordinates
            particles[p]->random_move(dr);
            Particle::update_distances(particles, particles[p]);
            //If there is no overlap in new position and it's inside the box
            if(particles[p]->hard_sphere(particles) && particles[p]->com[2] > particles[p]->d/2 + Base::wall &&
                particles[p]->com[2] < Base::zL - Base::wall - particles[p]->d/2 ){
                //Get new energy
                //MC::ewald3D.update_reciprocal(_old, particles[p]);
                energy::levin::update_f(_old, particles[p]);
                eNew = energy_function(particles, particles[p]);

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
                    energy::levin::update_f(particles[p], _old);
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

        template<typename F, typename FP>
        static void run(F&& energy_function, FP&& particle_energy_function, Particle** particles, double dr, int iter, bool sample){
            double energy_temp;
            int prevAccepted = 0;
            int outFreq = 10000;

            Analysis *xHist = new Analysis(0.1, Base::xL);
            Analysis *yHist = new Analysis(0.1, Base::yL);
            Analysis *zHist = new Analysis(0.1, Base::zL);

            char outName[] = ".txt";
            Base::eCummulative = energy_function(particles);
            Particle::oldEnergy = Base::eCummulative;
            printf("\nRunning MC-loop at temperature: %lf, Bjerrum length is %lf\n\n", Base::T, Base::lB);
            for(int i = 0; i < iter; i++){

                if(i % 100 == 0 && i >= 500000 && sample){
                    //rdf->sample_rdf(particles, histo, binWidth);
                    xHist->sampleHisto(particles, 0);
                    yHist->sampleHisto(particles, 1);
                    zHist->sampleHisto(particles, 2);
                }

                if(trans_move(particles, dr, particle_energy_function)){
                    prevAccepted++; 
                }    

                Base::totalMoves++;
                if(i % 100 && i > 10000 && !sample){
                    energy::valleau::update_charge_vector(particles);
                }

                if(i % outFreq == 0){
                    energy_temp = energy_function(particles);
                    //Particle::write_coordinates(outName , particles);
                    printf("Iteration: %d\n", i);
                    printf("Energy: %lf\n", energy_temp);
                    printf("Acceptance ratio: %lf\n", (double)Base::acceptedMoves/Base::totalMoves);
                    printf("Acceptance ratio for the last %i steps: %lf\n\n", outFreq, (double)prevAccepted/outFreq);
                    if(std::abs(energy_temp - Base::eCummulative)/std::abs(energy_temp) > std::pow(10, -12)){
                        printf("Error is too large!\n");
                        printf("Error: %.12lf\n", std::abs(energy_temp - Base::eCummulative)/std::abs(energy_temp));
                        exit(1);
                    }
                    prevAccepted = 0;
                }
            }
            if(sample){
                xHist->saveHisto(outName);
                yHist->saveHisto(outName);
                zHist->saveHisto(outName);
            }
            delete xHist;
            delete yHist;
            delete zHist;
        }
};

#endif