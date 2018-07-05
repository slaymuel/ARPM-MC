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
        //static Ewald3D ewald3D;
        static Ewald2D ewald2D;
        //static Direct direct;




        template <typename E>
        static int vol_move(Particle **particles, E energy_function){
            bool overlap = false;
            double vMax = 0.01;
            double lnNewVolume = std::log(Base::volume) + (ran2::get_random() - 0.5) * vMax;
            double newVolume = std::exp(lnNewVolume);
            double newLength = cbrt(newVolume);
            double volFrac = newLength / Base::xL;
            double oldxL = Base::xL;
            double oldyL = Base::yL;
            double oldzL = Base::zL;

            for(int i = 0; i < Particle::numOfParticles; i++){
                for(int j = i + 1; j < Particle::numOfParticles; j++){
                    if(particles[i]->com_distance(particles[j]) * volFrac < particles[i]->d/2 + particles[j]->d/2){
                        overlap = true;
                        break;
                    }
                }
                if(overlap){
                    break;
                }
            }
            if(!overlap){
                double oldEnergy = energy_function(particles);
                for(int i = 0; i < Particle::numOfParticles; i++){
                    particles[i]->com *= volFrac;
                    particles[i]->pos *= volFrac;
                }
                Base::xL = newLength;
                Base::yL = newLength;
                Base::zL = newLength;
                Particle::update_distances(particles);
                energy::ewald3D::set_alpha();
                energy::ewald3D::reset();
                energy::ewald3D::initialize(particles);
                double newEnergy = energy_function(particles);
                //(0.5 * Base::beta)
                double prob = exp(-(newEnergy - oldEnergy) - Base::beta * 100000 * 1e-30 * (newVolume - Base::volume) - 
                                                (Particle::numOfParticles + 1) * std::log(newVolume / Base::volume));

                if(ran2::get_random() > prob){  //Reject
                    
                    Base::xL = oldxL;
                    Base::yL = oldyL;
                    Base::zL = oldzL;
                    for(int i = 0; i < Particle::numOfParticles; i++){
                        particles[i]->com *= Base::xL / newLength;
                        particles[i]->pos *= Base::xL / newLength;
                    }
                    Particle::update_distances(particles);
                    energy::ewald3D::set_alpha();
                    energy::ewald3D::reset();
                    energy::ewald3D::initialize(particles);
                }
                else{                           //Accept
                    Base::eCummulative += newEnergy - oldEnergy;
                    Base::volume = newVolume;
                    Base::acceptedMoves++;
                    //printf("volume: %lf\n", Base::volume);
                    return 1;
                }
                return 0;
            }
            else{
                return 0;
            }
            /*
            double dVmax = 1.0;
            double accepted = 0;
            double dV = (ran2::get_random()-0.5)*dVmax;
            double newVolume = Base::volume + dV;
            double oldVolume = Base::volume;
            double newLength = cbrt(newVolume); //ted2
            double oldLength = Base::xL;
            double dLength = newLength / Base::xL;
            //double oldEnergy = energy_function(particles);
       
            for(int i = 0; i < Particle::numOfParticles; i++){
                particles[i]->com *= dLength;
                particles[i]->pos *= dLength;
            }
            Base::xL = newLength;
            Base::yL = newLength;
            Base::zL = newLength;

            //Particle::update_distances(particles);
            //double newEnergy = energy_function(particles);

            double dut = Base::beta * Base::P * 0.000000000000000000000000000001 * dV - Particle::numOfParticles * std::log(newVolume / Base::volume);

            if (std::exp(-dut) >= ran2::get_random()){  //accept
                //Base::eCummulative += newEnergy - oldEnergy;
                Base::volume = newVolume;
                accepted = 1;
                //printf("volume: %lf\n", Base::volume);
            }
            else{                           //Reject
                for(int i = 0; i < Particle::numOfParticles; i++){
                    particles[i]->com *= 1/dLength;
                    particles[i]->pos *= 1/dLength;
                    Base::xL = oldLength;
                    Base::yL = oldLength;
                    Base::zL = oldLength;
                    Base::volume = oldVolume;
                    //Particle::update_distances(particles);
                }                
            }
            return accepted;
            */
        }



        template <typename E>
        static int charge_rot_move(Particle **particles, E energy_function){
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
            eOld = energy_function(particles, particles[r]);

            particles[r]->random_charge_rot();
            Particle::update_distances(particles, particles[r]);
            energy::ewald3D::update_reciprocal(_old, particles[r]);
            //eNew = MC::direct.get_energy(particles);
            eNew = energy_function(particles, particles[r]);
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
                energy::ewald3D::update_reciprocal(particles[r], _old);
                particles[r]->chargeDisp = _old->chargeDisp;
                particles[r]->pos = _old->pos;
                Particle::update_distances(particles, particles[r]);
            }
            delete _old;
            return accepted;
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
            if(particles[p]->hard_sphere(particles)){// && particles[p]->com[2] > particles[p]->d/2 + Base::wall &&
                //particles[p]->com[2] < Base::zL - Base::wall - particles[p]->d/2 ){
                //Get new energy
                energy::ewald3D::update_reciprocal(_old, particles[p]);
                //energy::levin::update_f(_old, particles[p]);
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
                    energy::ewald3D::update_reciprocal(particles[p], _old);
                    //energy::levin::update_f(particles[p], _old);
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
            int rotAccepted = 0;
            int volAccepted = 0;
            int transAccepted = 0;
            int outFreq = 1000;
            double random = 0;

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

                random = ran2::get_random();
                if(random <= 0.7){
                    if(trans_move(particles, dr, particle_energy_function)){
                        prevAccepted++; 
                        transAccepted++;
                    }
                    Base::totalMoves++;
                }

                random = ran2::get_random();
                if(random <= 0.3){
                    if(charge_rot_move(particles, particle_energy_function)){
                        prevAccepted++;
                        rotAccepted++;
                    }
                    Base::totalMoves++;
                }

                random = ran2::get_random();
                if(random <= 0.1){
                    if(vol_move(particles, energy_function)){
                        prevAccepted++;
                        volAccepted++;
                    }
                    Base::totalMoves++;
                }

                
                if(i % 100 && i > 10000 && !sample){
                    //energy::valleau::update_charge_vector(particles);
                }

                if(i % outFreq == 0){
                    Base::volumes.push_back(Base::volume);
                    energy_temp = energy_function(particles);
                    //Particle::write_coordinates(outName , particles);
                    printf("Iteration: %d\n", i);
                    printf("Volume: %lf\n", Base::volume);
                    printf("Energy: %lf\n", energy_temp);
                    printf("Acceptance ratio: %lf\n", (double)Base::acceptedMoves/Base::totalMoves);
                    printf("Acceptance ratio for the last %i steps: %lf\n", outFreq, (double)prevAccepted/outFreq);
                    if(std::abs(energy_temp - Base::eCummulative)/std::abs(energy_temp) > std::pow(10, -12)){
                        printf("Error is too large!\n");
                        printf("Error: %.12lf\n", std::abs(energy_temp - Base::eCummulative)/std::abs(energy_temp));
                        exit(1);
                    }
                    printf("Trans moves: %d         Rot moves: %d       Vol moves: %d\n\n", transAccepted, rotAccepted, volAccepted);
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