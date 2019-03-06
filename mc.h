#ifndef MC_H
#define MC_H

#include "base.h"
#include "particle.h"
#include "particles.h"
#include "ran2_lib.cpp"
#include "ewald3D.h"
#include "ewald2D.h"
#include "analysis.h"
#include "valleau.h"
//#include "img_rep.h"
//#include "direct.h"
#include "levin.h"

class MC{
    public:
        Particles particles;

        void equilibrate();

        template <typename P>
        int grand_move(P particle_energy_function){
            int accept = 0, r;
            double newEnergy, prob;
            double chemPot = -16.0;//-10.7;
            double Donnan = 1200.0;//-22.5;
            double volumeP = (Base::zLBox - 1.0) * Base::xL * Base::yL;
            double volumeN = (Base::zLBox - 5.0) * Base::xL * Base::yL;
            
            if(ran2::get_random() <= 0.5){
                if(ran2::get_random() <= 0.5){  //Delete cation
                    r = ran2::get_random() * particles.numOfCations;
                    r = particles.cations[r];
                    if(r >= particles.numOfParticles){
                        printf("out of bounds cations\n");
                        exit(1);
                    }
                    newEnergy = particle_energy_function(particles, particles[r]);
                    prob = particles.numOfCations / volumeP * std::exp(Donnan * particles[r].q - chemPot + newEnergy);
                }

                else{   //Delete anion
                    r = ran2::get_random() * particles.numOfAnions;
                    r = particles.anions[r];
                    if(r >= particles.numOfParticles){
                        printf("out of bounds anion %i %i %lu\n", r, particles.numOfParticles, particles.anions.size());
                        exit(1);
                    }
                    newEnergy = particle_energy_function(particles, particles[r]);
                    prob =  particles.numOfAnions / volumeN * std::exp(Donnan * particles[r].q - chemPot + newEnergy);
                }
                
                if(ran2::get_random() < prob){  //Accept deletion
                    particles.remove(r);
                    Base::eCummulative -= newEnergy;
                    Base::acceptedMoves++;
                    accept = 1;
                }

                else{
                    accept = 0;
                }
            }


            //Add particle
            else{
                double q;
                if(ran2::get_random() <= 0.5){
                    q = -1.0;
                }
                else{
                    q = 1.0;
                }
                if(particles.add(q)){
                    newEnergy = particle_energy_function(particles, particles[particles.numOfParticles - 1]);
                    if(q > 0.0){    //Add cation
                        prob = volumeP / particles.numOfCations * std::exp(chemPot - Donnan * q - newEnergy);
                    }
                    else{   //add anion
                        prob =  volumeN / particles.numOfAnions * std::exp(chemPot - Donnan * q - newEnergy);
                    }
                    
                    if( ran2::get_random() < prob){ // Accept addition
                        Base::eCummulative += newEnergy;
                        Base::acceptedMoves++;
                        accept = 1;
                    }

                    else{ //Reject
                        particles.remove(particles.numOfParticles - 1);
                        accept = 0;
                    }
                }
                else{
                    accept = 0;
                }
            }
            return accept;
        }







        template <typename E>
        int vol_move(E energy_function){
            bool overlap = false;
            double vMax = 0.00025;
            double lnNewVolume = std::log(Base::volume) + (ran2::get_random() - 0.5) * vMax;

            //double newVolume = Base::volume + (ran2::get_random() - 0.5) * vMax;
            double newVolume = std::exp(lnNewVolume);
            //printf("Old volume %lf\n", Base::volume);
            //printf("New volume %lf\n", newVolume);
            double newLength = cbrt(newVolume);
            double lengthFrac = newLength / Base::xL;
            double oldLength = Base::xL;

            for(int i = 0; i < particles.numOfParticles; i++){
                for(int j = i + 1; j < particles.numOfParticles; j++){
                    if(particles[i].com_distance(particles[j]) * lengthFrac < particles[i].d/2 + particles[j].d/2){
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

                Base::xL = newLength;
                Base::yL = newLength;
                Base::zL = newLength;
                for(int i = 0; i < particles.numOfParticles; i++){
                    particles[i].oldCom = particles[i].com;
                    particles[i].oldPos = particles[i].pos;

                    particles[i].com *= newLength / oldLength;
                    particles[i].pos = particles[i].com + particles[i].chargeDisp;
                    particles[i].pbc_pos();
                }

                particles.update_distances();
                //energy::ewald3D::set_alpha();
                energy::ewald3D::reset();
                energy::ewald3D::initialize(particles);
                double newEnergy = energy_function(particles);
                //(0.5 * Base::beta)
                //printf("dU: %lf         dV: %lf         lnV: %lf\n", newEnergy - oldEnergy, 
                //                                                    Base::beta * 100000 * 1e-30 * (newVolume - Base::volume), 
                //                                                    (particles.numOfParticles + 1) * std::log(newVolume / Base::volume));
                //double prob = exp(-(newEnergy - oldEnergy) - Base::beta * (100000 * 1e-30 * (newVolume - Base::volume) - 
                //                                (particles.numOfParticles + 1) * std::log(newVolume / Base::volume)/Base::beta));
                double prob = exp(-(newEnergy - oldEnergy) - 0.000024305278638 * (newVolume - Base::volume) + //  0.0000243     0.005      0.00383374 0.000024305278638
                                                (particles.numOfParticles + 1) * std::log(newVolume / Base::volume));
                if(ran2::get_random() > prob && oldEnergy <= newEnergy){  //Reject
                    
                    Base::xL = oldLength;
                    Base::yL = oldLength;
                    Base::zL = oldLength;
                    for(int i = 0; i < particles.numOfParticles; i++){
                        particles[i].com = particles[i].oldCom;
                        particles[i].pos = particles[i].oldPos;
                    }
                    particles.update_distances();
                    //energy::ewald3D::set_alpha();
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
        }









        template <typename E>
        int charge_rot_move(E energy_function){
            double eOld;
            double eNew;
            double acceptProb;
            int r = ran2::get_random() * particles.numOfParticles;
            double dE;
            int accepted = 0;
            Particle _old;

            _old.pos = particles[r].pos;
            _old.com = particles[r].com;
            _old.chargeDisp = particles[r].chargeDisp;
            _old.q = particles[r].q;
            _old.index = particles[r].index;
            
            eOld = energy_function(particles, particles[r]);

            particles[r].random_charge_rot();
            particles.update_distances(particles[r]);
            //energy::ewald3D::update_reciprocal(_old, particles[r]);

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
                //energy::ewald3D::update_reciprocal(particles[r], _old);
                particles[r].chargeDisp = _old.chargeDisp;
                particles[r].pos = _old.pos;
                particles.update_distances(particles[r]);
            }
            //delete _old;
            return accepted;
        }












        template <typename E>
        int charge_disp_move(E energy_function){
            double eOld;
            double eNew;
            double acceptProb;
            int r = ran2::get_random() * particles.numOfParticles;
            double dE;
            int accepted = 0;
            Particle _old = new Particle(true);

            _old.pos = particles[r].pos;
            _old.com = particles[r].com;
            _old.chargeDisp = particles[r].chargeDisp;
            _old.q = particles[r].q;
            _old.index = particles[r].index;
            _old.b = particles[r].b;

            //eOld = MC::direct.get_energy(particles);
            eOld = energy_function(particles, particles[r]);

            particles[r].b += (ran2::get_random() - 0.5) * 0.1;
            particles[r].chargeDisp = particles[r].chargeDisp.normalized() * particles[r].b;
            particles[r].pos = particles[r].com + particles[r].chargeDisp;

            particles.update_distances(particles[r]);
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
                particles[r].chargeDisp = _old.chargeDisp;
                particles[r].pos = _old.pos;
                particles[r].b = _old.b;
                particles.update_distances(particles[r]);
            }
            //delete _old;
            return accepted;
        }











        template<typename E>
        int trans_move(double dr, E energy_function){
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
            Particle _old;

            int p =  random * particles.numOfParticles;
            //printf("random: %i\n", p);

            eOld = energy_function(particles, particles[p]);
            _old.pos = particles[p].pos;
            _old.com = particles[p].com;
            _old.q = particles[p].q;
            _old.index = particles[p].index;

            //Generate new trial coordinates
            particles[p].random_move(dr);
            //energy::imgrep::update_position(particles, particles[p]);

            //If there is no overlap in new position and it's inside the box
            if(particles.hard_sphere(particles[p]) && particles[p].wall_2d()){

                //Get new energy
                //energy::ewald3D::update_reciprocal(_old, particles[p]);
                //energy::levin::update_f(_old, particles[p]);
                eNew = energy_function(particles, particles[p]);

                //Accept move?
                dE = eNew - eOld;
                acceptProp = exp(-1.0 * dE);
                if(acceptProp > 1.0 || eNew < eOld){
                    acceptProp = 1.0;
                }

                double rand = ran2::get_random();

                if(rand <= acceptProp){  //Accept move
                    Base::eCummulative += dE; //Update cummulative energy
                    accepted = 1;
                    Base::acceptedMoves++;
                }
                else{   //Reject move
                    //energy::imgrep::update_position(particles, _old);
                    //energy::ewald3D::update_reciprocal(particles[p], _old);
                    //energy::levin::update_f(particles[p], _old);
                    particles[p].pos = _old.pos;
                    particles[p].com = _old.com;
                }
            }

            else{   //Reject move
                particles[p].pos = _old.pos;
                particles[p].com = _old.com;
                //energy::imgrep::update_position(particles, particles[p]);
            }

            //delete _old;
            return accepted;
        }









/*
        template<typename E>
        int trans_electron_move(double dr, E energy_function){
            double eOld = 0.0;
            double eNew = 0.0;
            double acceptProp = 0.0;
            double random = ran2::get_random();
            double dE = 0.0;
            int accepted= 0;

            Particle _old = new Particle(true);

            int p =  random * particles.numOfElectrons;
            p += particles.numOfParticles;

            //Calculate old energy
            eOld = energy_function(particles.particles, particles[p]);
            _old.pos = particles[p].pos;
            _old.com = particles[p].com;
            _old.q = particles[p].q;
            _old.index = particles[p].index;

            //Generate new trial coordinates
            //particles[p].random_move(dr);
            particles[p].com[0] += (ran2::get_random() * 2.0 - 1.0) * dr;
            particles[p].com[1] += (ran2::get_random() * 2.0 - 1.0) * dr;
            Particle::pbc_xy(particles[p].com);
            particles[p].pos = particles[p].com;
            Particle::pbc_xy(particles[p].pos);
            Particle::update_distances(particles, particles[p]);

            //If there is no overlap in new position and it's inside the box
            //if((particles[p].com[2] > Base::wall + Base::zL && particles[p].com[2] < 2 * Base::zL - Base::wall) || 
            //    (particles[p].com[2] > Base::wall - Base::zL && particles[p].com[2] < Base::wall)){

                //Get new energy
                eNew = energy_function(particles.particles, particles[p]);

                //Accept move?
                dE = eNew - eOld;
                //acceptProp = exp(-1*dE);
                //if(acceptProp > 1 || eNew < eOld){
                //    acceptProp = 1;
                //}

                double rand = ran2::get_random();

                if(dE < 0){//(rand <= acceptProp){  //Accept move
                    Base::eCummulative += dE; //Update cummulative energy
                    accepted = 1;
                    Base::acceptedMoves++;
                }
                else{   //Reject move
                    particles[p].pos = _old.pos;
                    particles[p].com = _old.com;
                    Particle::update_distances(particles, particles[p]);
                }


            //delete _old;
            return accepted;
        }
*/


        double sampleSurfPot(double z){
            double energy = 0;
            particles.add( {ran2::get_random() * Base::xL, 
                            ran2::get_random() * Base::yL, 
                            z} );
            //std::cout << particles[particles.numOfParticles - 1].pos << std::endl;
            //printf("\n");
            energy = energy::valleau::get_particle_pot(particles, particles[particles.numOfParticles - 1]);
            particles.remove(particles.numOfParticles - 1);
            //std::cout << particles[particles.numOfParticles - 1].pos << std::endl;

            return energy;
        }





        template<typename F, typename FP>
        void run(F&& energy_function, FP&& particle_energy_function, double dr, int iter, bool sample, std::string outputFile){
            double energy_temp, z;
            std::vector<double> _surfpot(2);
            _surfpot[0] = 0.0;
            _surfpot[1] = 0.0;
            double partRatio = (double)particles.numOfParticles/(particles.numOfParticles + particles.numOfElectrons);
            std::vector< std::pair<double, double> > surfPotential;
            int prevAccepted = 0, wallSamp = 0, midSamp = 0;
            int rotAccepted = 0;
            int rotTot = 0;
            int volTot = 0;
            int transTot = 0;
            int volAccepted = 0;
            int transAccepted = 0;
            int gcTot = 0;
            int gcAcc = 0;
            int elAcc = 0;
            int elTot = 0;
            int outFreq = 10000;
            int k = 0;
            int catCreated = 0;
            int catDestroyed = 0;
            int anCreated = 0;
            int anDestroyed = 0;
            double random = 0;
            double rN = 1.0 / particles.numOfParticles;
            double rE = 1.0 / (particles.numOfElectrons);
            char volOut[40] = "volumes_\0";
            char histOut[40];
            strcpy(histOut, outputFile.c_str());

            Analysis *xHist = new Analysis(0.05, Base::xL);
            Analysis *yHist = new Analysis(0.05, Base::yL);
            Analysis *zHist = new Analysis(0.05, Base::zLBox);
            Analysis *surfpot = new Analysis(0.05, Base::zLBox);

            strcat(volOut, outputFile.c_str());

            char outName[] = ".txt";
            Base::eCummulative = energy_function(particles);
            particles.oldEnergy = Base::eCummulative;
            printf("\nRunning MC-loop at temperature: %lf, Bjerrum length is %lf\n\n", Base::T, Base::lB);
            for(int i = 0; i <= iter; i++){

                if(i % 100 == 0 && i >= 10000 && sample){
                    //rdf->sample_rdf(particles, histo, binWidth);
                    xHist->sampleHisto(particles, 0);
                    yHist->sampleHisto(particles, 1);
                    zHist->sampleHisto(particles, 2);
                    
                    if(ran2::get_random() <= 0.5){  //Sample middle of box
                        z = 0.0;
                        _surfpot[0] += sampleSurfPot(z);
                        midSamp++;
                    }
                    else{       //Sample at walls
                        if(ran2::get_random() <= 0.5){
                            z = Base::zLBox / 2.0;
                        }
                        else{
                            z = -1.0 * Base::zLBox / 2.0;
                        }
                        _surfpot[1] += sampleSurfPot(z);
                        wallSamp++;
                    }
                    

/*
                    double random = ran2::get_random();
                    bool success;
                    if(random < 0.25){
                        //Add to right wall
                        success = particles.add(Base::zLBox / 2.0);
                    }
                    else if(random < 0.5){
                        //Add to left wall
                        success = particles.add(-1.0 * Base::zLBox / 2.0);
                    }
                    else{
                        //Add to middle
                        success = particles.add(0.0);
                    }
                    if(success){
                        surfpot->sampleSurfPot(particle_energy_function(particles, particles[particles.numOfParticles]), random < 0.5);
                        particles.remove(particles.numOfParticles - 1);
                    }*/
                }
                
                /*random = ran2::get_random();
                if(random <= rN && i > 1000000){
                    if(vol_move(particles, energy_function)){
                        prevAccepted++;
                        volAccepted++;
                    } 
                    volTot++;
                }
                else if(random < 0.5 + rN){*/
                
                    random = ran2::get_random();
                    //if(random <= partRatio){
                    if(random <= 0.7){
                        //random = ran2::get_random();
                        //if(random <= rE){
                            if(trans_move(dr, particle_energy_function)){
                                prevAccepted++; 
                                transAccepted++;
                            }
                            transTot++;
                        //}
                        /*else{
                            if(trans_move(particles, Base::zL, particle_energy_function)){
                                prevAccepted++; 
                                transAccepted++;
                            }
                        }*/
                    }
                    else{
                        if(grand_move(particle_energy_function)){
                            gcAcc++;
                        }
                        gcTot++;
                    }
                    /*else{
                        if(trans_electron_move(particles, 0.05, particle_energy_function)){
                            prevAccepted++;
                            elAcc++;
                        }
                        elTot++;
                    }*/
                /*}
                else{
                    if(charge_rot_move(particle_energy_function)){
                        prevAccepted++;
                        rotAccepted++;
                    }
                    rotTot++;
                }*/
                Base::totalMoves++;
                

                
                if(i % 100 == 0 && i > 10000 && !sample){
                    energy::valleau::update_charge_vector(particles);
                }

                if(i % outFreq == 0){
                    //Base::volumes.push_back(Base::volume);
                    //Base::volumes[k] = Base::volume;
                    //k++;
                    energy_temp = energy_function(particles);
                    //particles.write_coordinates(outName , particles);
                    printf("Iteration: %d\n", i);
                    printf("Volume: %lf, x-dimension: %lf\n", Base::volume, Base::xL);
                    printf("Energy: %lf\n", energy_temp);
                    printf("Acceptance ratio: %lf\n", (double)Base::acceptedMoves/Base::totalMoves);
                    printf("Acceptance ratio for the last %i steps: %lf\n", outFreq, (double)prevAccepted/outFreq);
                    if(std::abs(energy_temp - Base::eCummulative)/std::abs(energy_temp) > 1E-9){
                        printf("Error is too large!\n");
                        printf("Error: %.12lf\n", std::abs(energy_temp - Base::eCummulative)/std::abs(energy_temp));
                        exit(1);
                    }
                    printf("Error: %.12lf\n", std::abs(energy_temp - Base::eCummulative)/std::abs(energy_temp));
                    printf("Trans mv: %d, %.1lf  Rot mv: %d, %.1lf  Vol mv: %d, %.1lf GC mv: %d, %.1lf pnum: %i  nnum: %i surfPot: %lf wpot %lf mpot %lf\n\n", 
                                                                                                            transAccepted, (double) transAccepted/transTot * 100.0,
                                                                                                            rotAccepted,   (double) rotAccepted/rotTot * 100.0, 
                                                                                                            volAccepted,   (double) volAccepted/volTot * 100.0,
                                                                                                            gcAcc,   (double) gcAcc/gcTot * 100.0,
                                                                                                           particles.numOfCations, particles.numOfAnions,
                                                                                                           (_surfpot[1] / ((double) wallSamp) - _surfpot[0] / ((double) midSamp) ),
                                                                                                           _surfpot[1] / ((double) wallSamp),
                                                                                                           _surfpot[0] / ((double) midSamp));
                    
                    prevAccepted = 0;
                    //printf("size: %lu\n", Base::volumes.size());
                    /*
                    if(k == 1000){
                        FILE *f = fopen(volOut, "a");
                        fprintf(f, "");
                        if(f == NULL){
                            printf("Can't open file!\n");
                            exit(1);
                        }
                        for(int j = 0; j < k; j++){
                            fprintf(f, "%d      %lf\n",(i / outFreq) + j, Base::volumes[j]);
                        }
                        fclose(f);
                        k = 0;
                        //Base::volumes.clear();
                    }
                    */
                }
            }

            if(sample){
                xHist->saveHisto(histOut, particles);
                yHist->saveHisto(histOut, particles);
                zHist->saveHisto(histOut, particles);
                surfpot->saveSurfPot(histOut);
            }

            delete xHist;
            delete yHist;
            delete zHist;
        }
};

#endif
