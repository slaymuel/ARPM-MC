#include "mc.h"
#include <vector>

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
 /**
   * @brief Approximation of erfc-function
   * @param x Value for which erfc should be calculated 
   * @details Reference for this approximation is found in Abramowitz and Stegun, 
   *          Handbook of mathematical functions, eq. 7.1.26
   *
   * @f[
   *     \erf(x) = 1 - (a_1t + a_2t^2 + a_3t^3 + a_4t^4 + a_5t^5)e^{-x^2} + \epsilon(x)
   * @f]
   * @f[
   *     t = \frac{1}{1 + px}
   * @f]
   * @f[
   *     |\epsilon(x)| \le 1.5\times 10^{-7}
   * @f]
   * 
   * @warning Needs testing for x < 0.
   */
  template<typename T>
  T erfc_x( T x )
  {
    static_assert(std::is_floating_point<T>::value, "type must be floating point");
    if(x < 0)
    return ( 2.0 - erfc_x(-x) );
    T t = 1.0 / (1.0 + 0.3275911 * x);
    const T a1 = 0.254829592;
    const T a2 = -0.284496736;
    const T a3 = 1.421413741;
    const T a4 = -1.453152027;
    const T a5 = 1.061405429;
    return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5)))) * exp(-x * x);
}

/**
 * @brief Approximate 1 - erfc_x
 * @param x Value for which erf should be calculated 
 */
template<typename T>
T erf_x( T x ) { 
    return (1 - erfc_x(x)); 
}

template<typename T>
T ewald_F(T x){

}

double ewald(Particle *p, Particle **particles){
    int i = 0;
    double energy = 0;
    double alpha = 5/Base::xL;
    double r = 0;
    double qq = 0;
    int kx = 0;
    double kx2;
    int ky = 0;
    double ky2;
    int kz = 0;
    double kz2;
    double k2 = 0;
    int kNumMax = 400;
    int kNum= 0;
    double kMax = 8*PI/Base::xL;
    std::vector< std::vector<double> > kVec;

    for(kx = 0; kx < kMax; kx++){
        for(int ky = -kMax + 1; ky < kNumMax; ky++){
            kx2 = kx*kx;
            ky2 = ky*ky;
            k2 = kx2 + ky2;
            if(k2 != 0 && k2 < kMax){
                std::vector<double> vec;
                vec.push_back(2*PI*kx/Base::xL);
                vec.push_back(2*PI*ky/Base::yL);
                vec.push_back(0);
                kVec.push_back(vec);
                //kVec[kNum] = exp(-k2/());
                kNum++;
            }
        }
    }
    for(i = 0; i < Particle::numOfParticles; i++){
        qq = p->q * particles[i]->q;
        r = p->distance_xy(particles[i]);
        energy += erfc_x(alpha*r)/r;
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
        if(!particles[p]->hardSphere(particles) || particles[p]->pos[2] < particles[p]->d/2 || particles[p]->pos[2] > Base::zL - particles[p]->d/2){
            particles[p]->pos[0] = oldPos[0];
            particles[p]->pos[1] = oldPos[1];
            particles[p]->pos[2] = oldPos[2];
        }
        if(i % 50000 == 0){
            overlaps = Particle::get_overlaps(particles);
            //stepSize = log(abs(prevOverlaps - overlaps) + 3.0);
            //prevOverlaps = overlaps;
            printf("Overlaps: %d\r", overlaps);
            fflush(stdout);
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