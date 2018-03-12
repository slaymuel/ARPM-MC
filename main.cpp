//Create particles
//Check overlap
//Move particle, accept move if no overlap
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "particle.h"
#include "constants.h"

//Fortran ran2
//extern float ran2_();
extern "C"
{
    float ran2_(int*);
}

//Initializers
int Particle::numOfParticles = 0;

double Base::xL= 14;
double Base::yL= 14;
double Base::zL= 14;

//Globals
double xL = 14;
double yL = 14;
double zL = 14;
int numOfParticles = 1000;
int totalMoves = 0;
int acceptedMoves = 0;
int numberOfSamples = 0;
double kb = KB;
double T = 1;
double beta = 1/(1*T);
double eCum = 0;

double get_random(){
    int ran_input = -1*(int) time(NULL)*100*rand();
    float ran2_var = ran2_(&ran_input);
    return (double)ran2_var;
}

double getEnergy(Particle **particles){
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
            energy += pow(EC, 2)/(4 * VP * PI * 2 * kb * T * 1e-10) * (particles[i]->q * particles[j]->q)/dist;
            j++;
        }
    }
    return energy;
}

double getParticleEnergy(int pInd, Particle *p, Particle **particles){
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
            energy += pow(EC, 2)/(4 * VP * PI * 2 * kb * T * 1e-10) * (p->q * particles[i]->q)/dist;
            //printf("Coloumb: %lf\n", energy);
        }
    }
    return energy;
}

int mcmove(Particle **particles){
    int i = 0;
    double eOld = 0;
    double eNew = 0;
    double dist = 0;
    double acceptProp = 0;
    double random = get_random();
    double dr = 0.5;
    double dE = 0;
    int accepted= 0;

    int p =  random * numOfParticles;

    //Calculate old energy
    eOld = getParticleEnergy(p, particles[p], particles);

    //Generate new trial coordinates
    double oldPos[3] = {particles[p]->pos[0], particles[p]->pos[1], particles[p]->pos[2]};
    particles[p]->pos[0] = particles[p]->pos[0] + (get_random()*2 - 1) * dr;
    particles[p]->pos[1] = particles[p]->pos[1] + (get_random()*2 - 1) * dr;
    particles[p]->pos[2] = particles[p]->pos[2] + (get_random()*2 - 1) * dr;

    //Appy PBC
    particles[p]->pbc_xy();
    if(particles[p]->hardSphere(particles) && particles[p]->pos[2] > particles[p]->d/2 && particles[p]->pos[2] < Base::zL - particles[p]->d/2 ){
        //Get new energy
        eNew = getParticleEnergy(p, particles[p], particles);

        //Accept move?
        dE = eNew - eOld;
        acceptProp = exp(-1*beta*dE);

        if(acceptProp > 1 || eNew < eOld){
            acceptProp = 1;
        }

        double rand = get_random();
        if(rand < acceptProp){
            eCum += dE;
            accepted = 1;
            acceptedMoves++;
        }
        else{
            particles[p]->pos[0] = oldPos[0];
            particles[p]->pos[1] = oldPos[1];
            particles[p]->pos[2] = oldPos[2];
        }
    }

    else{
        particles[p]->pos[0] = oldPos[0];
        particles[p]->pos[1] = oldPos[1];
        particles[p]->pos[2] = oldPos[2];
    }
    return accepted;
}

void sampleHisto(Particle **particles, int *histo, double binWidth){
    int i = 0;
    int j = 0;
    double dist = 0;
    double xL2 = (xL * xL)/4;

    for(i = 0; i < Particle::numOfParticles; i++){
        j = i + 1;
        while(j < Particle::numOfParticles){
            dist = particles[i]->distance(particles[j]);
            if(dist < xL2){
                dist = sqrt(dist);
                histo[(int)(dist/binWidth)] = histo[(int)(dist/binWidth)] + 2;
            }
            j++;
        }
    }
    numberOfSamples++;
}

void saveRDF(int *histo, int bins, double binWidth){
    int i = 0;
    double dv = 0;
    double idealDen = 0;

    FILE *f = fopen("histo.txt", "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    //Number of particles in ideal gas with same density
    for(i = 0; i < bins; i++){
        dv = (double)4/3 * PI * (pow(i + 1, 3) - pow(i, 3)) * pow(binWidth, 3);
        idealDen = dv * Particle::numOfParticles/(xL*yL*zL);
        fprintf(f, "%lf     %lf\n", i * binWidth, histo[i]/(Particle::numOfParticles * numberOfSamples * idealDen));
    }
    fclose(f);
}

void readCoo(){
    int c;
    double a;

    FILE *file;
    file = fopen("output.xyz", "r");
    if (file) {
        while ((c = getc(file)) != EOF)
            fscanf(file, "%lf \n", &a);
            //putchar(c);
            printf("%lf", a);
        fclose(file);
    }
}

int main(int argc, char *argv[])
{
    if( argc == 2 ) {
        printf("input: %s\n", argv[1]);
      if(!strcmp(argv[1], "rdf")){
          readCoo();
          exit(1);
      }
    }
    int i = 0;
    int j = 0;
    int k = 0;
    int Digs = DECIMAL_DIG;
    int bins = 100;
    double binWidth = (xL/2)/bins;
    double dist = 0;
    int *histo;
    histo = (int*)malloc(bins * sizeof(int));

    //Particle variables
    double diameter = 1;
    double diameter2 = pow(diameter, 2);
    int overlaps = 0;
    double density = (double)numOfParticles/(xL*yL*zL)*pow(diameter, 2);
    int overlap = 1;
    //Seed
    srand(time(NULL));

    printf("Density is: %lf\n", density);

    //Pointers to hold particle objects
    Particle **particles = (Particle**) malloc(numOfParticles * sizeof(Particle));

    //Set random initial coordiantes
    for(i = 0; i < numOfParticles; i++){
        particles[i] = new Particle();
        particles[i]->pos = (double*)malloc(3*sizeof(double));
        particles[i]->pos[0] = (double) rand()/RAND_MAX * xL;
        particles[i]->pos[1] = (double) rand()/RAND_MAX * yL;
        particles[i]->pos[2] = (double) rand()/RAND_MAX * zL;
        particles[i]->d = diameter;
        particles[i]->index = i;

        if(i % 2 == 0){
            particles[i]->q = -1.0;
            strcpy(particles[i]->name, "Cl\0");
        }
        else{
            particles[i]->q = 1.0;
            strcpy(particles[i]->name, "Na\0");
        }
        //pbc(particles[i]);
        particles[i]->pbc();
    }
    
    //Calculate all distances
    for(i = 0; i < Particle::numOfParticles; i++){
        j = i + 1;
        while(j < Particle::numOfParticles){
            if(i != j){
                if(particles[i]->distance(particles[j]) < diameter2){
                    overlaps++;
                }
            }
            j++;
        }
    }

    printf("Overlaps before moves: %d\n", overlaps);

    //Move particles to prevent overlap
    while(k < 100000){
        overlaps = 0;
        for(i = 0; i < Particle::numOfParticles; i++){
            j = i + 1;
            while(j < Particle::numOfParticles){
                if(i != j){
                    if(particles[i]->distance(particles[j]) < diameter2){
                        //randomMove(particles[i]);
                        particles[i]->randomMove();
                        overlaps += 1;
                    }
                }
                j++;
            }
        }
        if(overlaps == 0){
            break;
        }
        k += 1;   
    }

    if(overlaps != 0){
        printf("%d overlaps persisted after random moves...\n", overlaps);
        exit(1);
    }
    printf("It took %d iterations to remove overlaps\n", k);

    //Write coordinates to file
    FILE *fi = fopen("output_before.xyz", "w");
    if(fi == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    fprintf(fi, "999\n");
    for(i = 0; i < Particle::numOfParticles; i++){
        fprintf(fi, "%s     %lf    %lf     %lf\n", particles[i]->name, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    }
    fclose(fi);

    //Update cumulative energy
    eCum = getEnergy(particles);

    // ///////////////////////////////         Main MC-loop          ////////////////////////////////////////
    printf("\nRunning main MC-loop at temperature: %lf\n", T);
    int prevAccepted = 0;
    for(i = 0; i < 1000000; i++){

        if(i % 10000 == 0 && i > 100000){
            sampleHisto(particles, histo, binWidth);
        }

        if(mcmove(particles)){
           prevAccepted++; 
        }
        totalMoves++;
        
        if(i % 100000 == 0 && i != 0){
            printf("Iteration: %d\n", i);
            printf("Energy: %lf\n", getEnergy(particles));
            printf("Error: ");
            printf("%.*e\n", Digs, fabs(getEnergy(particles) - eCum)/fabs(eCum));
            printf("Acceptance ratio: %lf\n", (double)acceptedMoves/totalMoves);
            printf("Acceptance ratio for the last 100000 steps: %lf\n\n", (double)prevAccepted/100000.0);
            if(fabs(getEnergy(particles) - eCum)/fabs(eCum) > pow(10, -12)){
                printf("Error is too large!\n");
                exit(1);
            }
            prevAccepted = 0;
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("Accepted moves: %d\n", acceptedMoves);
    printf("Rejected moves: %d\n", totalMoves - acceptedMoves);

    saveRDF(histo, bins, binWidth);

    //Write coordinates to file
    FILE *f = fopen("output.xyz", "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    fprintf(fi, "999\n");
    for(i = 0; i < Particle::numOfParticles; i++){
        fprintf(fi, "%s     %lf    %lf     %lf\n", particles[i]->name, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    }
    fclose(f);

    //Clean up allocated memory
    for(i = 0; i < Particle::numOfParticles; i++){
        free(particles[i]->pos);
        free(particles[i]);
    }
    free(particles);
    free(histo);
   return 0;
}