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

//Globals
double xL = 14;
double yL = 14;
double zL = 14;
int numOfParticles = 1000;
int totalMoves = 0;
int acceptedMoves = 0;
int numberOfSamples = 0;
double kb = 1;
double T = 1.5;
double beta = 1/(1*1.5);
double eCum = 0;

double get_random(){
    int ran_input = -1*(int) time(NULL)*100*rand();
    float ran2_var = ran2_(&ran_input);
    return (double)ran2_var;
}

Particle* pbc(Particle *p){
    //Translate particles according to periodic boundary conditions
    if(p->pos[0] > xL){
        p->pos[0] = p->pos[0] - xL;
    }
    if(p->pos[0] < 0){
        p->pos[0] = p->pos[0] + xL;
    }
    if(p->pos[1] > yL){
        p->pos[1] = p->pos[1] - yL;
    }
    if(p->pos[1] < 0){
        p->pos[1] = p->pos[1] + yL;
    }
    if(p->pos[2] < 0){
        p->pos[2] = p->pos[2] + zL;
    }
    if(p->pos[2] > zL){
        p->pos[2] = p->pos[2] - zL;
    }

    return p;
}

Particle* randomMove(Particle *p){
    //Move particle randomly
    p->pos[0] += (double) rand()/RAND_MAX*0.5 - 0.25;
    p->pos[1] += (double) rand()/RAND_MAX*0.5 - 0.25;
    p->pos[2] += (double) rand()/RAND_MAX*0.5 - 0.25;
    pbc(p);

    return p;
}

double distance(Particle *p1, Particle *p2){
    //Calculate distance between particles
    double xP1 = p1->pos[0];
    double yP1 = p1->pos[1];
    double zP1 = p1->pos[2];
    double xP2 = p2->pos[0];
    double yP2 = p2->pos[1];
    double zP2 = p2->pos[2];

    if(xP1 - xP2 < -1 * xL/2){
        xP2 = xP2 - xL;
    }
    if(xP1 - xP2 > xL/2){
        xP2 = xP2 + xL;
    }
    if(yP1 - yP2 < -1 * yL/2){
        yP2 = yP2 - yL;
    }
    if(yP1 - yP2 > yL/2){
        yP2 = yP2 + yL;
    }
    if(zP1 - zP2 < -1 * zL/2){
        zP2 = zP2 - zL;
    }
    if(zP1 - zP2 > zL/2){
        zP2 = zP2 + zL;
    }
    return (pow((xP1 - xP2), 2) + pow((yP1 - yP2), 2) + pow((zP1 - zP2), 2));
}

double getEnergy(Particle **particles){
    int i = 0;
    int j = 0;
    double dist = 0;
    double energy = 0;
    double r6 = 0;
    double sigma = 1;
    double epsilon = 1;

    for(i = 0; i < numOfParticles; i++){
        j = i + 1;
        while(j < numOfParticles){
            //LJ
            dist = distance(particles[j], particles[i]);
            r6 = pow(dist, 3);
            energy += sigma/(r6*r6) - sigma/r6;
            //Coloumb
            dist = sqrt(dist);
            energy += 1/(4 * PI * 1 * 2 * T) * (particles[i]->q * particles[j]->q)/dist;
            j++;
        }
    }
    return 4 * epsilon * energy;
}

double getParticleEnergy(int pInd, Particle *p, Particle **particles){
    int i = 0;
    double energy = 0;
    double dist = 0;
    double r6 = 0;
    double sigma = 1;
    double epsilon = 1;

    for(i = 0; i < numOfParticles; i++){
        if(i != pInd){
            //LJ
            dist = distance(p, particles[i]);
            r6 = pow(dist, 3);
            energy += sigma/(r6*r6) - sigma/r6;
            //printf("Lennard: %lf\n", energy);
            //Coloumb
            dist = sqrt(dist);
            energy += 1/(4 * PI * 1 * 2 * T) * (p->q * particles[i]->q)/dist;
            //printf("Coloumb: %lf\n", energy);
        }
    }
    return 4 * epsilon * energy;
}

int mcmove(Particle **particles){
    int i = 0;
    double eOld = 0;
    double eNew = 0;
    double dist = 0;
    double acceptProp = 0;
    double random = get_random();
    int p =  random * numOfParticles;
    double dr = 0.5;
    double dE = 0;
    int accepted= 0;
    double oldPos[3] = {particles[p]->pos[0], particles[p]->pos[1], particles[p]->pos[2]};
    beta = 1/(1*T);

    //Calculate old energy
    eOld = getParticleEnergy(p, particles[p], particles);
    //printf("Old energy: %lf\n", eOld);
    //Generate new trial coordinates
    particles[p]->pos[0] = particles[p]->pos[0] + (get_random()*2 - 1) * dr;
    particles[p]->pos[1] = particles[p]->pos[1] + (get_random()*2 - 1) * dr;
    particles[p]->pos[2] = particles[p]->pos[2] + (get_random()*2 - 1) * dr;

    //Appy PBC
    pbc(particles[p]);

    //Get new energy
    eNew = getParticleEnergy(p, particles[p], particles);
    //printf("New energy: %lf\n", eNew);

    //Accept move?
    dE = eNew - eOld;
    acceptProp = exp(-1*beta*dE);

    if(acceptProp > 1 || eNew < eOld){
        acceptProp = 1;
    }
    //printf("AcceptProp: %lf\n", acceptProp);

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
    return accepted;
}

void sampleHisto(Particle **particles, int *histo, double binWidth){
    int i = 0;
    int j = 0;
    double dist = 0;
    double xL2 = (xL * xL)/4;

    for(i = 0; i < numOfParticles; i++){
        j = i + 1;
        while(j < numOfParticles){
            dist = distance(particles[i], particles[j]);
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
        idealDen = dv * numOfParticles/(xL*yL*zL);
        fprintf(f, "%lf     %lf\n", i * binWidth, histo[i]/(numOfParticles * numberOfSamples * idealDen));
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
    beta = 1/(kb * T);
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
        particles[i]->d = 1.0;

        if(i % 2 == 0){
            particles[i]->q = -1.0;
            strcpy(particles[i]->name, "Cl\0");
        }
        else{
            particles[i]->q = 1.0;
            strcpy(particles[i]->name, "Na\0");

        }

        pbc(particles[i]);
    }
    
    //Calculate all distances
    for(i = 0; i < numOfParticles; i++){
        for(j = 0; j < numOfParticles; j++){
            if(i != j){
                if(distance(particles[i], particles[j]) < diameter2){
                    overlaps++;
                }
            }
        }
    }

    printf("Overlaps before moves: %d\n", overlaps);

    //Move particles to prevent overlap
    while(k < 100000){
        overlaps = 0;
        for(i = 0; i < numOfParticles; i++){
            for(j = 0; j < numOfParticles; j++){
                if(i != j){
                    if(distance(particles[i], particles[j]) < diameter2){
                        randomMove(particles[i]);
                        overlaps += 1;
                    }
                }
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

    // overlaps = 0;
    // for(i = 0; i < numOfParticles; i++){
    //     for(j = 0; j < numOfParticles; j++){
    //         if(i != j){
    //             if(distance(particles[i], particles[j]) < diameter2){
    //                 overlaps++;
    //             }
    //         }
    //     }
    // }

    //printf("Overlaps after moves: %d\n", overlaps);

    //int raninput = -1*(int) time(NULL);
    //float ran2_var = ran2_(&raninput);
    //printf("%f\n", ran2_var);

    //Write coordinates to file
    FILE *fi = fopen("output_before.xyz", "w");
    if(fi == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    fprintf(fi, "999\n");
    for(i = 0; i < numOfParticles; i++){
        fprintf(fi, "%s     %lf    %lf     %lf\n", particles[i]->name, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    }
    fclose(fi);

    //Update cumulative energy
    eCum = getEnergy(particles);

    // T = 0.3;
    // beta = 1/(10*T);
    // ///////////////////////////////         Main MC-loop          ////////////////////////////////////////
    printf("\nRunning main MC-loop at temperature: %lf\n", T);
    int prevAccepted = 0;
    for(i = 0; i < 500000; i++){

        // if(i % 1000000 == 0 && i != 0){
        //     if(i < 2000001){
        //         T -= 0.1;
        //         printf("Setting temperature to: %lf\n", T);
        //         beta = 1/(10*T);
        //     }
        // }

        if(i % 10000 == 0 && i > 100000){
            sampleHisto(particles, histo, binWidth);
        }

        if(mcmove(particles)){
           prevAccepted++; 
        }
        totalMoves++;
        //printf("%lf\n", fabs(getEnergy(pCoo) - eCum)/eCum);
        
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
    for(i = 0; i < numOfParticles; i++){
        fprintf(fi, "%s     %lf    %lf     %lf\n", particles[i]->name, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    }
    fclose(f);

    //Clean up allocated memory
    for(i = 0; i < numOfParticles; i++){
        free(particles[i]->pos);
        free(particles[i]);
    }
    free(particles);
    free(histo);
   return 0;
}