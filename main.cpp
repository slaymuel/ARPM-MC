//Create particles
//Check overlap
//Move particle, accept move if no overlap
#include "particle.h"
#include "constants.h"
#include "analysis.h"
#include "ran2_lib.cpp"
#include "mc.h"
#include "ewald3D.h"
//Initializers

Ewald3D MC::ewald;
Direct MC::direct;

int Particle::numOfParticles = 0;
int Analysis::numOfHisto = 0;
double Base::xL= 10;
double Base::yL= 10;
double Base::zL= 10;
double Base::T = 2000;
double Base::eCummulative = 0;
int Base::acceptedMoves = 0;
int Base::totalMoves = 0;

int numOfParticles = 2;//6728;
int totalMoves = 0;
int acceptedMoves = 0;
int numberOfSamples = 0;
double kb = KB;
double T = Base::T;
double beta = 1;

int main(int argc, char *argv[])
{
    if( argc == 2 ) {
        printf("input: %s\n", argv[1]);
      if(!strcmp(argv[1], "rdf")){
          //readCoo();
          exit(1);
      }
    }

    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int prevAccepted = 0;
    int Digs = 14;
    int bins = 100;
    double binWidth = (Base::xL/2)/bins;
    double dist = 0;
    int *histo;
    histo = (int*)malloc(bins * sizeof(int));

    Analysis *xHist = new Analysis(0.1, Base::xL);
    Analysis *yHist = new Analysis(0.1, Base::yL);
    Analysis *zHist = new Analysis(0.1, Base::zL);

    //Particle variables
    double diameter = 5;
    double diameter2 = pow(diameter, 2);
    int overlaps = 0;
    double density = (double)numOfParticles/(Base::xL * Base::yL * Base::zL) * pow(diameter, 3);
    int overlap = 1;
    double energy = 0;

    //Seed
    srand(time(NULL));

    printf("Density is: %lf\n", density);

    Particle **particles;
    particles = (Particle**) malloc(numOfParticles * sizeof(Particle*));
    particles[0] = new Particle();
    particles[1] = new Particle();
    particles[0]->pos = (double*)malloc(3*sizeof(double));
    particles[1]->pos = (double*)malloc(3*sizeof(double));

    particles[0]->pos[0] = 0;
    particles[0]->pos[1] = 0;
    particles[0]->pos[2] = 0;
    particles[1]->pos[0] = 1;
    particles[1]->pos[1] = 0;
    particles[1]->pos[2] = 0;

    particles[0]->d = 5;
    particles[0]->index = 0;
    particles[0]->q = -1.0;
    strcpy(particles[0]->name, "Cl\0");
    particles[1]->d = 5;
    particles[1]->index = 1;
    particles[1]->q = 1.0;
    strcpy(particles[1]->name, "Na\0");
    //char filename[] = "output_sim_newest.xyz";
    //particles = Particle::read_coordinates(filename);
    //particles = Particle::create_particles(numOfParticles);
    //printf("Outside: %lf\n", particles[0]->d);

    MC mc;
    MC::ewald.initialize();
    //mc.equilibrate(particles);
    printf("Initialized k-vectors.\n");
    overlaps = Particle::get_overlaps(particles);
    printf("Overlaps: %d\n", overlaps);
    // char name[] = "output_equilibrate_newest.xyz";
    // Particle::write_coordinates(name , particles);


    //Update cumulative energy
    Base::eCummulative = mc.getEnergy(particles);
    // ///////////////////////////////         Main MC-loop          ////////////////////////////////////////
    printf("\nRunning main MC-loop at temperature: %lf\n", T);
    for(i = 0; i < 1000000; i++){
        if(i % 1000 == 0 && i > 100){
            //sampleRDF(particles, zHisto, binWidth);
            //printf("Sampling...\n");
            // xHist->sampleHisto(particles, 0);
            // yHist->sampleHisto(particles, 1);
            // zHist->sampleHisto(particles, 2);
        }
        if(i % 100 == 0){
            if(mc.mcmove(particles, Base::zL)){
                prevAccepted++; 
            }
        }
        else{
            if(mc.mcmove(particles, 0.8)){
                prevAccepted++; 
            }    
        }

        Base::totalMoves++;
        
        if(i % 100000 == 0 && i != 0){
            energy = mc.getEnergy(particles);
            printf("Iteration: %d\n", i);
            printf("Energy: %lf\n", energy);
            printf("Error: ");
            printf("%.*e\n", Digs, fabs(energy - Base::eCummulative)/fabs(Base::eCummulative));
            printf("Acceptance ratio: %lf\n", (double)Base::acceptedMoves/Base::totalMoves);
            printf("Acceptance ratio for the last 100000 steps: %lf\n\n", (double)prevAccepted/100000.0);
            if(fabs(energy - Base::eCummulative)/fabs(energy) > pow(10, -12)){
                printf("Error is too large!\n");
                exit(1);
            }
            prevAccepted = 0;
        }
        
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    printf("Accepted moves: %d\n", Base::acceptedMoves);
    printf("Rejected moves: %d\n", Base::totalMoves - Base::acceptedMoves);

    //saveRDF(histo, bins, binWidth);
    // xHist->saveHisto();
    // yHist->saveHisto();
    // zHist->saveHisto();

    //Write coordinates to file
    // FILE *f = fopen("output_sim_newest.xyz", "w");
    // if(f == NULL){
    //     printf("Can't open file!\n");
    //     exit(1);
    // }

    // fprintf(f, "%d\n", numOfParticles);
    // fprintf(f, "\n");
    // for(i = 0; i < Particle::numOfParticles; i++){
    //     fprintf(f, "%s     %lf    %lf     %lf\n", particles[i]->name, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    // }
    // fclose(f);

    //Clean up allocated memory
    printf("Cleaning up...\n");
    for(i = 0; i < Particle::numOfParticles; i++){
        free(particles[i]->pos);
        free(particles[i]);
    }
    //for(i = 0; i < numOfCells; i++){
    //    free(grid[i]);
    //}
    //free(grid);
    free(particles);
    free(histo);
   return 0;
}