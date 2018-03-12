#include "analysis.h"

Analysis::Analysis(double binWidth){
    numberOfSamples = 0;
    this->binWidth = binWidth;
    bins = zL/binWidth;
    histo = (int*)malloc(bins*sizeof(int));
    pHisto = (int*)malloc(bins*sizeof(int));
    nHisto = (int*)malloc(bins*sizeof(int));
}
void Analysis::sampleHisto(Particle **particles, int d){
    int i = 0;
    double dist = 0;

    for(i = 0; i < Particle::numOfParticles; i++){
        histo[(int)(particles[i]->pos[d]/binWidth)]++;
        if(particles[i]->q < 0){
            nHisto[(int)(particles[i]->pos[d]/binWidth)]++;
        }
        else{
            pHisto[(int)(particles[i]->pos[d]/binWidth)]++;
        }
    }
    numberOfSamples++;
}

void Analysis::saveHisto(){
    int i = 0;
    double dv = 0;
    double idealDen = 0;

    FILE *f = fopen("histo_z.txt", "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    //Number of particles in ideal gas with same density
    for(i = 0; i < bins; i++){
        fprintf(f, "%lf     %lf\n", (double)i * binWidth, (double)histo[i]/(Particle::numOfParticles * numberOfSamples));
    }
    fclose(f);

    FILE *pf = fopen("pHisto_z.txt", "w");
    if(pf == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    //Number of particles in ideal gas with same density
    for(i = 0; i < bins; i++){
        fprintf(pf, "%lf     %lf\n", (double)i * binWidth, (double)pHisto[i]/(Particle::numOfParticles * numberOfSamples));
    }
    fclose(pf);

    FILE *nf = fopen("nHisto_z.txt", "w");
    if(nf == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    //Number of particles in ideal gas with same density
    for(i = 0; i < bins; i++){
        fprintf(nf, "%lf     %lf\n", (double)i * binWidth, (double)nHisto[i]/(Particle::numOfParticles * numberOfSamples));
    }
    fclose(f);
}
