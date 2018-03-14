#include "analysis.h"

Analysis::Analysis(double binWidth){
    numberOfSamples = 0;
    this->binWidth = binWidth;
    bins = zL/binWidth;
    histo = (int*)malloc(bins*sizeof(int));
    pHisto = (int*)malloc(bins*sizeof(int));
    nHisto = (int*)malloc(bins*sizeof(int));
    num = numOfHisto;
    numOfHisto++;
}
void Analysis::sampleHisto(Particle **particles, int d){
    int i = 0;
    double dist = 0;
    //printf("%d\n", Particle::numOfParticles);
    for(i = 0; i < Particle::numOfParticles; i++){
        histo[(int)(particles[i]->pos[d]/binWidth)]++;
        //printf("Current histo value: %d\n", histo[(int)(particles[i]->pos[d]/binWidth)]);
        if(particles[i]->q < 0){
            nHisto[(int)(particles[i]->pos[d]/binWidth)]++;
            //printf("Current histo value: %d\n", nHisto[(int)(particles[i]->pos[d]/binWidth)]);
        }
        else{
            pHisto[(int)(particles[i]->pos[d]/binWidth)]++;
            //printf("Current histo value: %d\n", pHisto[(int)(particles[i]->pos[d]/binWidth)]);
        }
    }
    numberOfSamples++;
}

void Analysis::saveHisto(){
    int i = 0;
    double dv = 0;
    double idealDen = 0;
    char pHisto_name[64];
    sprintf(pHisto_name, "pHisto_%d.txt", num);
    char nHisto_name[64];
    sprintf(nHisto_name, "nHisto_%d.txt", num);
    char histo_name[64];
    sprintf(histo_name, "histo_%d.txt", num);

    FILE *f = fopen(histo_name, "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    for(i = 0; i < bins; i++){
        fprintf(f, "%lf     %lf\n", (double)i * binWidth, (double)histo[i]/(Particle::numOfParticles * numberOfSamples));
    }
    fclose(f);

    FILE *pf = fopen(pHisto_name, "w");
    if(pf == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    for(i = 0; i < bins; i++){
        fprintf(pf, "%lf     %lf\n", (double)i * binWidth, (double)pHisto[i]/(Particle::numOfParticles * numberOfSamples));
    }
    fclose(pf);

    FILE *nf = fopen(nHisto_name, "w");
    if(nf == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    for(i = 0; i < bins; i++){
        fprintf(nf, "%lf     %lf\n", (double)i * binWidth, (double)nHisto[i]/(Particle::numOfParticles * numberOfSamples));
    }
    fclose(f);
}
