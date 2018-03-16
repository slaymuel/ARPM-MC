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


// void sampleRDF(Particle **particles, int *histo, double binWidth){
//     int i = 0;
//     int j = 0;
//     double dist = 0;
//     double xL2 = (Base::xL * Base::xL)/4;

//     for(i = 0; i < Particle::numOfParticles; i++){
//         j = i + 1;
//         while(j < Particle::numOfParticles){
//             dist = particles[i]->distance(particles[j]);
//             if(dist < xL2){
//                 dist = sqrt(dist);
//                 histo[(int)(dist/binWidth)] = histo[(int)(dist/binWidth)] + 2;
//             }
//             j++;
//         }
//     }
//     numberOfSamples++;
// }

// void saveRDF(int *histo, int bins, double binWidth){
//     int i = 0;
//     double dv = 0;
//     double idealDen = 0;

//     FILE *f = fopen("histo.txt", "w");
//     if(f == NULL){
//         printf("Can't open file!\n");
//         exit(1);
//     }

//     //Number of particles in ideal gas with same density
//     for(i = 0; i < bins; i++){
//         dv = (double)4/3 * PI * (pow(i + 1, 3) - pow(i, 3)) * pow(binWidth, 3);
//         idealDen = dv * Particle::numOfParticles/(Base::xL*Base::yL*Base::zL);
//         fprintf(f, "%lf     %lf\n", i * binWidth, histo[i]/(Particle::numOfParticles * numberOfSamples * idealDen));
//     }
//     fclose(f);
// }