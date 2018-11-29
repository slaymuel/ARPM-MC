#include "analysis.h"

Analysis::Analysis(double binWidth, double dLength){
    this->numberOfSamples = 0;
    this->binWidth = binWidth;
    this->bins = dLength/binWidth;//zL/binWidth;
    this->histo = (int*)malloc(this->bins * sizeof(int));
    this->pHisto = (int*)malloc(this->bins * sizeof(int));
    this->nHisto = (int*)malloc(this->bins * sizeof(int));

    for(int i = 0; i < this->bins; i++){
        this->histo[i] = 0;
        this->pHisto[i] = 0;
        this->nHisto[i] = 0;
    }
    this->num = numOfHisto;
    this->numOfHisto++;
}

void Analysis::sampleHisto(Particle **particles, int d){
    int i = 0;
    double dist = 0;
    for(i = 0; i < Particle::numOfParticles; i++){
        this->histo[(int)(particles[i]->pos[d]/binWidth)]++;
        if(particles[i]->q < 0){
            this->nHisto[(int)(particles[i]->pos[d]/binWidth)]++;
        }
        else if(particles[i]->q > 0){
            this->pHisto[(int)(particles[i]->pos[d]/binWidth)]++;
        }
    }
    this->numberOfSamples++;
}

void Analysis::saveHisto(char outName[]){
    int i = 0;
    double dv = 0;
    double idealDen = 0;

    char pHisto_name[64];
    sprintf(pHisto_name, "pHisto_%d_", this->num);
    strcat(pHisto_name, outName);

    char nHisto_name[64];
    sprintf(nHisto_name, "nHisto_%d_", this->num);
    strcat(nHisto_name, outName);

    char histo_name[64];
    sprintf(histo_name, "histo_%d_", this->num);
    strcat(histo_name, outName);
    

    FILE *f = fopen(histo_name, "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    for(i = 0; i < bins; i++){
        fprintf(f, "%lf     %.15lf\n", (double)i * this->binWidth, (double)this->histo[i] / (Base::xL * Base::yL * this->binWidth * this->numberOfSamples));
    }
    fclose(f);


    FILE *pf = fopen(pHisto_name, "w");
    if(pf == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    for(i = 0; i < bins; i++){
        fprintf(pf, "%lf     %.15lf\n", (double)i * this->binWidth, (double)this->pHisto[i] / (Base::xL * Base::yL * this->binWidth * this->numberOfSamples));
    }
    fclose(pf);


    FILE *nf = fopen(nHisto_name, "w");
    if(nf == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    for(i = 0; i < bins; i++){
        fprintf(nf, "%lf     %.15lf\n", (double)i * this->binWidth, (double)this->nHisto[i] / (Base::xL * Base::yL * this->binWidth * this->numberOfSamples));
    }
    fclose(nf);
}


void Analysis::sample_rdf(Particle **particles, int *histo, double binWidth){
    int i = 0;
    int j = 0;
    double dist = 0;
    double xL2 = (Base::xL * Base::xL)/4;

    for(i = 0; i < Particle::numOfParticles; i++){
        j = i + 1;
        while(j < Particle::numOfParticles){
            dist = particles[i]->com_distance(particles[j]);
            if(dist < xL2){
                if(particles[i]->q < 0 && particles[j]->q < 0){
                    dist = sqrt(dist);
                    this->histo[(int)(dist/binWidth)] += 2;
                }
            }
            j++;
        }
    }
    numberOfSamples++;
}

void Analysis::save_rdf(int *histo, int bins, double binWidth){
    int i = 0;
    double dv = 0;
    double idealDen = 0;

    FILE *f = fopen("rdf_ewald_ewald.txt", "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    //Number of particles in ideal gas with same density
    for(i = 0; i < bins; i++){
        //dv = (double)4.0/3.0 * PI * (pow(i + 1, 3) - pow(i, 3)) * pow(binWidth, 3);
        dv = (double)4.0 * PI * pow(i * binWidth, 2) * binWidth;
        idealDen = dv * 1/2 * Particle::numOfParticles/(Base::xL * Base::yL * Base::zL);
        fprintf(f, "%lf     %lf\n", i * binWidth, histo[i]/(Particle::numOfParticles * 1/2 * numberOfSamples * idealDen));
    }
    fclose(f);
}
