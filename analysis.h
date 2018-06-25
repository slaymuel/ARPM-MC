#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "base.h"
#include "particle.h"

class Analysis: public Base{
    protected:
        int *histo;
        double binWidth;
        int bins;
        
        int num;
        static int numOfHisto;

    public:
        int *pHisto;
        int *nHisto;
        int numberOfSamples;
        Analysis(double binWidth, double dLength);
        void sampleHisto(Particle **particles, int d);
        void sample_rdf(Particle **particles, int *histo, double binWidth);
        void save_rdf(int *histo, int bins, double binWidth);
        void saveHisto(char outName[]);
};

#endif