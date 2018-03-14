#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "base.cpp"
#include "particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

class Analysis: public Base{
    protected:
        int *histo;
        int *pHisto;
        int *nHisto;
        double binWidth;
        int bins;
        int numberOfSamples;
        int num;
        static int numOfHisto;

    public:
        Analysis(double binWidth);
        void sampleHisto(Particle **particles, int d);
        void saveHisto();
};

#endif