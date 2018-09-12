#include "sp3.h"

double split(double distance, double cutoff){
    double frac = distance/cutoff;
    double frac5 = frac * frac * frac * frac * frac;
    double frac6 = frac5 * frac;
    double frac7 = frac6 * frac;

    double s = 1 - 7/4 * frac + 21/4 * frac5 - 7 * frac6 + 5/2 * frac7;

    return s;
}

double get_self(Particle **particles, double cutoff){
    double energy = 0;
    for(int i = 0; i < Particle::numOfParticles; i++){
        energy += particles[i]->q * particles[i]->q;
    }
    energy *= -7/(8 * cutoff);

    return energy;
}

double sp3::get_short(Particle **particles, double cutoff){
    double energy = 0;
    int k = 0;
    double dist = 0;

    for(int i = 0; i < Particle::numOfParticles; i++){
        k = i + 1;
        while(k < Particle::numOfParticles){
            dist = sqrt(particles[i]->distance(particles[k]));
            energy += particles[i]->q * particles[k]->q * split(dist, cutoff) / dist;

        }
    }
    return energy;
}

double sp3::get_energy(Particle **particles){
    double eShort = 0;
    double eSelf = 0;
    double cutoff = 9;

    eShort = get_short(particles, cutoff);
    eSelf = get_self(particles, cutoff);
    return eShort + eSelf;
}
