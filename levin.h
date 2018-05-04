#ifndef Levin_H
#define Levin_H

#include "base.cpp"
#include <eigen3/Eigen/Dense>
#include "particle.h"
//#include <vector>
//#include <complex>
//#include <eigen3/Eigen/Dense>

class Levin{
    public:
        Levin();
        double get_energy();
        void initialize(Particle **particles);

    private:
        int kNumMax;
        double gamma;
        //std::vector< std::vector<double> > kVec;
        //std::complex<double> *rkVec;
        int kNum;
        int kMax;
        Eigen::MatrixXd kVectors;
        Eigen::ArrayXd kNorms;
        Eigen::ArrayXd f1;
        Eigen::ArrayXd f2;
        Eigen::ArrayXd f3;
        Eigen::ArrayXd f4;
        double get_polarization();
        
};

#endif