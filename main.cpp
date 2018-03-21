//Create particles
//Check overlap
//Move particle, accept move if no overlap
#include "particle.h"
#include "constants.h"
#include "analysis.h"
#include "ran2_lib.cpp"
#include "mc.h"
#include "boost/program_options.hpp"
#include <iterator>

//Initializers
Ewald3D MC::ewald3D;
Ewald2D MC::ewald2D;
Direct MC::direct;

int Particle::numOfParticles = 0;
int Analysis::numOfHisto = 0;
double Base::xL= 30;
double Base::yL= 30;
double Base::zL= 10;
double Base::T = 500;
double Base::lB;
double Base::eCummulative = 0;
int Base::acceptedMoves = 0;
int Base::totalMoves = 0;

namespace po = boost::program_options;


int main(int argc, char *argv[])
{
    printf("\n");

    // if( argc == 2 ) {
    //     printf("input: %s\n", argv[1]);
    //   if(!strcmp(argv[1], "rdf")){
    //       //readCoo();
    //       exit(1);
    //   }
    // }
    
    //Particle variables
    Particle **particles;
    int numOfParticles;//6728;
    double diameter = 5;
    double diameter2 = pow(diameter, 2);
    int overlaps = 0;
    double density = (double)numOfParticles/(Base::xL * Base::yL * Base::zL) * pow(diameter, 3);
    int overlap = 1;
    double energy = 0;
    
    //Simulation variables
    int totalMoves = 0;
    int acceptedMoves = 0;
    int numberOfSamples = 0;
    double kb = KB;
    double T = Base::T;
    double beta = 1;
    int prevAccepted = 0;
    Base::set_lB();

    int Digs = 14;

    //Analysis
    int bins = 100;
    double binWidth = (Base::xL/2)/bins;
    double dist = 0;
    int *histo;
    histo = (int*)malloc(bins * sizeof(int));
    Analysis *xHist = new Analysis(0.1, Base::xL);
    Analysis *yHist = new Analysis(0.1, Base::yL);
    Analysis *zHist = new Analysis(0.1, Base::zL);

    //Simulation parameters
    MC mc;
    MC::ewald3D.initialize();
    //MC::ewald2D.initialize();

    //Command line parser
    po::options_description desc("Command line options:");
    desc.add_options()
        ("help", "show this message")
        ("compression", po::value<int>(), "set compression level")
        ("np", po::value<int>(), "Number of particles")
        ("f", po::value<std::string>(), "Coordinates in xyz format");
    
    po::variables_map vm;        
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }
    if(vm.count("f")){
        std::cout << "Opening " << vm["f"].as<std::string>() << std::endl;
        std::string filename = vm["f"].as<std::string>().c_str();
        particles = Particle::read_coordinates(filename);
        return 0;
    }
    if(vm.count("np")){
        numOfParticles = vm["np"].as<int>();
        printf("%d particles will be created\n", numOfParticles);
        particles = Particle::create_particles(numOfParticles);
        mc.equilibrate(particles);
        return 0;
    }
    
    //Seed
    srand(time(NULL));

    printf("Density is: %lf\n", density);

    overlaps = Particle::get_overlaps(particles);
    if(overlaps > 0){
        printf("System contains overlaps!\n");
        exit(0);
    }

    // char name[] = "output_equilibrate_newest.xyz";
    // Particle::write_coordinates(name , particles);

    //Update cumulative energy
    Base::eCummulative = mc.get_energy(particles);

    // ///////////////////////////////         Main MC-loop          ////////////////////////////////////////
    printf("\nRunning main MC-loop at temperature: %lf\n", T);
    for(int i = 0; i < 5000; i++){
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
        
        if(i % 100 == 0 && i != 0){
            energy = mc.get_energy(particles);
            printf("Iteration: %d\n", i);
            printf("Energy: %lf\n", energy);
            //printf("Error: ");
            //printf("%.*e\n", Digs, fabs(energy - Base::eCummulative)/fabs(Base::eCummulative));
            printf("Acceptance ratio: %lf\n", (double)Base::acceptedMoves/Base::totalMoves);
            printf("Acceptance ratio for the last 100000 steps: %lf\n\n", (double)prevAccepted/100000.0);
            //if(fabs(energy - Base::eCummulative)/fabs(energy) > pow(10, -12)){
            //    printf("Error is too large!\n");
            //    exit(1);
            //}
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
    FILE *f = fopen("output_ewald.xyz", "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }

    fprintf(f, "%d\n", numOfParticles);
    fprintf(f, "\n");
    for(int i = 0; i < Particle::numOfParticles; i++){
        fprintf(f, "%s     %lf    %lf     %lf\n", particles[i]->name, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    }
    fclose(f);

    //Clean up allocated memory
    printf("Cleaning up...\n");
    for(int i = 0; i < Particle::numOfParticles; i++){
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