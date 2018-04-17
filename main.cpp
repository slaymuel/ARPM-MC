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
double Base::xL = 114.862525361502;
double Base::yL = 114.862525361502;
double Base::zL = 45;
double Base::T = 1000;
double Base::lB;
double Base::eCummulative = 0;
double Base::wall = 0;
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
    double density = 0;
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
    //int bins = 100;
    //double binWidth = (Base::xL/2)/bins;
    double dist = 0;
    // int *histo;
    // histo = (int*)malloc(bins * sizeof(int));

    //Simulation parameters
    MC mc;

    //Command line parser
    po::options_description desc("Command line options:");
    desc.add_options()
        ("help", "show this message")
        ("np", po::value<int>(), "Number of particles")
        ("f", po::value<std::string>(), "Coordinates in xyz format")
        ("f_jan", po::value<std::string>(), "Coordinates in Jan-format")
        ("density", po::value<double>(), "Specify density of the system.")
        ("wall", po::value<double>(), "Insert walls in the z-dimension.")
        ("box", po::value<std::vector<double> >()->multitoken(), "Box dimensions")
        ("rc", po::value<int>(), "Relative coordinates")
        ("T", po::value<double>(), "Temperature");
    
    po::variables_map vm;        
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(1);
    }
    if(vm.count("box")){
        std::vector<double> box(3);
        box = vm["box"].as<std::vector<double> >();
        Base::xL = box[0];
        Base::yL = box[1];
        Base::zL = box[2];
        printf("Set box dimensions to: %lf %lf %lf\n", Base::xL, Base::yL, Base::zL);
    }
    if(vm.count("wall")){
        Base::wall = vm["wall"].as<double>();
    }

    if(vm.count("f")){
        bool relative = false;
        if(vm.count("rc")){
            relative = true;
        }
        std::cout << "Opening " << vm["f"].as<std::string>() << std::endl;
        std::string filename = vm["f"].as<std::string>().c_str();
        particles = Particle::read_coordinates(filename, relative);
    }
    if(vm.count("f_jan")){
        particles = Particle::read_jan("coordp-2.dms", "coordn-2.dms");
    }
    if(vm.count("np")){
        numOfParticles = vm["np"].as<int>();
        printf("%d particles will be created\n", numOfParticles);
        particles = Particle::create_particles(numOfParticles);
        mc.equilibrate(particles);
        density = (double)numOfParticles/(Base::xL * Base::yL * Base::zL) * pow(diameter, 3);
    }
    if(vm.count("T")){
        Base::T = vm["T"].as<double>();
        T = Base::T;
    }
    if(vm.count("wall") && vm.count("density")){
        density = vm["density"].as<double>();
        Base::wall = vm["wall"].as<double>();
        numOfParticles = density/pow(diameter, 3) * (Base::xL * Base::yL * (Base::zL - 2 * Base::wall));
        density = (double)numOfParticles/(Base::xL * Base::yL * (Base::zL - 2 * Base::wall)) * pow(diameter, 3);
        //numOfParticles = vm["density"].as<double>()/pow(diameter, 3)*(Base::xL * Base::yL * Base::zL);
        printf("%d particles will be created\n", numOfParticles);
        printf("Density is: %lf\n", density);
        particles = Particle::create_particles(numOfParticles);
        mc.equilibrate(particles);
    }
    else{
        if(vm.count("density")){
            numOfParticles = vm["density"].as<double>()/pow(diameter, 3)*(Base::xL * Base::yL * Base::zL);
            printf("%d particles will be created\n", numOfParticles);
            particles = Particle::create_particles(numOfParticles);
            mc.equilibrate(particles);
            density = (double)numOfParticles/(Base::xL * Base::yL * Base::zL) * pow(diameter, 3);
        }
    }
    Analysis *xHist = new Analysis(0.1, Base::xL);
    Analysis *yHist = new Analysis(0.1, Base::yL);
    Analysis *zHist = new Analysis(0.1, Base::zL);

    MC::ewald3D.initialize(particles);
    MC::ewald2D.initialize();
    //Seed
    srand(time(NULL));

    // overlaps = Particle::get_overlaps(particles);
    // if(overlaps > 0){
    //     printf("System contains overlaps!\n");
    //     exit(0);
    // }
    printf("Box dimensions are x: %lf y: %lf z: %lf\n", Base::xL, Base::yL, Base::zL);
    char name[] = "output_equilibrate.gro";
    Particle::write_coordinates(name , particles);

    //Update cumulative energy
    //Base::eCummulative = mc.get_energy(particles);

    // ///////////////////////////////         Main MC-loop          ////////////////////////////////////////
    printf("\nRunning main MC-loop at temperature: %lf\n\n", T);
    for(int i = 0; i < 500000; i++){
        if(i % 1000 == 0 && i > 100){
            //sampleRDF(particles, zHisto, binWidth);
            //printf("Sampling...\n");
            xHist->sampleHisto(particles, 0);
            yHist->sampleHisto(particles, 1);
            zHist->sampleHisto(particles, 2);
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
            //printf("Acceptance ratio for the last 100000 steps: %lf\n\n", (double)prevAccepted/100000.0);
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
    //xHist->saveHisto();
    //yHist->saveHisto();
    //zHist->saveHisto();

    //Write coordinates to file
    char outName[] = "output_ewald.gro";
    Particle::write_coordinates(outName , particles);
    // FILE *f = fopen("output_ewald.gro", "w");
    // if(f == NULL){
    //     printf("Can't open file!\n");
    //     exit(1);
    // }

    // fprintf(f, "%d\n", numOfParticles);
    // fprintf(f, "\n");
    // for(int i = 0; i < Particle::numOfParticles; i++){
    //     fprintf(f, "%s     %lf    %lf     %lf\n", particles[i]->name, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    // }
    // fclose(f);

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
    //free(histo);
   return 0;
}