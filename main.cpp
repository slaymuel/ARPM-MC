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
#include "levin.h"
#include <ctime>
#include "direct.h"
#include "valleau.h"
#include "ewald3D.h"
#include "ewald2D.h"
#include "hard_sphere.h"
#include "imagitron.h"
#include "img_rep.h"
//Initializers
//Ewald3D MC::ewald3D;

//Direct MC::direct;
//Levin MC::levin;

int Analysis::numOfHisto = 0;
double Base::xL = 114.862525361502; //172
double Base::yL = 114.862525361502; //45
double Base::zL = 45;
double Base::T = 298;
double Base::lB;
double Base::beta;
double Base::P = 1.0;
double Base::volume;
bool Base::d2 = false;
double Base::eCummulative = 0;
double Base::wall = 0;
int Base::acceptedMoves = 0;
int Base::totalMoves = 0;
std::vector<double> Base::volumes(1001);
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
    double dr = 0.5;
    //Simulation variables
    int totalMoves = 0;
    int acceptedMoves = 0;
    int numberOfSamples = 0;
    int iter = 1000000;
    double kb = KB;
    double T = Base::T;
    double beta = 1;
    double random1;
    double random2;
    int prevAccepted = 0;
    bool nanometers = false;
    int Digs = 14;
    char outName[40] = "output_ewald.gro\0";
    char outName_charges[50] = "charges_\0";
    std::string outputFile;
    //Analysis
    int bins = 500;
    double binWidth = (Base::xL/2)/bins;
    double dist = 0;
    int *histo;
    histo = (int*)malloc(bins * sizeof(int));
    for(int i = 0; i < bins; i++){
        histo[i] = 0;
    }
    //Simulation parameters
    MC mc;

    //Levin levin;
    //levin.initialize();

    //Command line parser
    po::options_description desc("Command line options:");
    desc.add_options()
        ("help", "show this message")
        ("np", po::value<int>(), "Number of particles")
        ("pnum", po::value<int>(), "Number of positive particles")
        ("nnum", po::value<int>(), "Number of negative particles")
        ("f", po::value<std::string>(), "Coordinates in xyz format")
        ("f_gro", po::value<std::string>(), "Coordinates in gro-format")
        ("f_jan", po::value<std::string>(), "Coordinates in Jan-format")
        ("f_arpm_jan", po::value<std::string>(), "Coordinates in Jan-ARPM format")
        ("density", po::value<double>(), "Specify density of the system.")
        ("wall", po::value<double>()->default_value(0), "Insert walls in the z-dimension. Argument increases the box length as zL = 2 * wall")
        ("box", po::value<std::vector<double> >()->multitoken(), "Box dimensions")
        ("rc", po::value<int>(), "Relative coordinates")
        ("T", po::value<double>(), "Temperature")
        ("o", po::value<std::string>(), "Output filename")
        ("nm", po::value<int>(), "If input file is in nanometers")
        ("iter", po::value<int>(), "Number of iteration (MC-moves)")
        ("dr", po::value<double>(), "Size of MC step")
        ("2d", po::bool_switch()->default_value(false), "System only extends in the x and y dimenstions")
        ("overlap", po::bool_switch()->default_value(false), "Remove Overlaps")
        ("electrons", po::value<int>()->default_value(0), "Number of electrons")
        ("imgrep", po::bool_switch()->default_value(false), "Image replicates in 3DEwald");
    
    po::variables_map vm;        
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
    if(vm["2d"].as<bool>()){
        Base::d2 = true;
    }
    if(vm.count("dr")){
        dr = vm["dr"].as<double>();
    }
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
    }
    if(vm.count("wall")){
        Base::wall = vm["wall"].as<double>();
        Base::zL += 2 * Base::wall;
    }
    if(vm.count("nm")){
        nanometers = true;
    }
    if(vm.count("f")){
        bool relative = false;
        if(vm.count("rc")){
            relative = true;
        }
        std::cout << "Opening " << vm["f"].as<std::string>() << std::endl;
        std::string filename = vm["f"].as<std::string>().c_str();
        particles = Particle::read_coordinates(filename, relative, nanometers);
    }
    if(vm.count("f_gro")){
        std::cout << "Opening " << vm["f_gro"].as<std::string>() << std::endl;
        std::string filename = vm["f_gro"].as<std::string>().c_str();
        particles = Particle::read_coordinates_gro(filename);
    }
    if(vm.count("o")){
        outputFile = vm["o"].as<std::string>();
        strcpy(outName, outputFile.c_str());
    }
    if(vm.count("f_jan")){
        particles = Particle::read_jan("coordp_dense.dms", "coordn_dense.dms");
    }
    if(vm.count("f_arpm_jan")){
        particles = Particle::read_arpm_jan("coord");
    }
    if(vm.count("electrons")){
        Particle::numOfElectrons = vm["electrons"].as<int>();
    }
    if(vm.count("np") && vm["imgrep"].as<bool>()){
        numOfParticles = vm["np"].as<int>() * 5;
        if(numOfParticles % 2 != 0){
            printf("\033[31mIgnoring bad choice of number of particles.\033[30m\n");
            numOfParticles += 5;
        }
        printf("%d particles will be created\n", numOfParticles);
        particles = Particle::create_particles(numOfParticles/10, numOfParticles/10, Particle::numOfElectrons);
        energy::imgrep::set_positions(particles);
        //mc.equilibrate(particles);
        density = (double)numOfParticles/(Base::xL * Base::yL * (Base::zL - 2 * Base::wall)) * pow(diameter, 3);
    }
    else if(vm.count("np")){
        numOfParticles = vm["np"].as<int>();
        if(numOfParticles % 2 != 0){
            printf("\033[31mIgnoring bad choice of number of particles.\033[30m\n");
            numOfParticles += 1;
        }
        printf("%d particles will be created\n", numOfParticles);
        particles = Particle::create_particles(numOfParticles/2, numOfParticles/2, Particle::numOfElectrons);
        //mc.equilibrate(particles);
        density = (double)numOfParticles/(Base::xL * Base::yL * (Base::zL - 2 * Base::wall)) * pow(diameter, 3);
    }
    if(vm.count("T")){
        Base::T = vm["T"].as<double>();
        T = Base::T;
    }
    if(vm.count("wall") && vm.count("density")){
        density = vm["density"].as<double>();
        Base::wall = vm["wall"].as<double>();
        numOfParticles = density/pow(diameter, 3) * (Base::xL * Base::yL * (Base::zL - 2 * Base::wall));
        if(numOfParticles % 2 != 0){
            printf("Uneven number of particles, adding one....\n");
            numOfParticles += 1;
        }
        density = (double)numOfParticles/(Base::xL * Base::yL * (Base::zL - 2 * Base::wall)) * pow(diameter, 3);
        //numOfParticles = vm["density"].as<double>()/pow(diameter, 3)*(Base::xL * Base::yL * Base::zL);
        printf("%d particles will be created\n", numOfParticles);
        particles = Particle::create_particles(numOfParticles/2, numOfParticles/2, Particle::numOfElectrons);
        //mc.equilibrate(particles);
    }
    else{
        if(vm.count("density")){
            numOfParticles = vm["density"].as<double>()/pow(diameter, 3)*(Base::xL * Base::yL * Base::zL);
            if(numOfParticles % 2 != 0){
                numOfParticles += 1;
            }
            printf("%d particles will be created\n", numOfParticles);
            particles = Particle::create_particles(numOfParticles/2, numOfParticles/2, Particle::numOfElectrons);
        }
    }
    if(vm.count("iter")){
        iter= vm["iter"].as<int>();
    }
    if(vm.count("pnum")){
        particles = Particle::create_particles(vm["nnum"].as<int>(), vm["pnum"].as<int>(), vm["electrons"].as<int>());
        density = (double)Particle::numOfParticles/(Base::xL * Base::yL * Base::zL) * pow(diameter, 3);
    }
    printf("\033[34mBox dimensions are x: %lf y: %lf z: %lf\033[30m\n", Base::xL, Base::yL, Base::zL);
    printf("num of particles including images: %d\n", Particle::numOfParticles);

    Particle::distances = (double**) malloc((Particle::numOfParticles + Particle::numOfElectrons) * sizeof(double*));
    for(int i = 0; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
        Particle::distances[i] = (double*) malloc((Particle::numOfParticles + Particle::numOfElectrons) * sizeof(double*));
    }
    
    //Seed
    srand(time(NULL));

    if(density == 0){
        density = (double)numOfParticles/(Base::xL * Base::yL * (Base::zL - 2 * Base::wall)) * pow(diameter, 3);
    }
    printf("\033[34mDensity is: %lf\033[30m\n", density);

    Particle::update_distances(particles);

    if(vm["overlap"].as<bool>()){
        printf("Removing overlaps...\n");
        mc.equilibrate(particles);
        //Particle::place_particles(particles);
    }

    char name[] = "output_equilibrate_ewald.gro";
    char name_charges[] = "output_equilibrate_charges_ewald.gro";

    Particle::write_coordinates(name, particles);
    Particle::write_charge_coordinates(name_charges, particles);

    Particle::update_distances(particles);
    //exit(1);
    //Analysis *xHist = new Analysis(0.05, Base::xL);
    //Analysis *yHist = new Analysis(0.05, Base::yL);
    //Analysis *zHist = new Analysis(0.05, Base::zL);
    //Analysis *rdf = new Analysis(0.1, Base::zL);
    Base::set_lB();
    Base::set_beta();
    Base::volume = Base::xL * Base::yL * Base::zL;
    Base::volumes.push_back(Base::volume);

    energy::ewald2D::set_alpha();
    energy::ewald2D::initialize();
    energy::ewald3D::set_alpha();
    energy::ewald3D::initialize(particles);

    energy::levin::initialize(particles);
    energy::valleau::initialize();

    // overlaps = Particle::get_overlaps(particles);
    // if(overlaps > 0){
    //     printf("System contains overlaps!\n");
    //     exit(0);
    // }

    //Update cumulative energy
    int energyOut = 0;

    // ///////////////////////////////         Main MC-loop          ////////////////////////////////////////
    //Ewald3D ewald;

    strcat(outName_charges, outputFile.c_str());
    char volOut[40] = "volumes_\0";
    strcat(volOut, outputFile.c_str());
    FILE *f = fopen(volOut, "w");
    fprintf(f, "");
    fclose(f);

    std::string energyFunction = "ewald";

    if(energyFunction == "valleau"){
        //MC::run(&energy::hs::get_energy, &energy::hs::get_particle_energy, particles, dr, iter, false);

        MC::run(&energy::valleau::get_energy, &energy::valleau::get_particle_energy, particles, 0.1, 1000000, false, outputFile);
        energy::valleau::update_potential();
        
        MC::run(&energy::valleau::get_energy, &energy::valleau::get_particle_energy, particles, dr, 1000000, false, outputFile);
        energy::valleau::update_potential();
        
        //MC::run(&energy::valleau::get_energy, &energy::valleau::get_particle_energy, particles, 15, 1000000, false);
        //energy::valleau::update_potential();

        MC::run(&energy::valleau::get_energy, &energy::valleau::get_particle_energy, particles, dr, iter, true, outputFile);
        
        //MC::run(&energy::levin::get_energy, &energy::levin::get_particle_energy, particles, dr, iter, true, outputFile);
    }

    if(energyFunction == "ewald"){
        MC::run(&energy::ewald3D::get_energy, &energy::ewald3D::get_particle_energy, particles, dr, iter, true, outputFile);
    }

    if(energyFunction == "ewald2d"){
        MC::run(&energy::ewald2D::get_energy, &energy::ewald2D::get_particle_energy, particles, dr, iter, true, outputFile);
    }

    if(energyFunction == "test"){
        MC::run(&energy::valleau::get_energy, &energy::valleau::get_particle_energy, particles, dr, 0, true, outputFile);
        MC::run(&energy::levin::get_energy, &energy::levin::get_particle_energy, particles, dr, 0, true, outputFile);
    }

    if(energyFunction == "levin"){
        MC::run(&energy::levin::get_energy, &energy::levin::get_particle_energy, particles, dr, iter, true, outputFile);
    }

    if(energyFunction == "electron"){
        energy::imagitron::initialize();
        MC::run(&energy::imagitron::get_energy, &energy::imagitron::get_particle_energy, particles, dr, iter, true, outputFile);
    }
    if(energyFunction == "direct"){
        MC::run(&energy::direct::get_energy, &energy::direct::get_particle_energy, particles, dr, iter, true, outputFile);
    }
/*
    Base::eCummulative = MC::ewald3D.get_energy(particles);
    //Base::eCummulative = energy::valleau::get_energy(particles);
    double avgTime = 0;
    printf("\nRunning main MC-loop at temperature: %lf, Bjerrum length is %lf\n\n", T, Base::lB);

    for(int i = 0; i < iter; i++){
        if(i % 100 == 0 && i >= 1){
            //rdf->sample_rdf(particles, histo, binWidth);
            xHist->sampleHisto(particles, 0);
            yHist->sampleHisto(particles, 1);
            zHist->sampleHisto(particles, 2);
        }

        random1 = ran2::get_random();
        //random2 = ran2::get_random();
        //if(random2 <= 0.3){
        //     if(mc.charge_rot_move(particles)){
        //        prevAccepted++;
        //     }
        //}
        //else{
            if(random1 <= 0.1){
                if(mc.trans_move(particles, Base::xL)){
                    prevAccepted++; 
                }
            }

            else{
                if(mc.trans_move(particles, 0.2)){
                    prevAccepted++; 
                }    
            }
        //}


        Base::totalMoves++;
        if(i % 10000 == 0 && i != 0) {
            //energy = energy::valleau::get_energy(particles);
            energy = MC::ewald3D.get_energy(particles);
            //energy = MC::direct.get_energy(particles);
            energyOut++;
            //Particle::write_coordinates(outName , particles);
            printf("Iteration: %d\n", i);
            printf("Energy: %lf\n", energy);
            //printf("Error: ");
            //printf("%lf\n", fabs(energy - Base::eCummulative)/fabs(Base::eCummulative));
            printf("Acceptance ratio: %lf\n", (double) Base::acceptedMoves / Base::totalMoves);
            printf("Acceptance ratio for the last 10000 steps: %lf\n\n", (double) prevAccepted / 10000.0);
            if(fabs(energy - Base::eCummulative) / fabs(energy) > pow(10, -12)) {
                printf("Error is too large!\n");
                printf("Error: %.20lf\n", fabs(energy - Base::eCummulative) / fabs(energy));
                //exit(1);
            }
            fflush(stdout);
            prevAccepted = 0;
            //printf("Time: %lf\n", omp_get_wtime() - stime);
        }
        //avgTime += omp_get_wtime() - stime;
        //printf("Real time: %lf\n", avgTime/(i + 1);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////
*/
    printf("Accepted moves: %d\n", Base::acceptedMoves);
    printf("Rejected moves: %d\n", Base::totalMoves - Base::acceptedMoves);

    //rdf->save_rdf(histo, bins, binWidth);
    //xHist->saveHisto(outName);
    //yHist->saveHisto(outName);
    //zHist->saveHisto(outName);

    //Write coordinates to file
    printf("Saving output coordinates to: %s\n", outName);
    Particle::write_coordinates(outName, particles);
    Particle::write_charge_coordinates(outName_charges, particles);

    

    //Clean up allocated memory
    printf("Cleaning up...\n");
    for(int i = 0; i < Particle::numOfParticles; i++){
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
