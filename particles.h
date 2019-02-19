#pragma once
#include "particle.h"

class Particles{
    public:
    std::vector<Particle> particles;
    int numOfParticles, numOfCations, numOfAnions, numOfElectrons;
    std::vector< std::vector<double> > distances;
    double oldEnergy;

    Particle& operator[] (int x) {
        return particles[x];
    }
    



    Particles(){
        numOfParticles = 0;
        numOfAnions = 0;
        numOfCations = 0;
        numOfElectrons = 0;
    }






    void initialize(){
        numOfParticles = particles.size();
        distances.resize(numOfParticles);
        for(int i = 0; i < numOfParticles; i++){
            distances[i].resize(numOfParticles);

            if(particles[i].q < 0){
                numOfAnions++;
            }
            else{
                numOfCations++;
            }
        }
        printf("Cations: %i, Anions: %i\n", numOfCations, numOfAnions);
    }





    bool add(){
        Particle temp = new Particle();
        temp.d = 5.0;
        temp.com[0] = (double) rand()/RAND_MAX * (-Base::xL) + Base::xL / 2.0;
        temp.com[1] = (double) rand()/RAND_MAX * (-Base::yL) + Base::yL / 2.0;
        //particles[i]->com[2] = (double) rand()/RAND_MAX * (Base::zL - particles[i]->d - 2 * Base::wall) + particles[i]->d/2.0 + Base::wall;
        temp.com[2] = (double) rand()/RAND_MAX * -2.0 * (Base::zLBox / 2.0 - temp.d / 2.0) + Base::zLBox / 2.0 - temp.d / 2.0;
        temp.index = -1;
        if(hard_sphere(temp)){
            particles.push_back(temp);
            return true;
        }

        else{
            return false;
        }

    }






    void remove(int ind){
        numOfParticles--;
        if(particles[ind].q > 0){
            numOfCations--;
        }
        else{
            numOfAnions--;
        }
        particles.erase(particles.begin() + ind);
        for(auto &row : distances){
            row.erase(std::next(row.begin(), ind));
        }
        distances.erase(std::next(distances.begin(), ind));
    }





    void read_coordinates_gro(std::string name){
        int i = 0;
        int j = 0;
        double x, y, z;
        int c, ind;
        double nano;
        std::string molecule;
        std::string atom;
        std::string line;
        std::ifstream infile(name);

        while (std::getline(infile, line))
        {
            if(i == 1){
                std::istringstream iss(line);
                if (!(iss >> c)) {
                    printf("The second line in the input file should be the total number of atoms!\n");
                    exit(1); 
                } // error
            }

            if(i > 1){
                std::istringstream iss(line);
                if (!(iss >> molecule >> atom >> ind >> x >> y >> z)) {
                    printf("Done reading input file\n");
                    break;
                    //exit(1);
                }
                //particles[j] = new Particle();
                particles.push_back(new Particle());

                particles[j].d = 5.0;

                particles[j].index = j;
                particles[j].com[0] = x * 10.0 - Base::xL / 2.0;
                particles[j].com[1] = y * 10.0 - Base::yL / 2.0;
                particles[j].com[2] = z * 10.0 - Base::zLBox / 2.0;    

                //Get random charge displacement vector
                particles[j].chargeDisp[0] = (double) rand()/RAND_MAX * 2.0 - 1.0;
                particles[j].chargeDisp[1] = (double) rand()/RAND_MAX * 2.0 - 1.0;
                particles[j].chargeDisp[2] = (double) rand()/RAND_MAX * 2.0 - 1.0;

                if(atom == "Cl"){
                    strcpy(particles[j].name, "Cl\0");
                    particles[j].b = 0;
                    particles[j].q = -1.0;
                }

                else if(atom == "Na"){
                    strcpy(particles[j].name, "Na\0");
                    particles[j].b = 0;
                    particles[j].q = 1.0;            
                }

                else{
                    printf("Atom is not Na or Cl!\n");
                    exit(1);
                }

                particles[j].chargeDisp = particles[j].b * particles[j].chargeDisp.normalized();
                //Calculate position of the charge
                particles[j].pos = particles[j].com + particles[j].chargeDisp;
                j++;
            }
            i++;
        }
        if(i == 0){
            printf("Could not open file: %s\n", name.c_str());
        }
        if(c != j){
            printf("Did not read all particles from file, is the second line really the number of particles?\n");
            exit(1);
        }
        printf("%d particles read from file.\n", j);
    }






    void create_particles(int nNum, int pNum, int eNum){

        int i = 0;
        int num = nNum + pNum;
        double norm = 0;

        for(i = 0; i < num; i++){
            particles.push_back(new Particle());
            particles[i].index = i;
            particles[i].d = 5.0;    //Diameter of particles

            if(i < nNum){
                particles[i].q = -1.0;
                particles[i].b = 0.0; //Length of charge displacement vector
                strcpy(particles[i].name, "Cl\0");
            }

            else{
                particles[i].q = 1.0;
                particles[i].b = 0.0; //Length of charge displacement vector
                strcpy(particles[i].name, "Na\0");
            }
            

            //Get random center of mass coordinates
            particles[i].com[0] = (double) rand()/RAND_MAX * (-Base::xL) + Base::xL / 2.0;
            particles[i].com[1] = (double) rand()/RAND_MAX * (-Base::yL) + Base::yL / 2.0;
            //particles[i]->com[2] = (double) rand()/RAND_MAX * (Base::zL - particles[i]->d - 2 * Base::wall) + particles[i]->d/2.0 + Base::wall;
            particles[i].com[2] = (double) rand()/RAND_MAX * -2.0 * (Base::zLBox / 2.0 - particles[i].d / 2.0) + Base::zLBox / 2.0 - particles[i].d / 2.0;

            //Get random charge displacement vector
            particles[i].chargeDisp[0] = (double) rand()/RAND_MAX * 2 - 1;
            particles[i].chargeDisp[1] = (double) rand()/RAND_MAX * 2 - 1;
            particles[i].chargeDisp[2] = (double) rand()/RAND_MAX * 2 - 1;

            particles[i].chargeDisp = particles[i].b * particles[i].chargeDisp.normalized();
            //Calculate position of the charge
            particles[i].pos = particles[i].com + particles[i].chargeDisp;
            
            if(particles[i].com[2] < -Base::zLBox / 2.0 + particles[i].d/2 || particles[i].com[2] > Base::zLBox / 2.0 - particles[i].d/2){
                printf("%lf %lf %lf Particle in forbidden area...\n", particles[i].com[0], particles[i].com[1], particles[i].com[2]);
                exit(1);
            }
        }
        printf("\033[34mCreated %d particles.\033[30m\n", num);
    }








    void create_electrons(int num){
        int i, j = 0;
        i = numOfParticles;

        for(j = 0; j < num; j++){

            particles.push_back(new Particle(true));
            particles[i].index = i;
            particles[i].d = 0;    //Diameter of particles
            particles[i].b = 0; //Length of charge displacement vector

            //Get random center of mass coordinates
            particles[i].com[0] = (double) rand()/RAND_MAX * Base::xL;
            particles[i].com[1] = (double) rand()/RAND_MAX * Base::yL;
            if(j < num / 2){
                particles[i].com[2] = Base::zL;//(double) rand()/RAND_MAX * (Base::zL - 2 * Base::wall) + Base::wall + Base::zL;
            }
            else{
                particles[i].com[2] = 0;//(double) rand()/RAND_MAX * (Base::zL - 2 * Base::wall) + Base::wall - Base::zL;
            }

            //Get random charge displacement vector
            particles[i].chargeDisp << 0, 0, 0;
            particles[i].pos = particles[i].com;
            
            if(particles[i].com[2] > 2 * Base::zL){
                printf("%lf %lf %lf Electron in forbidden area...\n", particles[i].com[0], particles[i].com[1], particles[i].com[2]);
                exit(1);
            }

            particles[i].q = -1.0;
            strcpy(particles[i].name, "e\0");
            i++;
        }

        printf("\033[34mCreated %d electrons.\033[30m\n", j);
    }








    void read_jan(std::string pName, std::string nName){
        int i = 0;
        int j = 0;
        double x, y, z;
        int c;
        std::string line;
        std::ifstream infileP(pName);

        //Read positive particles
        while (std::getline(infileP, line)){
            if(i < 1){
                std::istringstream iss(line);
                if (!(iss >> c)) {
                    printf("Error reading file...\n");
                    exit(1); 
                } // error
            }
            if(i >= 1){
                
                std::istringstream iss(line);
                if (!(iss >> x >> y >> z)) {
                    break; 
                } // error
                particles.push_back(new Particle());
                //particles[j]->pos = (double*) malloc(3 * sizeof(double));

                particles[j].pos[0] = x;
                particles[j].pos[1] = y;
                particles[j].pos[2] = z - Base::zLBox / 2.0 - 2.5;
                particles[j].com = particles[j].pos;

                particles[j].d = 5;
                particles[j].b = 0;
                particles[j].index = j;
                
                strcpy(particles[j].name, "Na\0");
                particles[j].q = 1.0;

                j++;
            }
            i++;
        }

        //Read negative particles
        i = 0;
        std::ifstream infileN(nName);
        while (std::getline(infileN, line)){
            if(i >= 0){
                std::istringstream iss(line);
                if (!(iss >> x >> y >> z)) {
                    break; 
                } // error
                particles.push_back(new Particle());
                //particles[j]->pos = (double*) malloc(3 * sizeof(double));

                particles[j].pos[0] = x;
                particles[j].pos[1] = y;
                particles[j].pos[2] = z - Base::zLBox / 2.0 - 2.5;
                particles[j].com = particles[j].pos;

                particles[j].d = 5;
                particles[j].b = 0;
                particles[j].index = j;
                
                strcpy(particles[j].name, "Cl\0");
                particles[j].q = -1.0;

                j++;
            }
            i++;
        }

        printf("%d particles read from file.\n", j);
    }








    void write_coordinates(char name[40]){
        int i = 0;
        FILE *f = fopen(name, "w");
        if(f == NULL){
            printf("Can't open file!\n");
            exit(1);
        }
        fprintf(f, "Generated by Slaymulator.\n");
        fprintf(f, "%d\n", numOfParticles + numOfElectrons);
        for(i = 0; i < numOfParticles + numOfElectrons; i++){
            fprintf(f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "ion", particles[i].name, i, (particles[i].com[0] + Base::xL / 2.0)/10.0, (particles[i].com[1] + Base::yL / 2.0)/10.0, (particles[i].com[2] + Base::zLBox / 2.0)/10.0);
        }
        fprintf(f, "%lf    %lf     %lf\n", Base::xL/10.0, Base::yL/10.0, Base::zL/10.0);
        fclose(f);
    }







    //void write_charge_coordinates(char name[], Particle **particles){
    void write_charge_coordinates(char name[]){
        int i = 0;
        FILE *f = fopen(name, "w");
        if(f == NULL){
            printf("Can't open file!\n");
            exit(1);
        }
        fprintf(f, "Generated by Slaymulator.\n");
        fprintf(f, "%d\n", numOfParticles);
        for(i = 0; i < numOfParticles; i++){
            fprintf(f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "cha", "H", i, (particles[i].pos[0] + Base::xL / 2.0) / 10.0, (particles[i].pos[1] + Base::yL / 2.0) / 10.0, (particles[i].pos[2] + Base::zLBox / 2.0) / 10.0);
        }
        fprintf(f, "%lf    %lf     %lf\n", Base::xL/10.0, Base::yL/10.0, Base::zL/10.0);
        fclose(f);
    }







    bool hard_sphere(Particle &p){
        int i = 0;
        int j = 0;

        for(i = 0; i < numOfParticles; i++){
            if(Base::wall > 0 || Base::d2){
                if(i != p.index){
                    if(p.com_distance_xy(particles[i]) < (p.d + particles[i].d) / 2.0){
                        return false;
                    }
                }
            }
            else{
                if(i != p.index){
                    if(p.com_distance(particles[i]) < (p.d + particles[i].d) / 2.0){
                        return false;
                    }
                }
            }
        }
        return true;
    }








    int get_overlaps(){
        int overlaps = 0;
        int i = 0;
        int j = 0;

        for(i = 0; i < numOfParticles; i++){
            j = i + 1;
            while(j < numOfParticles){
                if(i != j){
                    if(particles[i].com_distance_xy(particles[j]) < (particles[i].d/2 + particles[j].d/2)){
                        overlaps++;
                    }
                }
                j++;
            }
        }
        return overlaps;
    }









    void update_distances(Particle &p){
        if(Base::wall > 0 || Base::d2){
            for(int i = p.index + 1; i < numOfParticles + numOfElectrons; i++){
                distances[p.index][i] = p.distance_xy(particles[i]);
            }
            for(int i = 0; i < p.index; i++){
                distances[i][p.index] = p.distance_xy(particles[i]);
            }
        }
        else{
            for(int i = p.index + 1; i < numOfParticles + numOfElectrons; i++){
                distances[p.index][i] = p.distance(particles[i]);
            }
            for(int i = 0; i < p.index; i++){
                distances[i][p.index] = p.distance(particles[i]);
            }
        }
    }








    void update_distances(){
        int k = 0;
        if(Base::wall > 0 || Base::d2){
            printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING ONLY USING PBC IN TWO DIMENSIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
            for(int i = 0; i < numOfParticles + numOfElectrons; i++){
                //printf("i = %d\n", i);
                k = i + 1;
                while(k < numOfParticles + numOfElectrons){
                    //printf("k = %d\n", k);
                    distances[i][k] = particles[i].distance_xy(particles[k]);
                    k++;
                }
            }
        }
        else{
            for(int i = 0; i < numOfParticles + numOfElectrons; i++){
                k = i + 1;
                while(k < numOfParticles + numOfElectrons){
                    distances[i][k] = particles[i].distance(particles[k]);
                    k++;
                }
            }
        }
    }

};