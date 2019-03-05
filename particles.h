#pragma once
#include "particle.h"

class Particles{
    public:
    std::vector<Particle> particles;
    std::vector<int> cations;
    std::vector<int> anions;
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







    /*bool add(double z){
        Particle* temp = new Particle();
        temp->d = 5.0;
        temp->com[0] = (double) rand()/RAND_MAX * (-Base::xL) + Base::xL / 2.0;
        temp->com[1] = (double) rand()/RAND_MAX * (-Base::yL) + Base::yL / 2.0;
        temp->com[2] = z;
        temp->pos = temp->com;
        temp->q = 1.0;
        temp->index = numOfParticles;

        if(hard_sphere(*temp)){
            this->distances.resize(numOfParticles + 1);

            for(int i = 0; i < numOfParticles; i++){
                this->distances[i].push_back( particles[i].distance_xy(*temp) );
            }

            particles.push_back(*temp);
            numOfParticles++;
            numOfCations++;

            delete temp;
            return true;
        }

        else{
            delete temp;
            return false;
        }

    }*/

    void add(std::vector<double> pos){
        Particle* temp = new Particle();
        temp->d = 5.0;
        temp->com[0] = pos[0];
        temp->com[1] = pos[1];
        //particles[i]->com[2] = (double) rand()/RAND_MAX * (Base::zL - particles[i]->d - 2 * Base::wall) + particles[i]->d/2.0 + Base::wall;
        temp->com[2] = pos[2]; // 25.0 - Base::zLBox / 2.0;//
        temp->pos = temp->com;
        temp->index = numOfParticles;
        //temp->index = numOfCations;
        std::vector<double> tempVec(numOfParticles);

        numOfParticles++;
        temp->q = 1.0;
        temp->b = 0.0;
        temp->chargeDisp.setZero();
        strcpy(temp->name, "Na\0");
        cations.push_back(temp->index);
        numOfCations++;
            
        this->distances.resize(numOfParticles);
        std::for_each(distances.begin(), distances.end(), [=](std::vector<double> &row){ row.resize(numOfParticles); });
        particles.push_back(std::move(*temp));

        update_distances(*temp);
        delete temp;

        if(cations.size() != numOfCations || anions.size() != numOfAnions){
            printf("Sizes are wrong in add, cat: vec %lu, int %i an: vec %lu, int %i!\n", cations.size(), numOfCations, anions.size(), numOfAnions);
            exit(1);
        }
    }


    bool add(double q){
        Particle* temp = new Particle();
        temp->d = 5.0;
        temp->com[0] = (double) rand()/RAND_MAX * (-Base::xL) + Base::xL / 2.0;
        temp->com[1] = (double) rand()/RAND_MAX * (-Base::yL) + Base::yL / 2.0;
        //particles[i]->com[2] = (double) rand()/RAND_MAX * (Base::zL - particles[i]->d - 2 * Base::wall) + particles[i]->d/2.0 + Base::wall;
        temp->com[2] = (double) rand()/RAND_MAX * -2.0 * (Base::zLBox / 2.0 - temp->d / 2.0) + Base::zLBox / 2.0 - temp->d / 2.0; // 25.0 - Base::zLBox / 2.0;//
        temp->pos = temp->com;
        temp->index = numOfParticles;
        //temp->index = numOfCations;
        std::vector<double> tempVec(numOfParticles);
        if(hard_sphere(*temp)){
            this->numOfParticles++;
            temp->q = q;
            temp->b = 0.0;
            temp->chargeDisp.setZero();
            if(q > 0){
                strcpy(temp->name, "Na\0");
                cations.push_back(temp->index);
                numOfCations++;
            }
            else{
                strcpy(temp->name, "Cl\0");
                anions.push_back(temp->index);
                numOfAnions++;
            }
            

            /*
            this->distances.resize(numOfParticles);
            std::for_each(distances.begin(), distances.end(), [=](std::vector<double> &row){ row.resize(numOfParticles); });
            */
           /*
            this->distances.reserve(numOfParticles);
            std::for_each(distances.begin(), distances.end(), [=](std::vector<double> &row){ row.reserve(numOfParticles); });
            */
           
            if(distances.size() < this->numOfParticles){
                //printf("%lu\n", distances.size());
                this->distances.push_back(std::vector<double>(numOfParticles - 1));
                std::for_each(distances.begin(), distances.end(), [=](std::vector<double> &row){ row.push_back(0.0); });
            }
            particles.push_back(std::move(*temp));
          
            update_distances(*temp);
            delete temp;

            if(cations.size() != numOfCations || anions.size() != numOfAnions){
                printf("Sizes are wrong in add, cat: vec %lu, int %i an: vec %lu, int %i!\n", cations.size(), numOfCations, anions.size(), numOfAnions);
                exit(1);
            }

            return true;
        }

        else{
            delete temp;
            return false;
        }
    }





    void remove(std::size_t ind){
        numOfParticles--;
        int i;
        //printf("%lu\n", ind);

        if(particles.at(ind).q > 0){

            for(i = 0; i < cations.size(); i++){
                if(cations[i] == ind){
                    break;
                }
            }
            if(i == cations.size()){
                printf("cant find in cations...\n");
                exit(1);   
            }
            cations.erase(cations.begin() + i);
            
            for( ; i < cations.size(); i++){
                cations[i]--;
            }  

            for(i = 0 ; i < anions.size(); i++){
                if(anions[i] > ind){
                    anions[i]--;
                }
            } 

            numOfCations--; 
        }

        else{
            for(i = 0; i < anions.size(); i++){
                if(anions[i] == ind){
                    break;
                }
            }

            if(i == anions.size()){
                printf("cant find in anions...\n");
                exit(1);   
            }

            anions.erase(anions.begin() + i);

            for( ; i < anions.size(); i++){
                anions[i]--;
            }  
            for( ; i < cations.size(); i++){
                if(cations[i] > ind){
                    cations[i]--;
                }
            }  

            numOfAnions--;
        }
        

        particles.erase(particles.begin() + ind);
        distances.erase(distances.begin() + ind);
        for(int i = 0; i < distances.size(); i++){
            distances[i].erase(distances[i].begin() + ind);
        }
        for(int i = ind; i < numOfParticles; i++){
            particles.at(i).index--;
        }

        if(cations.size() != numOfCations || anions.size() != numOfAnions){
            printf("Sizes are wrong, cat in remove: vec %lu, int %i an: vec %lu, int %i!\n", cations.size(), numOfCations, anions.size(), numOfAnions);
            exit(1);
        }
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
                    anions.push_back(j);
                }

                else if(atom == "Na"){
                    strcpy(particles[j].name, "Na\0");
                    particles[j].b = 0;
                    particles[j].q = 2.0;      
                    cations.push_back(j);      
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
                anions.push_back(i);
            }

            else{
                particles[i].q = 2.0;
                particles[i].b = 0.0; //Length of charge displacement vector
                strcpy(particles[i].name, "Na\0");
                cations.push_back(i);
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
            /*if(i < 1){
                std::istringstream iss(line);
                if (!(iss >> c)) {
                    printf("Error reading file...\n");
                    exit(1); 
                } // error
            }*/
            if(i >= 0){
                
                std::istringstream iss(line);
                if (!(iss >> x >> y >> z)) {
                    break; 
                } // error
                particles.push_back(new Particle());
                //particles[j]->pos = (double*) malloc(3 * sizeof(double));
                particles[j].d = 5;
                particles[j].pos[0] = x;
                particles[j].pos[1] = y;
                particles[j].pos[2] = z - Base::zLBox / 2.0;
                particles[j].com = particles[j].pos;

                particles[j].b = 0;
                particles[j].index = j;
                cations.push_back(j);

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
                particles[j].pos[2] = z - Base::zLBox / 2.0;
                particles[j].com = particles[j].pos;

                particles[j].d = 5;
                particles[j].b = 0;
                particles[j].index = j;
                anions.push_back(j);

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
            fprintf(f, "%5d%-5s%5s%5d %.10lf %.10lf %.10lf\n", 1, "ion", particles[i].name, i, (particles[i].com[0] + Base::xL / 2.0)/10.0, (particles[i].com[1] + Base::yL / 2.0)/10.0, (particles[i].com[2] + Base::zLBox / 2.0)/10.0);
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
                distances.at(p.index).at(i) = p.distance_xy(particles[i]);
            }
            for(int i = 0; i < p.index; i++){
                distances.at(i).at(p.index) = p.distance_xy(particles[i]);
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