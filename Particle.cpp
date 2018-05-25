#include "particle.h"
#include <iostream>
#include <math.h>
#include "ran2_lib.cpp"

//Constructor
Particle::Particle(bool dummie){
    //Keep track of the number of particles
    if(!dummie){
        numOfParticles++;
    }
}
Particle::Particle(){
    //Keep track of the number of particles
    numOfParticles++;
}
void Particle::pbc(){
    //Translate particles according to periodic boundary conditions
    if(this->pos[0] > xL){
        this->pos[0] = this->pos[0] - xL;
    }
    if(this->pos[0] < 0){
        this->pos[0] = this->pos[0] + xL;
    }
    if(this->pos[1] > yL){
        this->pos[1] = this->pos[1] - yL;
    }
    if(this->pos[1] < 0){
        this->pos[1] = this->pos[1] + yL;
    }
    if(this->pos[2] > zL){
        this->pos[2] = this->pos[2] - zL;
    }
    if(this->pos[2] < 0){
        this->pos[2] = this->pos[2] + zL;
    }
}

void Particle::pbc_xy(){
    //Translate particles according to periodic boundary conditions in the xy-directions
    if(this->pos[0] > xL){
        this->pos[0] = this->pos[0] - xL;
    }
    if(this->pos[0] < 0){
        this->pos[0] = this->pos[0] + xL;
    }
    if(this->pos[1] > yL){
        this->pos[1] = this->pos[1] - yL;
    }
    if(this->pos[1] < 0){
        this->pos[1] = this->pos[1] + yL;
    }
}

// void Particle::randomMove(){
//     this->pos[0] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pos[1] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pos[2] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pbc();
// }

void Particle::randomMove(double stepSize){
    this->pos[0] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;
    this->pos[1] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;
    this->pos[2] += ((double)rand()/RAND_MAX * 2 - 1.0) * stepSize;

    this->pbc();
}

double Particle::distance(Particle *p){
    //Calculate distance between particles
    // double xP1 = p->pos[0];
    // double yP1 = p->pos[1];
    // double zP1 = p->pos[2];
    // double xP2 = this->pos[0];
    // double yP2 = this->pos[1];
    // double zP2 = this->pos[2];
    double xP12 = p->pos[0] - this->pos[0];
    double yP12 = p->pos[1] - this->pos[1];
    double zP12 = p->pos[2] - this->pos[2];

    if(xP12 < -1 * xL/2){
        //xP2 = xP2 - xL;
        xP12 += xL;
    }
    if(xP12 > xL/2){
        //xP2 = xP2 + xL;
        xP12 -= xL;
    }
    if(yP12 < -1 * yL/2){
        //yP2 = yP2 - yL;
        yP12 += yL;
    }
    if(yP12 > yL/2){
        //yP2 = yP2 + yL;
        yP12 -= yL;
    }
    if(zP12 < -1 * zL/2){
        //zP2 = zP2 - zL;
        zP12 += zL;
    }
    if(zP12 > zL/2){
        //zP2 = zP2 + zL;
        zP12 -= zL;
    }
    //return (pow((xP1 - xP2), 2) + pow((yP1 - yP2), 2) + pow((zP1 - zP2), 2));
    return sqrt(xP12 * xP12 + yP12 * yP12 + zP12 * zP12);
}

double Particle::distance_xy(Particle *p){
    //Calculate distance between particles
    double xP1 = p->pos[0];
    double yP1 = p->pos[1];
    double zP1 = p->pos[2];
    double xP2 = this->pos[0];
    double yP2 = this->pos[1];
    double zP2 = this->pos[2];

    if(xP1 - xP2 < -1 * xL/2){
        xP2 = xP2 - xL;
    }
    if(xP1 - xP2 > xL/2){
        xP2 = xP2 + xL;
    }
    if(yP1 - yP2 < -1 * yL/2){
        yP2 = yP2 - yL;
    }
    if(yP1 - yP2 > yL/2){
        yP2 = yP2 + yL;
    }
    return (pow((xP1 - xP2), 2) + pow((yP1 - yP2), 2) + pow((zP1 - zP2), 2));
}

double Particle::distance_z(Particle *p){
    //Calculate distance between particles
    double distance = 0;
    distance = this->pos[2] - p->pos[2];

    return distance;
}

// double Particle::com_distance(Particle *p){


// }

int Particle::hardSphere(Particle **particles){
    int i = 0;
    int j = 0;

    // for(i = 0; i < Particle::numOfParticles; i++){
    //     if(i != this->index){
    //         if(this->distance(particles[i]) < (this->d+particles[i]->d)/2*(this->d+particles[i]->d)/2){
    //             return 0;
    //         }
    //     }
    // }
    for(i = this->index + 1; i < Particle::numOfParticles; i++){
        if(Particle::distances[this->index][particles[i]->index] < (this->d+particles[i]->d)/2){
            return 0;
        }
    }
    for(i = 0; i < this->index; i++){
        if(Particle::distances[particles[i]->index][this->index] < (this->d+particles[i]->d)/2){
            return 0;
        }
    }
    return 1;
}

void Particle::place_particles(Particle **particles){
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    double random;
    double diameter = particles[0]->d;
    int numOfColumns = Base::xL/diameter - 1;
    int numOfRows = Base::yL/diameter - 1;
    int index = 0;
    int numOfCells = numOfRows * numOfColumns * (Base::zL/diameter - 1);
    double **grid;
    double *temp;
    temp = (double*)malloc(3*sizeof(double));
    grid = (double**) malloc(numOfCells*sizeof(double*));

    printf("Fitting %d particles in box.\n", Particle::numOfParticles);
    printf("Number of Cells: %d\n", numOfCells);

    //Create grid
    for(j = 0; j < Base::zL/diameter - 1; j++){
        for(l = 0; l < numOfRows; l++){
            for(k = 0; k < numOfColumns; k++){    
                index = l * numOfColumns + k + j * numOfColumns * numOfRows;
                grid[index] = (double*) malloc(3*sizeof(double));
                grid[index][0] = k * 5 + diameter/2;
                grid[index][1] = l * 5 + diameter/2;
                grid[index][2] = j * 5 + diameter/2;
                //printf("%lf %lf %lf\n", grid[index][0],grid[index][1],grid[index][2]);
                
                //grid[l][k][j] = ;
            }
        }
    }
    printf("Randomizing cells...\n");
    //Randomize cells
    for(i = 0; i < numOfCells; i++){
        random = ran2::get_random() * (numOfCells - 1);
        temp[0] = grid[i][0];
        temp[1] = grid[i][1];
        temp[2] = grid[i][2];
        grid[i][0] = grid[(int)random][0];
        grid[i][1] = grid[(int)random][1];
        grid[i][2] = grid[(int)random][2];
        grid[(int)random][0] = temp[0];
        grid[(int)random][1] = temp[1];
        grid[(int)random][2] = temp[2];
    }
    printf("Setting coordinates...\n");
    //Set particle coordinates
    for(i = 0; i < Particle::numOfParticles; i++){
        particles[i]->pos[0] = grid[i][0];
        particles[i]->pos[1] = grid[i][1];
        particles[i]->pos[2] = grid[i][2];
        particles[i]->d = diameter;
        particles[i]->index = i;

        if(i % 2 == 0){
            particles[i]->q = -1.0;
            strcpy(particles[i]->name, "Cl\0");
        }
        else{
            particles[i]->q = 1.0;
            strcpy(particles[i]->name, "Na\0");
        }       
    }
}

Particle** Particle::create_particles(int num){

    int i = 0;
    double norm = 0;
    Particle **particles;
    particles = (Particle**) malloc(num * sizeof(Particle*));

    for(i = 0; i < num; i++){
        particles[i] = new Particle();
        particles[i]->index = i;
        particles[i]->d = 5;    //Diameter of particles
        particles[i]->b = 0; //Length of charge displacement vector
        particles[i]->pos = (double*)malloc(3 * sizeof(double));

        //Get random center of mass coordinates
        particles[i]->com[0] = (double) rand()/RAND_MAX * Base::xL;
        particles[i]->com[1] = (double) rand()/RAND_MAX * Base::yL;
        particles[i]->com[2] = (double) rand()/RAND_MAX * (Base::zL - particles[i]->d - 2 * Base::wall) + particles[i]->d/2.0 + Base::wall;

        //Get random charge displacement vector
        particles[i]->chargeDisp[0] = (double) rand()/RAND_MAX * 2 - 1;
        particles[i]->chargeDisp[1] = (double) rand()/RAND_MAX * 2 - 1;
        particles[i]->chargeDisp[2] = (double) rand()/RAND_MAX * 2 - 1;

        //Normalize charge displacement vector and set length to b
        norm = sqrt(particles[i]->chargeDisp[0] * particles[i]->chargeDisp[0] +
                    particles[i]->chargeDisp[1] * particles[i]->chargeDisp[1] +
                    particles[i]->chargeDisp[2] * particles[i]->chargeDisp[2]);
        particles[i]->chargeDisp[0] = particles[i]->b * particles[i]->chargeDisp[0]/norm;
        particles[i]->chargeDisp[1] = particles[i]->b * particles[i]->chargeDisp[1]/norm;
        particles[i]->chargeDisp[2] = particles[i]->b * particles[i]->chargeDisp[2]/norm;

        //Calculate position of the charge
        particles[i]->pos[0] = particles[i]->com[0] + particles[i]->chargeDisp[0];
        particles[i]->pos[1] = particles[i]->com[1] + particles[i]->chargeDisp[1];
        particles[i]->pos[2] = particles[i]->com[2] + particles[i]->chargeDisp[2];
        
        if(particles[i]->pos[2] < particles[i]->d/2 || particles[i]->pos[2] > Base::zL - particles[i]->d/2){
            printf("%lf %lf %lf Particle in forbidden area...\n", particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
            exit(1);
        }

        if(i % 2 == 0){
            particles[i]->q = -1.0;
            strcpy(particles[i]->name, "Cl\0");
        }
        else{
            particles[i]->q = 1.0;
            strcpy(particles[i]->name, "Na\0");
        }
    }
    printf("Created %d particles.\n", num);
    return particles;
}

Particle **Particle::create_dummies(Particle **particles){
    Particle **dummies;
    dummies = (Particle**) malloc(2 * Particle::numOfParticles * (sizeof(Particle*)));

    for(int i = 0; i < Particle::numOfParticles; i++){
        dummies[i] = new Particle();
        dummies[i]->pos[0] = particles[i]->pos[0];
        dummies[i]->pos[1] = particles[i]->pos[1];
        dummies[i]->pos[2] = particles[i]->pos[2] + Base::zL;
        strcpy(dummies[i]->name, "DUM\0");

        dummies[i + 1] = new Particle();
        dummies[i + 1]->pos[0] = particles[i]->pos[0];
        dummies[i + 1]->pos[1] = particles[i]->pos[1];
        dummies[i + 1]->pos[2] = particles[i]->pos[2] - Base::zL;
        strcpy(dummies[i + 1]->name, "DUM\0");
    }

    return dummies;
}

int Particle::get_overlaps(Particle **particles){
    int overlaps = 0;
    int i = 0;
    int j = 0;

    for(i = 0; i < Particle::numOfParticles; i++){
        j = i + 1;
        while(j < Particle::numOfParticles){
            if(i != j){
                if(particles[i]->distance_xy(particles[j]) < 25){
                    overlaps++;
                }
            }
            j++;
        }
    }
    return overlaps;
}

Particle** Particle::read_coordinates(std::string name, bool relative = false, bool nanometers = false){
    int i = 0;
    int j = 0;
    double x, y, z;
    int c;
    double nano;
    if(nanometers){
        nano = 10;
    }
    else{
        nano = 1;
    }
    std::string line;
    std::ifstream infile(name);
    Particle** particles;
    while (std::getline(infile, line))
    {
        if(i < 1){
            std::istringstream iss(line);
            if (!(iss >> c)) {
                printf("Error reading file...\n");
                exit(1); 
            } // error
            particles = (Particle**) malloc(c * sizeof(Particle*));
        }
        if(i > 1){
            std::istringstream iss(line);
            if (!(iss >> name >> x >> y >> z)) {
                break; 
            } // error
            particles[j] = new Particle();
            particles[j]->pos = (double*) malloc(3 * sizeof(double));

            if(relative){
                particles[j]->pos[0] = Base::xL * x;
                particles[j]->pos[1] = Base::yL * y;
                particles[j]->pos[2] = Base::zL * z;
            }
            else{
                particles[j]->pos[0] = x * nano;
                particles[j]->pos[1] = y * nano;
                particles[j]->pos[2] = z * nano + Base::wall;    
            }

            particles[j]->d = 5;
            particles[j]->index = j;
            
            if(j % 2 == 0){
                strcpy(particles[j]->name, "Cl\0");
                particles[j]->q = -1.0;
            }
            else{
                strcpy(particles[j]->name, "Na\0");
                particles[j]->q = 1.0;            
            }
            j++;
        }
        i++;
    }
    printf("%d particles read from file.\n", j);
    return particles;
}

Particle** Particle::read_jan(std::string pName, std::string nName){
    int i = 0;
    int j = 0;
    double x, y, z;
    int c;
    std::string line;
    std::ifstream infileP(pName);
    Particle** particles;

    //Read positive particles
    while (std::getline(infileP, line)){
        if(i < 1){
            std::istringstream iss(line);
            if (!(iss >> c)) {
                printf("Error reading file...\n");
                exit(1); 
            } // error
            particles = (Particle**) malloc(2 * c * sizeof(Particle*));
        }
        if(i >= 1){
            
            std::istringstream iss(line);
            if (!(iss >> x >> y >> z)) {
                break; 
            } // error
            particles[j] = new Particle();
            particles[j]->pos = (double*) malloc(3 * sizeof(double));

            particles[j]->pos[0] = x + Base::xL/2;
            particles[j]->pos[1] = y + Base::yL/2;
            particles[j]->pos[2] = z + Base::wall - 2.5;    

            particles[j]->d = 5;
            particles[j]->index = j;
            
            strcpy(particles[j]->name, "Na\0");
            particles[j]->q = 1.0;

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
            particles[j] = new Particle();
            particles[j]->pos = (double*) malloc(3 * sizeof(double));

            particles[j]->pos[0] = x + Base::xL/2;
            particles[j]->pos[1] = y + Base::yL/2;
            particles[j]->pos[2] = z + Base::wall - 2.5;    

            particles[j]->d = 5;
            particles[j]->index = j;
            
            strcpy(particles[j]->name, "Cl\0");
            particles[j]->q = -1.0;

            j++;
        }
        i++;
    }

    printf("%d particles read from file.\n", j);
    return particles;
}

void Particle::write_coordinates(char name[], Particle **particles){
    int i = 0;
    FILE *f = fopen(name, "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    fprintf(f, "File created by .....\n");
    fprintf(f, "%d\n", Particle::numOfParticles);
    for(i = 0; i < Particle::numOfParticles; i++){
        fprintf(f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "liq", particles[i]->name, i, particles[i]->pos[0]/10, particles[i]->pos[1]/10, particles[i]->pos[2]/10);
    }
    fprintf(f, "%lf    %lf     %lf\n", Base::xL/10, Base::yL/10, Base::zL/10);
    fclose(f);
}

void Particle::update_distances(Particle *p, Particle **particles){
    for(int i = p->index + 1; i < Particle::numOfParticles; i++){
        Particle::distances[p->index][i] = p->distance(particles[i]);
    }
    for(int i = 0; i < p->index; i++){
        Particle::distances[i][p->index] = p->distance(particles[i]);
    }
}

void Particle::update_distances(Particle **particles){
    int k = 0;
    for(int i = 0; i < Particle::numOfParticles; i++){
        k = i + 1;
        while(k < Particle::numOfParticles){
            Particle::distances[i][k] = particles[i]->distance(particles[k]);
            k++;
        }
    }
}