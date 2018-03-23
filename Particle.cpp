#include "particle.h"
#include <iostream>
#include <math.h>
#include "ran2_lib.cpp"

//Constructor
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
    if(zP1 - zP2 < -1 * zL/2){
        zP2 = zP2 - zL;
    }
    if(zP1 - zP2 > zL/2){
        zP2 = zP2 + zL;
    }
    return (pow((xP1 - xP2), 2) + pow((yP1 - yP2), 2) + pow((zP1 - zP2), 2));
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

    return distance * distance;
}

int Particle::hardSphere(Particle **particles){
    int i = 0;
    int j = 0;

    for(i = 0; i < Particle::numOfParticles; i++){
        if(i != this->index){
            if(this->distance_xy(particles[i]) < (this->d+particles[i]->d)/2*(this->d+particles[i]->d)/2){
                return 0;
            }
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
    Particle **particles;
    particles = (Particle**) malloc(num * sizeof(Particle*));

    for(i = 0; i < num; i++){
        particles[i] = new Particle();
        particles[i]->d = 5;
        particles[i]->pos = (double*)malloc(3*sizeof(double));
        particles[i]->pos[0] = (double) rand()/RAND_MAX * Base::xL;
        particles[i]->pos[1] = (double) rand()/RAND_MAX * Base::yL;
        particles[i]->pos[2] = (double) rand()/RAND_MAX * (Base::zL - 5 - 2 * Base::wall) + 5.0/2.0 + Base::wall;
        
        particles[i]->index = i;

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

Particle** Particle::read_coordinates(std::string name){
    int i = 0;
    int j = 0;
    double x, y, z;
    int c;
    std::string line;
    std::ifstream infile(name);
    Particle** particles;
    while (std::getline(infile, line))
    {
        if(i < 1){
            std::istringstream iss(line);
            if (!(iss >> c)) {
                break; 
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
            particles[j]->pos[0] = x;
            particles[j]->pos[1] = y;
            particles[j]->pos[2] = z;
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
        fprintf(f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "liq", particles[i]->name, i, particles[i]->pos[0], particles[i]->pos[1], particles[i]->pos[2]);
    }
    fprintf(f, "%lf    %lf     %lf\n", Base::xL, Base::yL, Base::zL);
    fclose(f);
}
