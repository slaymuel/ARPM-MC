#include "particle.h"
#include <iostream>
#include <math.h>
#include "ran2_lib.cpp"

int Particle::numOfParticles = 0;
int Particle::numOfElectrons = 0;
double **Particle::distances;
double Particle::oldEnergy;

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
    if(this->com[0] > xL){
        this->com[0] = this->com[0] - xL;
    }
    if(this->com[0] < 0){
        this->com[0] = this->com[0] + xL;
    }
    if(this->com[1] > yL){
        this->com[1] = this->com[1] - yL;
    }
    if(this->com[1] < 0){
        this->com[1] = this->com[1] + yL;
    }
    if(this->com[2] > zL){
        this->com[2] = this->com[2] - zL;
    }
    if(this->com[2] < 0){
        this->com[2] = this->com[2] + zL;
    }
}

void Particle::pbc_pos(){
    //Translate charges according to periodic boundary conditions
 
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

void Particle::pbc(Eigen::Vector3d& x){
    //Translate particles according to periodic boundary conditions
    if(x[0] > Base::xL){
        x[0] -= Base::xL;
    }
    if(x[0] < 0){
        x[0] += Base::xL;
    }
    if(x[1] > Base::yL){
        x[1] -= Base::yL;
    }
    if(x[1] < 0){
        x[1] += Base::yL;
    }
    if(x[2] > Base::zL){
        x[2] -= Base::zL;
    }
    if(x[2] < 0){
        x[2] += Base::zL;
    }
}

void Particle::pbc_xy(Eigen::Vector3d& x){
    //Translate particles according to periodic boundary conditions in the xy-directions
    if(x[0] > Base::xL){
        x[0] -= Base::xL;
    }
    if(x[0] < 0){
        x[0] += Base::xL;
    }
    if(x[1] > Base::yL){
        x[1] -= Base::yL;
    }
    if(x[1] < 0){
        x[1] += Base::yL;
    }
}

// void Particle::randomMove(){
//     this->pos[0] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pos[1] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pos[2] += (double) rand()/RAND_MAX*0.5 - 0.25;
//     this->pbc();
// }

void Particle::random_move(double stepSize){
    this->com[0] += (ran2::get_random()*2.0 - 1.0) * stepSize;
    this->com[1] += (ran2::get_random()*2.0 - 1.0) * stepSize;
    this->com[2] += (ran2::get_random()*2.0 - 1.0) * stepSize;

    // this->pos[0] = this->com[0] + this->chargeDisp[0];
    // this->pos[1] = this->com[1] + this->chargeDisp[1];
    // this->pos[2] = this->com[2] + this->chargeDisp[2];
    if(Base::wall > 0 || Base::d2){
        pbc_xy(this->com);
        this->pos = this->com + this->chargeDisp;
        pbc_xy(this->pos);
    }
    else{
        pbc(this->com);
        this->pos = this->com + this->chargeDisp;
        pbc(this->pos);
    }
}

void Particle::random_charge_rot(){
    double da = 0.8;
    double random = ran2::get_random();
    if(random >= 0.1){
        this->chargeDisp[0] += (ran2::get_random() * 2.0 - 1.0) * da;
        this->chargeDisp[1] += (ran2::get_random() * 2.0 - 1.0) * da;
        this->chargeDisp[2] += (ran2::get_random() * 2.0 - 1.0) * da;
    }
    else{
        this->chargeDisp[0] = ran2::get_random() * 2.0 - 1.0;
        this->chargeDisp[1] = ran2::get_random() * 2.0 - 1.0;
        this->chargeDisp[2] = ran2::get_random() * 2.0 - 1.0;        
    }
    this->chargeDisp = this->b * this->chargeDisp.normalized();
    this->pos = this->com + this->chargeDisp;
    pbc(this->pos);
}

double Particle::distance(Particle *p){
    //Calculate distance between particles
    Eigen::Vector3d disp;
    disp = p->pos - this->pos;

    if(disp[0] < -1 * xL/2){
        disp[0] += xL;
    }
    if(disp[0] > xL/2){
        disp[0] -= xL;
    }
    if(disp[1] < -1 * yL/2){
        disp[1] += yL;
    }
    if(disp[1] > yL/2){
        disp[1] -= yL;
    }
    if(disp[2] < -1 * zL/2){
        disp[2] += zL;
    }
    if(disp[2] > zL/2){
        disp[2] -= zL;
    }
    return disp.norm();
}

double Particle::distance_xy(Particle *p){
    //Calculate distance between particles
    double xP1 = p->pos[0];
    double yP1 = p->pos[1];
    double zP1 = p->pos[2];
    double xP2 = this->pos[0];
    double yP2 = this->pos[1];
    double zP2 = this->pos[2];

    if(xP1 - xP2 < -1.0 * xL/2.0){
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
    return sqrt(pow((xP1 - xP2), 2) + pow((yP1 - yP2), 2) + pow((zP1 - zP2), 2));
}

double Particle::distance_z(Particle *p){
    //Calculate distance between particles
    double distance = 0;
    distance = this->pos[2] - p->pos[2];

    return distance;
}

double Particle::com_distance(Particle *p){
    Eigen::Vector3d disp;
    disp = p->com - this->com;

    if(disp[0] < -1 * xL/2){
        disp[0] += xL;
    }
    if(disp[0] > xL/2){
        disp[0] -= xL;
    }
    if(disp[1] < -1 * yL/2){
        disp[1] += yL;
    }
    if(disp[1] > yL/2){
        disp[1] -= yL;
    }
    if(disp[2] < -1 * zL/2){
        disp[2] += zL;
    }
    if(disp[2] > zL/2){
        disp[2] -= zL;
    }
    return disp.norm();
}

double Particle::com_distance_xy(Particle *p){
    Eigen::Vector3d disp;
    disp = p->com - this->com;

    if(disp[0] < -1 * xL/2){
        disp[0] += xL;
    }
    if(disp[0] > xL/2){
        disp[0] -= xL;
    }
    if(disp[1] < -1 * yL/2){
        disp[1] += yL;
    }
    if(disp[1] > yL/2){
        disp[1] -= yL;
    }
    return disp.norm();
}

int Particle::hard_sphere(Particle **particles){
    int i = 0;
    int j = 0;

    for(i = 0; i < Particle::numOfParticles; i++){
        if(Base::wall > 0 || Base::d2){
            if(i != this->index){
                if(this->com_distance_xy(particles[i]) < (this->d + particles[i]->d) / 2){
                    return 0;
                }
            }
        }
        else{
            if(i != this->index){
                if(this->com_distance(particles[i]) < (this->d + particles[i]->d) / 2){
                    return 0;
                }
            }
        }
    }
    // for(i = this->index + 1; i < Particle::numOfParticles; i++){
    //     if(Particle::distances[this->index][particles[i]->index] < (this->d+particles[i]->d)/2){
    //         return 0;
    //     }
    // }
    // for(i = 0; i < this->index; i++){
    //     if(Particle::distances[particles[i]->index][this->index] < (this->d+particles[i]->d)/2){
    //         return 0;
    //     }
    // }
    return 1;
}

void Particle::place_particles(Particle **particles){
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    double random;
    double diameter = particles[0]->d;
    int numOfColumns = Base::xL/diameter - 1.0;
    int numOfRows = Base::yL/diameter - 1.0;
    int index = 0;
    int numOfCells = numOfRows * numOfColumns * (Base::zL/diameter - 1.0);
    double **grid;
    double *temp;
    temp = (double*)malloc(3 * sizeof(double));
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


Particle** Particle::create_particles(int nNum, int pNum, int eNum){

    int i = 0;
    int num = nNum + pNum;
    double norm = 0;
    Particle **particles;
    particles = (Particle**) malloc((num + eNum) * sizeof(Particle*));

    for(i = 0; i < num; i++){
        particles[i] = new Particle();
        particles[i]->index = i;
        particles[i]->d = 5;    //Diameter of particles

        if(i < nNum){
            particles[i]->q = -1.0;
            particles[i]->b = 0; //Length of charge displacement vector
            strcpy(particles[i]->name, "Cl\0");
        }
        else{
            particles[i]->q = 1.0;
            particles[i]->b = 0; //Length of charge displacement vector
            strcpy(particles[i]->name, "Na\0");
        }
        

        //Get random center of mass coordinates
        particles[i]->com[0] = (double) rand()/RAND_MAX * Base::xL;
        particles[i]->com[1] = (double) rand()/RAND_MAX * Base::yL;
        particles[i]->com[2] = (double) rand()/RAND_MAX * (Base::zL - particles[i]->d - 2 * Base::wall) + particles[i]->d/2.0 + Base::wall;

        //Get random charge displacement vector
        particles[i]->chargeDisp[0] = (double) rand()/RAND_MAX * 2 - 1;
        particles[i]->chargeDisp[1] = (double) rand()/RAND_MAX * 2 - 1;
        particles[i]->chargeDisp[2] = (double) rand()/RAND_MAX * 2 - 1;

        particles[i]->chargeDisp = particles[i]->b * particles[i]->chargeDisp.normalized();
        //Calculate position of the charge
        particles[i]->pos = particles[i]->com + particles[i]->chargeDisp;
        
        if(particles[i]->com[2] < particles[i]->d/2 || particles[i]->com[2] > Base::zL - particles[i]->d/2){
            printf("%lf %lf %lf Particle in forbidden area...\n", particles[i]->com[0], particles[i]->com[1], particles[i]->com[2]);
            exit(1);
        }
    }
    printf("\033[34mCreated %d particles.\033[30m\n", Particle::numOfParticles);
    create_electrons(particles, eNum);

    return particles;
}

void Particle::create_electrons(Particle** particles, int num){
    int i, j = 0;
    i = Particle::numOfParticles;

    for(j = 0; j < num; j++){

        particles[i] = new Particle(true);
        particles[i]->index = i;
        particles[i]->d = 0;    //Diameter of particles
        particles[i]->b = 0; //Length of charge displacement vector

        //Get random center of mass coordinates
        particles[i]->com[0] = (double) rand()/RAND_MAX * Base::xL;
        particles[i]->com[1] = (double) rand()/RAND_MAX * Base::yL;
        //if(j < num / 2){
            particles[i]->com[2] = Base::zL + 2.5;//(double) rand()/RAND_MAX * (Base::zL - 2 * Base::wall) + Base::wall + Base::zL;
        //}
        //else{
        //    particles[i]->com[2] = 0;//(double) rand()/RAND_MAX * (Base::zL - 2 * Base::wall) + Base::wall - Base::zL;
        //}

        //Get random charge displacement vector
        particles[i]->chargeDisp << 0, 0, 0;
        particles[i]->pos = particles[i]->com;
        
        if(particles[i]->com[2] > 2 * Base::zL){
            printf("%lf %lf %lf Electron in forbidden area...\n", particles[i]->com[0], particles[i]->com[1], particles[i]->com[2]);
            exit(1);
        }

        particles[i]->q = -1.0;
        strcpy(particles[i]->name, "e\0");
        i++;
    }

    printf("\033[34mCreated %d electrons.\033[30m\n", j);
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
                if(particles[i]->com_distance_xy(particles[j]) < (particles[i]->d/2 + particles[j]->d/2)){
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
            //particles[j]->pos = (double*) malloc(3 * sizeof(double));
            particles[j]->chargeDisp.setZero();
            if(relative){
                particles[j]->pos[0] = Base::xL * x;
                particles[j]->pos[1] = Base::yL * y;
                particles[j]->pos[2] = Base::zL * z;
            }
            else{
                particles[j]->com[0] = x * nano;
                particles[j]->com[1] = y * nano;
                particles[j]->com[2] = z * nano + Base::wall; 
                particles[j]->pos = particles[j]->com + particles[j]->chargeDisp;
                //pbc(particles[j]->pos);   
            }

            particles[j]->d = 5;
            particles[j]->b = 0;
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

Particle** Particle::read_coordinates_gro(std::string name){
    int i = 0;
    int j = 0;
    double x, y, z;
    int c, ind;
    double nano;
    std::string molecule;
    std::string atom;
    std::string line;
    std::ifstream infile(name);
    Particle** particles;
    while (std::getline(infile, line))
    {
        if(i == 1){
            std::istringstream iss(line);
            if (!(iss >> c)) {
                printf("The second line in the input file should be the total number of atoms!\n");
                exit(1); 
            } // error
            particles = (Particle**) malloc(c * sizeof(Particle*));
        }
        if(i > 1){
            std::istringstream iss(line);
            if (!(iss >> molecule >> atom >> ind >> x >> y >> z)) {
                printf("Done reading input file\n");
                break;
                //exit(1);
            }
            particles[j] = new Particle();

            particles[j]->com[0] = x * 10;
            particles[j]->com[1] = y * 10;
            particles[j]->com[2] = z * 10;    

            particles[j]->pos[0] = x * 10;
            particles[j]->pos[1] = y * 10;
            particles[j]->pos[2] = z * 10;  

            particles[j]->d = 5;
            particles[j]->b = 0;
            particles[j]->index = j;

            if(atom == "Cl"){
                strcpy(particles[j]->name, "Cl\0");
                particles[j]->q = -1.0;
            }
            else if(atom == "Na"){
                strcpy(particles[j]->name, "Na\0");
                particles[j]->q = 1.0;            
            }
            else{
                printf("Atom is not Na or Cl!\n");
                exit(1);
            }
            j++;
        }
        i++;
    }
    if(c != j){
        printf("Did not read all particles from file, is the second line really the number of particles?\n");
        exit(1);
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
            //particles[j]->pos = (double*) malloc(3 * sizeof(double));

            particles[j]->pos[0] = x + Base::xL/2;
            particles[j]->pos[1] = y + Base::yL/2;
            particles[j]->pos[2] = z + Base::wall - 2.5;
            particles[j]->com = particles[j]->pos;

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
            //particles[j]->pos = (double*) malloc(3 * sizeof(double));

            particles[j]->pos[0] = x + Base::xL/2;
            particles[j]->pos[1] = y + Base::yL/2;
            particles[j]->pos[2] = z + Base::wall - 2.5;
            particles[j]->com = particles[j]->pos;

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

Particle** Particle::read_arpm_jan(std::string fileName){
    int i = 0;
    int j = 0;
    double x, y, z;
    int c;
    std::string line;
    std::ifstream infile(fileName);
    Particle** particles;

    //Read positive particles
    while (std::getline(infile, line)){
        if(i < 1){
            std::istringstream iss(line);
            if (!(iss >> c)) {
                printf("Error reading file...\n");
                exit(1); 
            } // error
            particles = (Particle**) malloc(c * sizeof(Particle*));
        }
        
        if(i >= 1){
            std::istringstream iss(line);
            if (!(iss >> x >> y >> z)) {
                printf("File in wrong format...\n");
                break; 
            } // error

            //Create Cation
            particles[j] = new Particle();
            //particles[j]->pos = (double*) malloc(3 * sizeof(double));

            particles[j]->com[0] = x + Base::xL/2;
            particles[j]->com[1] = y + Base::yL/2;
            particles[j]->com[2] = z + Base::zL/2;//Base::wall - 2.5;   

            std::getline(infile, line);
            iss.str(line);
            if (!(iss >> x >> y >> z)) {
                printf("File in wrong format...\n");
                break; 
            } // error

            particles[j]->pos[0] = x + Base::xL/2;
            particles[j]->pos[1] = y + Base::yL/2;
            particles[j]->pos[2] = z + Base::zL/2;

            particles[j]->chargeDisp = particles[j]->com - particles[j]->pos;
            particles[j]->chargeDisp.normalize();
            particles[j]->d = 5;
            particles[j]->index = j;
            strcpy(particles[j]->name, "Na\0");
            particles[j]->q = 1.0;
            particles[j]->b = 1.0;
            particles[j]->pos = particles[j]->com;
            j++;

            //Create Anion
            particles[j] = new Particle();
            std::getline(infile, line);
            iss.str(line);
            if (!(iss >> x >> y >> z)) {
                printf("File in wrong format...\n");
                break; 
            } // error
            particles[j]->com[0] = x + Base::xL/2;
            particles[j]->com[1] = y + Base::yL/2;
            particles[j]->com[2] = z + Base::zL/2;

            std::getline(infile, line);
            iss.str(line);
            if (!(iss >> x >> y >> z)) {
                printf("File in wrong format...\n");
                break; 
            } // error

            particles[j]->pos[0] = x + Base::xL/2;
            particles[j]->pos[1] = y + Base::yL/2;
            particles[j]->pos[2] = z + Base::zL/2;
            particles[j]->chargeDisp = particles[j]->com - particles[j]->pos;
            particles[j]->chargeDisp.normalize();
            particles[j]->d = 5;
            particles[j]->index = j;
            strcpy(particles[j]->name, "Cl\0");
            particles[j]->q = -1.0;
            particles[j]->pos = particles[j]->com;
            particles[j]->b = 1.0;
            j++;
        }
        i++;
    }


    printf("%d particles read from file.\n", j);
    //exit(1);
    return particles;
}

void Particle::write_coordinates(char name[40], Particle **particles){
    int i = 0;
    FILE *f = fopen(name, "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    fprintf(f, "Generated by Slaymulator.\n");
    fprintf(f, "%d\n", Particle::numOfParticles + Particle::numOfElectrons);
    for(i = 0; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
        fprintf(f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "ion", particles[i]->name, i, particles[i]->com[0]/10, particles[i]->com[1]/10, particles[i]->com[2]/10);
    }
    fprintf(f, "%lf    %lf     %lf\n", Base::xL/10, Base::yL/10, Base::zL/10);
    fclose(f);
}

void Particle::write_charge_coordinates(char name[], Particle **particles){
    int i = 0;
    FILE *f = fopen(name, "w");
    if(f == NULL){
        printf("Can't open file!\n");
        exit(1);
    }
    fprintf(f, "Generated by Slaymulator.\n");
    fprintf(f, "%d\n", Particle::numOfParticles);
    for(i = 0; i < Particle::numOfParticles; i++){
        fprintf(f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "cha", "H", i, particles[i]->pos[0]/10, particles[i]->pos[1]/10, particles[i]->pos[2]/10);
    }
    fprintf(f, "%lf    %lf     %lf\n", Base::xL/10, Base::yL/10, Base::zL/10);
    fclose(f);
}

void Particle::update_distances(Particle **particles, Particle *p){
    if(Base::wall > 0 || Base::d2){
        for(int i = p->index + 1; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
            Particle::distances[p->index][i] = p->distance_xy(particles[i]);
        }
        for(int i = 0; i < p->index; i++){
            Particle::distances[i][p->index] = p->distance_xy(particles[i]);
        }
    }
    else{
        for(int i = p->index + 1; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
            Particle::distances[p->index][i] = p->distance(particles[i]);
        }
        for(int i = 0; i < p->index; i++){
            Particle::distances[i][p->index] = p->distance(particles[i]);
        }
    }
}

void Particle::update_distances(Particle **particles){
    int k = 0;
    if(Base::wall > 0 || Base::d2){
        printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING ONLY USING PBC IN TWO DIMENSIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
        for(int i = 0; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
            k = i + 1;
            while(k < Particle::numOfParticles + Particle::numOfElectrons){
                Particle::distances[i][k] = particles[i]->distance_xy(particles[k]);
                k++;
            }
        }
    }
    else{
        for(int i = 0; i < Particle::numOfParticles + Particle::numOfElectrons; i++){
            k = i + 1;
            while(k < Particle::numOfParticles + Particle::numOfElectrons){
                Particle::distances[i][k] = particles[i]->distance(particles[k]);
                k++;
            }
        }
    }
}