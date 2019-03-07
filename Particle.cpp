#include "particle.h"
#include <iostream>
#include <math.h>
#include "ran2_lib.cpp"
#include <random>


//Constructor
Particle::Particle(bool dummie){
    //Keep track of the number of particles
    if(!dummie){
        //numOfParticles++;
    }
    this->chargeDisp.setZero();
    this->b = 0.0;;
}

Particle::Particle(){

}





/*
void Particle::pbc(){
    //Translate particles according to periodic boundary conditions
    if(this->com[0] > Base::xL){
        this->com[0] = this->com[0] - Base::xL;
    }
    if(this->com[0] < 0){
        this->com[0] = this->com[0] + Base::xL;
    }
    if(this->com[1] > Base::yL){
        this->com[1] = this->com[1] - Base::yL;
    }
    if(this->com[1] < 0){
        this->com[1] = this->com[1] + Base::yL;
    }
    if(this->com[2] > zL){
        this->com[2] = this->com[2] - zL;
    }
    if(this->com[2] < 0){
        this->com[2] = this->com[2] + zL;
    }
}
*/




/*
void Particle::pbc_pos(){
    //Translate charges according to periodic boundary conditions
 
    if(this->pos[0] > Base::xL){
        this->pos[0] = this->pos[0] - Base::xL;
    }
    if(this->pos[0] < 0){
        this->pos[0] = this->pos[0] + Base::xL;
    }
    if(this->pos[1] > Base::yL){
        this->pos[1] = this->pos[1] - Base::yL;
    }
    if(this->pos[1] < 0){
        this->pos[1] = this->pos[1] + Base::yL;
    }
    if(this->pos[2] > zL){
        this->pos[2] = this->pos[2] - zL;
    }
    if(this->pos[2] < 0){
        this->pos[2] = this->pos[2] + zL;
    }
}

*/



void Particle::pbc(Eigen::Vector3d& x){
    //Translate particles according to periodic boundary conditions
    if(x[0] > Base::xL / 2.0){
        x[0] -= Base::xL;
    }
    if(x[0] < -Base::xL / 2.0){
        x[0] += Base::xL;
    }
    if(x[1] > Base::yL / 2.0){
        x[1] -= Base::yL;
    }
    if(x[1] < -Base::yL / 2.0){
        x[1] += Base::yL;
    }
    if(x[2] > Base::zLBox / 2.0){
        x[2] -= Base::zLBox;
    }
    if(x[2] < -Base::zLBox / 2.0){
        x[2] += Base::zLBox;
    }
}




void Particle::pbc_xy(Eigen::Vector3d& x){
    //Translate particles according to periodic boundary conditions in the xy-directions
    if(x[0] > Base::xLHalf){
        x[0] -= Base::xL;
    }
    if(x[0] < -Base::xLHalf){
        x[0] += Base::xL;
    }
    if(x[1] > Base::yLHalf){
        x[1] -= Base::yL;
    }
    if(x[1] < -Base::yLHalf){
        x[1] += Base::yL;
    }
}






void Particle::random_move(double stepSize){
    
    this->com[0] += (ran2::get_random() * 2.0 - 1.0) * stepSize;
    this->com[1] += (ran2::get_random() * 2.0 - 1.0) * stepSize;
    this->com[2] += (ran2::get_random() * 2.0 - 1.0) * stepSize;
    //Eigen::Vector3d temp = this->com;
    // this->pos[0] = this->com[0] + this->chargeDisp[0];
    // this->pos[1] = this->com[1] + this->chargeDisp[1];
    // this->pos[2] = this->com[2] + this->chargeDisp[2];
    if(Base::wall > 0 || Base::d2){
        
        pbc_xy(this->com);
        /*if(temp != this->com){
            std::cout << "Before: " << temp << std::endl;
            std::cout << "After: " << this->com << std::endl;
        }*/
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
    double da = 1.0;
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







double Particle::distance(Particle &p){
    //Calculate distance between particles
    Eigen::Vector3d disp;
    disp = p.pos - this->pos;

    if(disp[0] < -Base::xL / 2.0){
        disp[0] += Base::xL;
    }
    if(disp[0] > Base::xL / 2.0){
        disp[0] -= Base::xL;
    }
    if(disp[1] < -Base::yL / 2.0){
        disp[1] += Base::yL;
    }
    if(disp[1] > Base::yL / 2.0){
        disp[1] -= Base::yL;
    }
    if(disp[2] < -Base::zLBox / 2.0){
        disp[2] += zLBox;
    }
    if(disp[2] > Base::zLBox / 2.0){
        disp[2] -= zLBox;
    }
    return disp.norm();
}







double Particle::distance_xy(Particle &p){
    //Calculate distance between particles
    double xP1 = p.pos[0];
    double yP1 = p.pos[1];
    double zP1 = p.pos[2];
    double xP2 = this->pos[0];
    double yP2 = this->pos[1];
    double zP2 = this->pos[2];

    Eigen::Vector3d disp;
    disp = p.pos - this->pos;
    /*
    if(xP1 - xP2 < -1.0 * Base::xL / 2.0){
        xP2 = xP2 - Base::xL;
    }
    if(xP1 - xP2 > Base::xL / 2.0){
        xP2 = xP2 + Base::xL;
    }
    if(yP1 - yP2 < -1.0 * Base::yL / 2.0){
        yP2 = yP2 - Base::yL;
    }
    if(yP1 - yP2 > Base::yL / 2.0){
        yP2 = yP2 + Base::yL;
    }
    return sqrt(pow((xP1 - xP2), 2.0) + pow((yP1 - yP2), 2.0) + pow((zP1 - zP2), 2.0));*/
    if(disp[0] < -Base::xLHalf){
        disp[0] += Base::xL;
    }
    if(disp[0] > Base::xLHalf){
        disp[0] -= Base::xL;
    }
    if(disp[1] < -Base::yLHalf){
        disp[1] += Base::yL;
    }
    if(disp[1] > Base::yLHalf){
        disp[1] -= Base::yL;
    }
    return disp.norm();
}







double Particle::distance_z(Particle *p){
    //Calculate distance between particles
    double distance = 0;
    distance = this->pos[2] - p->pos[2];

    return distance;
}







double Particle::com_distance(Particle &p){
    Eigen::Vector3d disp;
    disp = p.com - this->com;

    if(disp[0] < -Base::xL / 2.0){
        disp[0] += Base::xL;
    }
    if(disp[0] > Base::xL / 2.0){
        disp[0] -= Base::xL;
    }
    if(disp[1] < -Base::yL / 2.0){
        disp[1] += Base::yL;
    }
    if(disp[1] > Base::yL / 2.0){
        disp[1] -= Base::yL;
    }
    if(disp[2] < -Base::zLBox / 2.0){
        disp[2] += zLBox;
    }
    if(disp[2] > Base::zLBox / 2.0){
        disp[2] -= zLBox;
    }
    return disp.norm();
}







double Particle::com_distance_xy(Particle &p){
    Eigen::Vector3d disp;
    disp = p.com - this->com;

    if(disp[0] < -Base::xLHalf){
        disp[0] += Base::xL;
    }
    if(disp[0] > Base::xLHalf){
        disp[0] -= Base::xL;
    }
    if(disp[1] < -Base::yLHalf){
        disp[1] += Base::yL;
    }
    if(disp[1] > Base::yLHalf){
        disp[1] -= Base::yL;
    }
    return disp.norm();
}

/*
void Particle::place_particles(Particle **particles){
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    double random;
    double diameter = particles[0]->d;
    int numOfColumns = Base::xL/diameter;
    int numOfRows = Base::yL/diameter;
    int index = 0;
    int numOfCells = numOfRows * numOfColumns * (Base::zL/diameter);
    double **grid;
    double *temp;
    temp = (double*)malloc(3 * sizeof(double));
    grid = (double**) malloc(numOfCells * sizeof(double*));

    printf("Fitting %d particles in box.\n", Particle::numOfParticles);
    printf("Number of Cells: %d\n", numOfCells);

    //Create grid
    for(j = 0; j < Base::zL/diameter; j++){
        for(l = 0; l < numOfRows; l++){
            for(k = 0; k < numOfColumns; k++){    
                index = l * numOfColumns + k + j * numOfColumns * numOfRows;
                grid[index] = (double*) malloc(3 * sizeof(double));
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
        random = ran2::get_random() * (numOfCells);
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
        particles[i]->com[0] = grid[i][0];
        particles[i]->com[1] = grid[i][1];
        particles[i]->com[2] = grid[i][2];
        particles[i]->pos = particles[i]->com;
        //particles[i]->d = diameter;
        //particles[i]->index = i;

        //if(i % 2 == 0){
        //    particles[i]->q = -1.0;
        //    strcpy(particles[i]->name, "Cl\0");
        //}
        //else{
        //    particles[i]->q = 1.0;
        //    strcpy(particles[i]->name, "Na\0");
        //}      
    }
}
*/

/*
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
*/

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


bool Particle::wall_2d(){
    //particles[p]->com[2] > particles[p]->d/2 + Base::wall && particles[p]->com[2] < Base::zL - Base::wall - particles[p]->d/2
    
    bool inside;

    if(this->q < 0){
        (this->com[2] > -Base::zLBox / 2.0 +       this->d / 2.0 && this->com[2] < Base::zLBox / 2.0 - this->d / 2.0) ?      inside = true : inside = false;
    }

    else{
        (this->com[2] > -Base::zLBox / 2.0 - 2.0 + this->d / 2.0 && this->com[2] < Base::zLBox / 2.0 - this->d / 2.0 + 2.0) ? inside = true : inside = false;
    }

    return inside;
}
