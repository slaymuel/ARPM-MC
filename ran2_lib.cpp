#ifndef RAN_2_H
#define RAN_2_H

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

extern "C"
{
    float ran2_(int*);
}

class ran2{
    public:
        static double get_random(){
            int ran_input = -1*(int) time(NULL)*100*rand();
            float ran2_var = ran2_(&ran_input);
            return (double)ran2_var;
        }
};

#endif