#pragma once

#include "base.h"
#include "particle.h"

namespace energy{ namespace imgrep {

    void set_positions(Particle **particles);
    void update_position(Particle **particles, Particle *p);
    double get_energy(Particle **particles);
    double get_particle_energy(Particle **particles, Particle *p);




}}