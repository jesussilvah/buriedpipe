#pragma once

#include <vector>

#include "vec2.hpp"

// A particle which is a disk
struct Particle {
  // Only position, velocity and acceleration
  // are expressed in the 'virtual world',
  // i.e., with reduced coordinates
  vec2r pos;  // Position
  vec2r vel;  // Velocity
  vec2r acc;  // Acceleration

  double rot;   // Angular position
  double vrot;  // Angular velocity
  double arot;  // Angular acceleration

  double radius;
  double inertia;
  double mass;

  vec2r force;
  double moment;

  Particle();  // Ctor
};


