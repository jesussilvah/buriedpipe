#pragma once

#include <vector>

#include "mat4.hpp"
#include "vec2.hpp"

// The pipe modeled as a core-shell
struct Pipe {
  std::vector<vec2r> pos;     // Positions (reduced coordinates)
  std::vector<vec2r> vel;     // Velocities (reduced coordinates)
  std::vector<vec2r> acc;     // Accelerations (reduced coordinates)
  std::vector<vec2r> u;       // unit vectors at each side
  std::vector<double> L;      // lengths for each side (real world unit)
  std::vector<double> angle;  // angle for each vertex
  std::vector<vec2r> force;   // node forces

  double radius;       // this the external radius
  double node_radius;  // half the thickness
  double node_mass;    // lumped mass

  double L0;  // initial side length (same for each bar)
  double angle0; // initial bar-to-bar angle (same for each node)

  // mechanical parameters
  double k_stretch;
  double k_bend;
  double yieldMoment;

  Pipe();  // Ctor

  void reset(size_t N);
  void addNode(vec2r & POS, vec2r & VEL, vec2r & ACC);
  void update(mat4r& h);
};
