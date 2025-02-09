#include "Pipe.hpp"

Pipe::Pipe()
    : radius(0.0),
      node_radius(0.0),
      node_mass(0.0),
      L0(0.0),
      angle0(0.0),
      k_stretch(0.0),
      k_bend(0.0),
      yieldMoment(0.0) {}

void Pipe::reset(size_t N) {
  pos.clear();
  vel.clear();
  acc.clear();
  u.clear();
  L.clear();
  angle.clear();
  force.clear();
  // L0 = 2.0 * M_PI * (radius - node_radius) / (double)N;
  L0 = 2.0 * (radius - node_radius) * sin(M_PI / (double)N);
  angle0 = 2.0 * M_PI / (double)N;
}

void Pipe::addNode(vec2r& POS, vec2r& VEL, vec2r& ACC) {
  pos.push_back(POS);
  vel.push_back(VEL);
  acc.push_back(ACC);
  u.push_back(vec2r::zero());
  L.push_back(L0);
  angle.push_back(angle0);
  force.push_back(vec2r::zero());
}

void Pipe::update(mat4r& h) {
  // unit vectors and lengths of the sides
  for (size_t i = 0; i < pos.size(); i++) {
    size_t i_next = (i + 1) % pos.size();
    u[i] = h * (pos[i_next] - pos[i]);
    L[i] = u[i].normalize();
  }

  // angles
  for (size_t i = 0; i < pos.size(); i++) {
    size_t i_next = (i + 1) % pos.size();
    size_t i_prev = (i - 1 + pos.size()) % pos.size();
    // angle[i] = angle(u[i_prev], u[i_next]);
    angle[i] = atan2(u[i_prev].x * u[i_next].y - u[i_next].x * u[i_prev].y,
                     u[i_prev].x * u[i_next].x + u[i_prev].y * u[i_next].y);
  }
}