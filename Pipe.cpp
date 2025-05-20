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

void Pipe::build_and_init(double ext_radius, double thickness, size_t N) {
  radius = ext_radius;
  node_radius = 0.5 * thickness;
  pos.clear();
  vel.clear();
  acc.clear();
  u.clear();
  L.clear();
  angle.clear();
  force.clear();
  moment.clear();
  L0 = 2.0 * (radius - node_radius) * sin(M_PI / (double)N);
  angle0 = -2.0 * M_PI / (double)N;
}

void Pipe::addNode(vec2r& POS, vec2r& VEL, vec2r& ACC) {
  pos.push_back(POS);
  vel.push_back(VEL);
  acc.push_back(ACC);
  u.push_back(vec2r::zero());
  L.push_back(L0);
  angle.push_back(angle0);
  force.push_back(vec2r::zero());
  moment.push_back(0.0);
}

// update bar lengths and angles
void Pipe::update(mat4r& h) {
  // unit vectors and lengths of the sides (bar)
  for (size_t i = 0; i < pos.size(); i++) {
    size_t i_next = (i + 1) % pos.size();
    u[i] = h * (pos[i_next] - pos[i]);
    L[i] = u[i].normalize();
  }

  // node angles
  for (size_t i = 0; i < pos.size(); i++) {
    size_t i_prev = (i - 1 + pos.size()) % pos.size();
    angle[i] = angleBetweenVectors(u[i], u[i_prev]);
  }
}

// update the shape in real world coordinates relative to pipe center
void Pipe::updateShape(mat4r& h) {
  center.reset();
  for (size_t i = 0; i < pos.size(); i++) {
    center += pos[i];
  }
  center /= (double)pos.size();

  CartesianPos.clear();
  CartesianPos.resize(pos.size());
  for (size_t i = 0; i < pos.size(); i++) {
    CartesianPos[i] = pos[i] - center;
    CartesianPos[i] = h * CartesianPos[i];
  }
  center = h * center;

  PolarPos.clear();
  PolarPos.resize(pos.size());
  for (size_t i = 0; i < pos.size(); i++) {
    PolarPos[i][0] = norm(CartesianPos[i]);
    double A = atan2(CartesianPos[i].y, CartesianPos[i].x);
    if (A < 0.0) {
      A += 2 * M_PI;
    }
    PolarPos[i][1] = A;
  }
}
