#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Interaction.hpp"
#include "InteractionPipe.hpp"
#include "Loading.hpp"
#include "Particle.hpp"
#include "PeriodicCell.hpp"
#include "Pipe.hpp"

class BuriedPipe {
 public:
  Pipe pipe;
  std::vector<Particle> Particles;
  std::vector<Interaction> Interactions;
  std::vector<InteractionPipe> InteractionsPipe;
  Loading Load;
  PeriodicCell Cell;
  // double Sigxx, Sigxy, Sigyx, Sigyy;  //  Internal stress
  mat4r Sig;

  struct force_backup_t {
    size_t i{0}, j{0};
    double fn{0.0}, ft{0.0};
  };

  // Parameters
  double t{0.0};
  double tmax{5.0};
  double dt{1e-6};
  double interCloseC{0.0}, interClose{0.01}, dVerlet;
  double interOutC{0.0}, interOut{0.1};
  double interHistC{0.0}, interHist{0.25};
  double density{2700.0};
  double kn{1e4};
  double kt{1e4};
  double dampRate{0.95};
  double vBarrier{0.0};
  double numericalDampingCoeff{0.0};
  double vBarrierExpo{2.0};
  double mu{0.8};
  int iconf{0};
  int constrainedInFrame{0};

  double radiusMin{0.0};
  double radiusMax{0.0};
  double radiusMean{0.0};
  double fnMax{0.0};

  // Methods
  BuriedPipe();
  void setSample();
  void integrate();
  
  void computeForces_particle_particle();
  void computeForces_particle_pipe();
  void computeForces_internal_pipe();
  void accelerations();
  
  void ResetCloseList(double dmax);
  void ResetCloseListPipe(double dmax);
  
  void saveConf(int i);
  void saveConf(const char* name);
  void loadConf(const char* name);
};
