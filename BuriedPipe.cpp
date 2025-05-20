// Periodic boundary conditions with disks and a deformable pipe

#include "BuriedPipe.hpp"

BuriedPipe::BuriedPipe() {
  // default values?
}

void BuriedPipe::saveConf(int i) {
  char fname[256];
  snprintf(fname, sizeof(fname), "conf%d", i);
  saveConf(fname);
}

void BuriedPipe::saveConf(const char *fname) {
  std::ofstream conf(fname);

  conf << "BuriedPipe 2025" << std::endl; // format: progName version-date
  conf << "t " << t << std::endl;
  conf << "tmax " << tmax << std::endl;
  conf << "dt " << dt << std::endl;
  conf << "interClose " << interClose << std::endl;
  conf << "interOut " << interOut << std::endl;
  conf << "interHist " << interHist << std::endl;
  conf << "dVerlet " << dVerlet << std::endl;
  conf << "constrainedInFrame " << constrainedInFrame << std::endl;
  conf << "density " << density << std::endl;
  conf << "kn " << kn << std::endl;
  conf << "kt " << kt << std::endl;
  conf << "dampRate " << dampRate << std::endl;
  conf << "numericalDampingCoeff " << numericalDampingCoeff << std::endl;
  conf << "mu " << mu << std::endl;
  conf << "iconf " << iconf << std::endl;
  conf << "h " << Cell.h << std::endl;
  conf << "vh " << Cell.vh << std::endl;
  conf << "ah " << Cell.ah << std::endl;
  conf << "hmass " << Cell.mass << std::endl;
  conf << "Load " << Load.StoredCommand << std::endl;

  conf << std::setprecision(15);

  conf << "Pipe " << pipe.pos.size() << ' ' << pipe.radius << ' ' << pipe.node_radius << ' ' << pipe.node_mass << ' '
       << pipe.L0 << ' ' << pipe.angle0 << ' ' << pipe.k_stretch << ' ' << pipe.k_bend << ' ' << pipe.yieldMoment
       << std::endl;
  for (size_t i = 0; i < pipe.pos.size(); i++) {
    conf << pipe.pos[i] << " " << pipe.vel[i] << " " << pipe.acc[i] << " " << pipe.force[i] << " " << pipe.moment[i]
         << " " << pipe.L[i] << " " << pipe.angle[i] << std::endl;
  }

  conf << "Particles " << Particles.size() << std::endl;
  for (size_t i = 0; i < Particles.size(); i++) {
    conf << Particles[i].pos << " " << Particles[i].vel << " " << Particles[i].acc << " " << Particles[i].rot << " "
         << Particles[i].vrot << " " << Particles[i].arot << " " << Particles[i].radius << " " << Particles[i].inertia
         << " " << Particles[i].mass << std::endl;
  }

  size_t nb_inter_saved = 0;
  for (size_t i = 0; i < Interactions.size(); i++) {
    if (fabs(Interactions[i].fn) < 1.0e-15) {
      continue;
    }
    nb_inter_saved++;
  }
  conf << "Interactions " << nb_inter_saved << std::endl;
  for (size_t i = 0; i < Interactions.size(); i++) {
    if (fabs(Interactions[i].fn) < 1.0e-15) {
      continue;
    }
    conf << Interactions[i].i << " " << Interactions[i].j << " " << Interactions[i].fn << " " << Interactions[i].ft
         << " " << Interactions[i].damp << std::endl;
  }

  size_t nb_inter_pipe_saved = 0;
  for (size_t i = 0; i < InteractionsPipe.size(); i++) {
    if (fabs(InteractionsPipe[i].fn) < 1.0e-15) {
      continue;
    }
    nb_inter_pipe_saved++;
  }
  conf << "InteractionsPipe " << nb_inter_pipe_saved << std::endl;
  for (size_t i = 0; i < InteractionsPipe.size(); i++) {
    if (fabs(InteractionsPipe[i].fn) < 1.0e-15) {
      continue;
    }
    conf << InteractionsPipe[i].i << " " << InteractionsPipe[i].inode << " " << InteractionsPipe[i].proj_div_L << " "
         << InteractionsPipe[i].fn << " " << InteractionsPipe[i].ft << " " << InteractionsPipe[i].damp << std::endl;
  }
}

void BuriedPipe::loadConf(const char *name) {
  std::ifstream conf(name);
  if (!conf.is_open()) {
    std::cerr << "Cannot read " << name << std::endl;
  }

  // Check header
  std::string prog;
  conf >> prog;
  if (prog != "BuriedPipe") {
    std::cerr << "This is not file for BuriedPipe executable!" << std::endl;
  }
  std::string date;
  conf >> date;
  if (date != "2025") {
    std::cerr << "The version-date should be 2025!" << std::endl;
  }

  fnMax = 0.0;

  std::string token;
  conf >> token;
  while (conf.good()) {
    if (token == "t") {
      conf >> t;
    } else if (token == "tmax") {
      conf >> tmax;
    } else if (token == "dt") {
      conf >> dt;
    } else if (token == "interClose") {
      conf >> interClose;
    } else if (token == "interOut") {
      conf >> interOut;
    } else if (token == "interHist") {
      conf >> interHist;
    } else if (token == "dVerlet") {
      conf >> dVerlet;
    } else if (token == "density") {
      conf >> density;
    } else if (token == "constrainedInFrame") {
      conf >> constrainedInFrame;
    } else if (token == "kn") {
      conf >> kn;
    } else if (token == "kt") {
      conf >> kt;
    } else if (token == "dampRate") {
      conf >> dampRate;
    } else if (token == "vBarrier") {
      conf >> vBarrier;
      vBarrier = fabs(vBarrier);
    } else if (token == "vBarrierExpo") {
      conf >> vBarrierExpo;
      vBarrierExpo = fabs(vBarrierExpo);
    } else if (token == "numericalDampingCoeff") {
      conf >> numericalDampingCoeff;
    } else if (token == "mu") {
      conf >> mu;
    } else if (token == "iconf") {
      conf >> iconf;
    } else if (token == "h") {
      conf >> Cell.h;
    } else if (token == "vh") {
      conf >> Cell.vh;
    } else if (token == "ah") {
      conf >> Cell.ah;
    } else if (token == "hmass") {
      conf >> Cell.mass;
    } else if (token == "Load") {
      std::string command;
      conf >> command;
      if (command == "BiaxialCompression") {
        double pressure, velocity;
        conf >> pressure >> velocity;
        Load.BiaxialCompression(pressure, velocity);
      } else if (command == "IsostaticCompression") {
        double pressure;
        conf >> pressure;
        Load.IsostaticCompression(pressure);
      } else if (command == "VelocityControl") {
        double vhxx, vhxy, vhyx, vhyy;
        conf >> vhxx >> vhxy >> vhyx >> vhyy;
        Load.VelocityControl(vhxx, vhxy, vhyx, vhyy);
      } else if (command == "SimpleShear") {
        double pressure, gammaDot;
        conf >> pressure >> gammaDot;
        Load.SimpleShear(pressure, gammaDot);
      } else {
        std::cerr << "Unknown command for loading: " << command << std::endl;
      }
    } else if (token == "Pipe") {
      size_t nb;
      conf >> nb >> pipe.radius >> pipe.node_radius >> pipe.node_mass >> pipe.L0 >> pipe.angle0 >> pipe.k_stretch >>
          pipe.k_bend >> pipe.yieldMoment;
      pipe.build_and_init(pipe.radius, 2.0 * pipe.node_radius, nb);
      vec2r POS, VEL, ACC, FORCE;
      double MOMENT, LEN, ANGLE;
      for (size_t i = 0; i < nb; i++) {
        conf >> POS >> VEL >> ACC >> FORCE >> MOMENT >> LEN >> ANGLE;
        pipe.addNode(POS, VEL, ACC);
        pipe.force[i] = FORCE;
        pipe.moment[i] = MOMENT;
        pipe.L[i] = LEN;
        pipe.angle[i] = ANGLE;
      }
    } else if (token == "Particles") {
      size_t nb;
      conf >> nb;
      Particles.clear();
      Particle P;
      radiusMin = 1e20;
      radiusMax = -1e20;
      radiusMean = 0.0;
      for (size_t i = 0; i < nb; i++) {
        conf >> P.pos >> P.vel >> P.acc >> P.rot >> P.vrot >> P.arot >> P.radius >> P.inertia >> P.mass;
        if (P.radius > radiusMax) {
          radiusMax = P.radius;
        }
        if (P.radius < radiusMin) {
          radiusMin = P.radius;
        }
        radiusMean += P.radius;
        Particles.push_back(P);
      }
      if (Particles.size() > 0) {
        radiusMean /= (double)(Particles.size());
      }
    } else if (token == "Interactions") {
      size_t nb;
      conf >> nb;
      Interactions.clear();
      Interaction I;
      for (size_t i = 0; i < nb; i++) {
        conf >> I.i >> I.j >> I.fn >> I.ft >> I.damp;
        if (I.fn > fnMax) {
          fnMax = I.fn;
        }
        Interactions.push_back(I);
      }
    } else if (token == "InteractionsPipe") {
      size_t nb;
      conf >> nb;
      InteractionsPipe.clear();
      InteractionPipe I;
      for (size_t i = 0; i < nb; i++) {
        conf >> I.i >> I.inode >> I.proj_div_L >> I.fn >> I.ft >> I.damp;
        if (I.fn > fnMax) {
          fnMax = I.fn;
        }
        InteractionsPipe.push_back(I);
      }
    } else {
      std::cerr << "Unknown token: " << token << std::endl;
    }

    conf >> token;
  }

  pipe.update(Cell.h);
}

void BuriedPipe::setSample() {
  std::ifstream file("prepa.txt");
  if (!file.is_open()) {
    std::cerr << "Cannot read 'prepa.txt'" << std::endl;
  }

  auto trim = [](const std::string &str) -> std::string {
    size_t first = str.find_first_not_of(' ');
    size_t last = str.find_last_not_of(' ');
    return (first == std::string::npos) ? "" : str.substr(first, (last - first + 1));
  };

  std::string line;
  Particle P;
  P.rot = 0.0;

  int ngw = 15;
  std::cout << std::endl;
  file >> ngw;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << ngw << std::endl;
  double step = 1.0 / (2.0 * ngw); // in range 0 to 1

  double radius = 1e-3;
  file >> radius;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << radius << std::endl;
  double deltaR = 0.2e-3;
  file >> deltaR;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << deltaR << std::endl;
  vec2r from(0.0, 0.0);
  vec2r to(2.0 * (ngw + 1) * radius, 2.0 * (ngw + 1) * radius);
  Cell.Define(to.x, from.y, from.x, to.y);
  dVerlet = 0.95 * (radius - deltaR);

  double pipeRadius = 20e-3;
  file >> pipeRadius;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << pipeRadius << std::endl;

  double pipeThickness = 2e-3;
  file >> pipeThickness;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << pipeThickness << std::endl;

  int pipeNbNodes = 36;
  file >> pipeNbNodes;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << pipeNbNodes << std::endl;
  pipe.build_and_init(pipeRadius, pipeThickness, pipeNbNodes);

  double pipe_h = 1.0; // out-of-plane size
  file >> pipe_h;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << pipe_h << std::endl;

  double sigY = 1.0;
  file >> sigY;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << sigY << std::endl;

  double Young = 1.0;
  file >> Young;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << Young << std::endl;

  double Poisson = 1.0;
  file >> Poisson;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << Poisson << std::endl;

  double pipeDensity = 3000;
  file >> pipeDensity;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << pipeDensity << std::endl;

  pipe.radius = pipeRadius;
  pipe.node_radius = 0.5 * pipeThickness;
  pipe.node_mass = M_PI * (pipeRadius * pipeRadius - (pipeRadius - pipeThickness) * (pipeRadius - pipeThickness)) *
                   pipeDensity / static_cast<double>(pipeNbNodes);

  pipe.k_stretch = Young * (pipeThickness * pipe_h) / (1.0 - Poisson * Poisson);
  pipe.k_bend = Young * (pipeThickness * pipeThickness * pipeThickness) / (12.0 * (1.0 - Poisson * Poisson));

  double Ic = (pipeThickness * pipeThickness * pipeThickness) * pipe_h / 12.0;
  pipe.yieldMoment = 4.0 * Ic * sigY / pipeThickness;

  std::cout << "   Pipe k_stretch   = " << pipe.k_stretch << std::endl;
  std::cout << "   Pipe k_bend      = " << pipe.k_bend << std::endl;
  std::cout << "   Pipe yieldMoment = " << pipe.yieldMoment << std::endl;

  mat4r hinv = Cell.h.get_inverse();
  double da = 2.0 * M_PI / static_cast<double>(pipeNbNodes);
  vec2r zero(0.0, 0.0);
  vec2r mid(0.5, 0.5);
  double Rmean = pipeRadius - 0.5 * pipeThickness;
  for (int i = 0; i < pipeNbNodes; i++) {
    double a = i * da;
    vec2r B(Rmean * cos(a), Rmean * sin(a));
    pipe.pos.push_back(mid + hinv * B);
    pipe.vel.push_back(zero);
    pipe.acc.push_back(zero);
    pipe.u.push_back(zero);
    pipe.force.push_back(zero);
    pipe.moment.push_back(0.0);
    pipe.L.push_back(pipe.L0);
    pipe.angle.push_back(pipe.angle0);
  }

  int i = 0;
  double massTot = 0.0;
  while (P.pos.y <= 1.0) {
    P.radius = radius - deltaR * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
    P.mass = M_PI * P.radius * P.radius * density;
    massTot += P.mass;
    P.inertia = 0.5 * P.mass * P.radius * P.radius;
    int column = i % ngw;
    int row = i / ngw;
    if (row % 2 == 0) { // even row
      P.pos.x = step + 2 * column * step;
    } else { // odd row
      P.pos.x = 2 * step + 2 * column * step;
    }
    P.pos.y = step + 2 * row * step;
    double dd = norm(Cell.h * (P.pos - mid)) - 1.01 * P.radius;
    if (P.pos.y <= 1. - step && dd > pipe.radius) {
      Particles.push_back(P);
    }
    i++;
  }

  Cell.mass = massTot / sqrt((double)Particles.size()); // cell mass = mean mass of one particle

  // Set ramdom velocities to particles
  double vmax = 0.1;
  file >> vmax;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << vmax << std::endl;
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].vel.x = vmax * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
    Particles[i].vel.y = vmax * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
  }

  double press = 1000.0;
  file >> press;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << press << std::endl;
  Load.IsostaticCompression(press);

  // Parametres de simu

  density = 2700.0;
  file >> density;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << density << std::endl;

  kn = 1.e4;
  file >> kn;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << kn << ", so kn/p = " << kn / press << std::endl;

  double ktkn = 1.0;
  file >> ktkn;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << ktkn << std::endl;
  kt = ktkn * kn;

  dampRate = 0.95;
  file >> dampRate;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << dampRate << std::endl;

  mu = 0.8;
  file >> mu;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << mu << std::endl;

  dt = 1e-6;
  file >> dt;
  std::getline(file, line);
  line = trim(line);
  double mass_mini = M_PI * (radius - deltaR) * (radius - deltaR) * density;
  double kn_maxi = std::max(kn, pipe.k_stretch);
  std::cout << line << " -> " << dt << ", so dt_crit/dt = " << (M_PI * sqrt(mass_mini / kn_maxi)) / dt << std::endl;

  t = 0.0;
  iconf = 0;
  interCloseC = 0.0;
  interOutC = 0.0;
  interHistC = 0.0;
  tmax = 5.0;
  file >> tmax;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << tmax << std::endl;

  interClose = 0.01;
  file >> interClose;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << interClose << std::endl;

  interOut = 0.1;
  file >> interOut;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << interOut << std::endl;

  interHist = 0.25;
  file >> interHist;
  std::getline(file, line);
  line = trim(line);
  std::cout << line << " -> " << interHist << std::endl;

  constrainedInFrame = 1;
  numericalDampingCoeff = 0.7;

  return;
}

// The integration loop (Velocity Verlet)
void BuriedPipe::integrate() {
  double dt_2 = 0.5 * dt;
  double dt2_2 = 0.5 * dt * dt;

  ResetCloseList(dVerlet);
  ResetCloseListPipe(dVerlet);
  saveConf(iconf);

  std::ofstream fileOut("output.txt");
  while (t <= tmax) {

    if (Load.ServoFunction != nullptr) {
      Load.ServoFunction(Load, Cell);
    }

    for (size_t i = 0; i < Particles.size(); i++) {
      Particles[i].pos += dt * Particles[i].vel + dt2_2 * Particles[i].acc;
      if (constrainedInFrame == 1) {
        while (Particles[i].pos.x < 0.0) {
          Particles[i].pos.x += 1.0;
        }
        while (Particles[i].pos.x > 1.0) {
          Particles[i].pos.x -= 1.0;
        }
        while (Particles[i].pos.y < 0.0) {
          Particles[i].pos.y += 1.0;
        }
        while (Particles[i].pos.y > 1.0) {
          Particles[i].pos.y -= 1.0;
        }
      }
      Particles[i].vel += dt_2 * Particles[i].acc;

      Particles[i].rot += dt * Particles[i].vrot + dt2_2 * Particles[i].arot;
      Particles[i].vrot += dt_2 * Particles[i].arot;
    }

    // pipe nodes
    for (size_t i = 0; i < pipe.pos.size(); i++) {
      pipe.pos[i] += dt * pipe.vel[i] + dt2_2 * pipe.acc[i];
      pipe.vel[i] += dt_2 * pipe.acc[i];
    }

    if (Load.Drive.xx == ForceDriven) {
      Cell.h.xx += dt * Cell.vh.xx + dt2_2 * Cell.ah.xx;
      Cell.vh.xx += dt_2 * Cell.ah.xx;
    } else {
      Cell.h.xx += dt * Load.vh.xx;
      Cell.vh.xx = Load.vh.xx;
      Cell.ah.xx = 0.0;
    }

    if (Load.Drive.xy == ForceDriven) {
      Cell.h.xy += dt * Cell.vh.xy + dt2_2 * Cell.ah.xy;
      Cell.vh.xy += dt_2 * Cell.ah.xy;
    } else {
      Cell.h.xy += dt * Load.vh.xy;
      Cell.vh.xy = Load.vh.xy;
      Cell.ah.xy = 0.0;
    }

    if (Load.Drive.yx == ForceDriven) {
      Cell.h.yx += dt * Cell.vh.yx + dt2_2 * Cell.ah.yx;
      Cell.vh.yx += dt_2 * Cell.ah.yx;
    } else {
      Cell.h.yx += dt * Load.vh.yx;
      Cell.vh.yx = Load.vh.yx;
      Cell.ah.yx = 0.0;
    }

    if (Load.Drive.yy == ForceDriven) {
      Cell.h.yy += dt * Cell.vh.yy + dt2_2 * Cell.ah.yy;
      Cell.vh.yy += dt_2 * Cell.ah.yy;
    } else {
      Cell.h.yy += dt * Load.vh.yy;
      Cell.vh.yy = Load.vh.yy;
      Cell.ah.yy = 0.0;
    }

    pipe.update(Cell.h);
    accelerations();

    // Limit the velocity of cell boundaries (vertical and horizontal only)
    // useful only for preparation
    if (vBarrier > 0.0) {
      double xxratio = pow(fabs(Cell.vh.xx / vBarrier), vBarrierExpo);
      Cell.ah.xx *= (1.0 - xxratio) / (1.0 + xxratio);
      double yyratio = pow(fabs(Cell.vh.yy / vBarrier), vBarrierExpo);
      Cell.ah.yy *= (1.0 - yyratio) / (1.0 + yyratio);
    }

    vec2r vmean;
    for (size_t i = 0; i < Particles.size(); i++) {
      Particles[i].vel += dt_2 * Particles[i].acc;
      vmean += Particles[i].vel;
      Particles[i].vrot += dt_2 * Particles[i].arot;
    }
    for (size_t i = 0; i < pipe.pos.size(); i++) {
      pipe.vel[i] += dt_2 * pipe.acc[i];
      vmean += pipe.vel[i];
    }
    vmean /= (double)(Particles.size() + pipe.pos.size());
    // substract the mean reduced velocity
    for (size_t i = 0; i < Particles.size(); i++) {
      Particles[i].vel -= vmean;
    }
    for (size_t i = 0; i < pipe.pos.size(); i++) {
      pipe.vel[i] -= vmean;
    }

    if (Load.Drive.xx == ForceDriven) {
      Cell.vh.xx += dt_2 * Cell.ah.xx;
    }
    if (Load.Drive.xy == ForceDriven) {
      Cell.vh.xy += dt_2 * Cell.ah.xy;
    }
    if (Load.Drive.yx == ForceDriven) {
      Cell.vh.yx += dt_2 * Cell.ah.yx;
    }
    if (Load.Drive.yy == ForceDriven) {
      Cell.vh.yy += dt_2 * Cell.ah.yy;
    }

    interHistC += dt;
    interOutC += dt;
    interCloseC += dt;
    t += dt;

    // ----

    if (interCloseC >= interClose - dt_2) {
      ResetCloseList(dVerlet);
      ResetCloseListPipe(dVerlet);
      interCloseC = 0.0;
    }

    if (interOutC >= interOut - dt_2) {
      fileOut << t << " " << Cell.h << " " << Sig << std::endl;
      interOutC = 0.0;
    }

    if (interHistC >= interHist - dt_2) {
      iconf++;
      std::cout << "iconf = " << iconf << ", Time = " << t << std::endl;
      saveConf(iconf);
      interHistC = 0.0;
    }
  }

  return;
}

void BuriedPipe::ResetCloseList(double dmax) {
  // store ft because the list will be cleared before being rebuilt
  std::vector<force_backup_t> force_backup;
  force_backup_t I;
  for (size_t k = 0; k < Interactions.size(); k++) {
    I.i = Interactions[k].i;
    I.j = Interactions[k].j;
    I.fn = Interactions[k].fn;
    I.ft = Interactions[k].ft;
    force_backup.push_back(I);
  }

  // now rebuild the list
  Interactions.clear();
  for (size_t i = 0; i < Particles.size(); i++) {
    for (size_t j = i + 1; j < Particles.size(); j++) {

      vec2r sij = Particles[j].pos - Particles[i].pos;
      sij.x -= floor(sij.x + 0.5);
      sij.y -= floor(sij.y + 0.5);
      vec2r branch = Cell.h * sij;

      double sum = dmax + Particles[i].radius + Particles[j].radius;
      if (norm2(branch) <= sum * sum) {
        double m = (Particles[i].mass * Particles[j].mass) / (Particles[i].mass + Particles[j].mass);
        double Damp = dampRate * 2.0 * sqrt(kn * m);
        Interactions.push_back(Interaction(i, j, Damp));
      }
    }
  }

  // retrieve ft values
  size_t k, kold = 0;
  for (k = 0; k < Interactions.size(); ++k) {
    while (kold < force_backup.size() && force_backup[kold].i < Interactions[k].i) {
      ++kold;
    }
    if (kold == force_backup.size()) {
      break;
    }

    while (kold < force_backup.size() && force_backup[kold].i == Interactions[k].i &&
           force_backup[kold].j < Interactions[k].j) {
      ++kold;
    }
    if (kold == force_backup.size()) {
      break;
    }

    if (force_backup[kold].i == Interactions[k].i && force_backup[kold].j == Interactions[k].j) {
      Interactions[k].fn = force_backup[kold].fn;
      Interactions[k].ft = force_backup[kold].ft;
      ++kold;
    }
  }
}

void BuriedPipe::ResetCloseListPipe(double dmax) {
  // Here, we suppose that the pipe is far enough from any periodic boundary
  // So, no periodic images is considered in this function

  // store ft because the list will be cleared before being rebuilt
  std::vector<force_backup_t> force_backup;
  force_backup_t I;
  for (size_t k = 0; k < InteractionsPipe.size(); k++) {
    I.i = InteractionsPipe[k].i;
    I.j = InteractionsPipe[k].inode;
    I.fn = InteractionsPipe[k].fn;
    I.ft = InteractionsPipe[k].ft;
    force_backup.push_back(I);
  }

  // now rebuild the list
  InteractionsPipe.clear();
  pipe.update(Cell.h);
  for (size_t i = 0; i < Particles.size(); i++) {
    // if (to far) {continue;}
    for (size_t in = 0; in < pipe.pos.size(); in++) {
      vec2r n = pipe.u[in].quarterRightTurned(); //(pipe.u[in].y, -pipe.u[in].x);
      vec2r Sij = Cell.h * (pipe.pos[in] - Particles[i].pos);
      double distn = fabs(Sij * n) - Particles[i].radius - pipe.node_radius;
      double projt = -Sij * pipe.u[in];
      if (projt >= -pipe.node_radius - dmax && projt < pipe.L[in] + pipe.node_radius + dmax && distn <= dmax) {
        double m = (Particles[i].mass * pipe.node_mass) / (Particles[i].mass + pipe.node_mass);
        double Damp = dampRate * 2.0 * sqrt(kn * m);
        InteractionsPipe.push_back(InteractionPipe(i, in, Damp));
      }
    }
  }

  // retrieve ft values
  size_t k, kold = 0;
  for (k = 0; k < InteractionsPipe.size(); ++k) {
    while (kold < force_backup.size() && force_backup[kold].i < InteractionsPipe[k].i) {
      ++kold;
    }
    if (kold == force_backup.size()) {
      break;
    }

    while (kold < force_backup.size() && force_backup[kold].i == InteractionsPipe[k].i &&
           force_backup[kold].j < InteractionsPipe[k].inode) {
      ++kold;
    }
    if (kold == force_backup.size()) {
      break;
    }

    if (force_backup[kold].i == InteractionsPipe[k].i && force_backup[kold].j == InteractionsPipe[k].inode) {
      InteractionsPipe[k].fn = force_backup[kold].fn;
      InteractionsPipe[k].ft = force_backup[kold].ft;
      ++kold;
    }
  }
}

void BuriedPipe::computeForces_particle_particle() {
  size_t i, j;
  for (size_t k = 0; k < Interactions.size(); k++) {
    i = Interactions[k].i;
    j = Interactions[k].j;

    vec2r sij = Particles[j].pos - Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r branch = Cell.h * sij;

    double sum = Particles[i].radius + Particles[j].radius;
    if (norm2(branch) <= sum * sum) { // it means that i and j are in contact

      vec2r n = branch;
      double len = n.normalize();
      vec2r t(-n.y, n.x);

      // real relative velocities
      vec2r vel = Particles[j].vel - Particles[i].vel;
      vec2r realVel = Cell.h * vel + Cell.vh * sij;
      double dn = len - Particles[i].radius - Particles[j].radius;
      double Bi = Particles[i].radius + 0.5 * dn;
      double Bj = Particles[j].radius + 0.5 * dn;

      // Normal force (elastic + viscuous)
      double vn = realVel * n;
      double fne = -kn * dn;
      double fnv = -Interactions[k].damp * vn;
      Interactions[k].fn = fne + fnv;

      // Tangential force (friction)
      double vijt = realVel * t - Particles[i].vrot * Bi - Particles[j].vrot * Bj;

      double ft = Interactions[k].ft - kt * dt * vijt;
      double ftest = mu * fne;
      if (fabs(ft) > ftest) {
        ft = (ft > 0.0) ? ftest : -ftest;
      }
      Interactions[k].ft = ft;

      // Resultant force and moment
      vec2r f = Interactions[k].fn * n + Interactions[k].ft * t;
      Particles[i].force -= f;
      Particles[j].force += f;
      Particles[i].moment -= ft * Bi;
      Particles[j].moment -= ft * Bj;

      // Internal stress
      Sig.xx += f.x * branch.x;
      Sig.xy += f.x * branch.y;
      Sig.yx += f.y * branch.x;
      Sig.yy += f.y * branch.y;
    } else {
      Interactions[k].fn = 0.0;
      Interactions[k].ft = 0.0;
    }
  } // Loop over particle interactions
}

void BuriedPipe::computeForces_particle_pipe() {
  // ====================================
  // Loop over the interactions with pipe
  // ====================================
  size_t i;
  size_t inode;
  for (size_t k = 0; k < InteractionsPipe.size(); k++) {
    i = InteractionsPipe[k].i;
    inode = InteractionsPipe[k].inode;

    vec2r sij = pipe.pos[inode] - Particles[i].pos;
    vec2r branch = Cell.h * sij;
    // remember we do not consider periodicity here (we are far from the
    // boundaries)

    double proj = -branch * pipe.u[inode];
    InteractionsPipe[k].proj_div_L = proj / pipe.L[inode];
    vec2r n = pipe.u[inode].quarterLeftTurned(); //(-pipe.u[inode].y, pipe.u[inode].x);
    double dn = fabs(branch * n) - Particles[i].radius - pipe.node_radius;

    if (dn < 0.0 && proj >= 0.0 && proj <= pipe.L[inode]) { // contact with flat side

      double xi_next = proj / pipe.L[inode];
      double xi = 1.0 - xi_next;

      size_t inode_next = (inode + 1) % pipe.pos.size();
      sij = (xi * pipe.pos[inode] + xi_next * pipe.pos[inode_next]) - Particles[i].pos;
      branch = Cell.h * sij;

      // real relative velocities
      vec2r vel = (xi * pipe.vel[inode] + xi_next * pipe.vel[inode_next]) - Particles[i].vel;
      vec2r realVel = Cell.h * vel + Cell.vh * sij;
      double Bi = Particles[i].radius + 0.5 * dn;

      // Normal force (elastic + viscuous)
      double vn = realVel * n;
      double fne = -kn * dn;
      double fnv = -InteractionsPipe[k].damp * vn;
      InteractionsPipe[k].fn = fne + fnv;

      // Tangential force (friction)
      double vijt = realVel * pipe.u[inode] - Particles[i].vrot * Bi;

      double ft = InteractionsPipe[k].ft - kt * dt * vijt;
      double ftest = mu * fne;
      if (fabs(ft) > ftest) {
        ft = (ft > 0.0) ? ftest : -ftest;
      }
      InteractionsPipe[k].ft = ft;

      // Resultant force and moment
      vec2r f = InteractionsPipe[k].fn * n + InteractionsPipe[k].ft * pipe.u[inode];
      Particles[i].force -= f;
      pipe.force[inode] += xi * f;
      pipe.force[inode_next] += xi_next * f;
      Particles[i].moment -= ft * Bi;

      // Internal stress
      Sig.xx += f.x * branch.x;
      Sig.xy += f.x * branch.y;
      Sig.yx += f.y * branch.x;
      Sig.yy += f.y * branch.y;

    } else if (proj < 0.0) { // contact with the vextex disk (starting one, not
                             // the end)

      double sum = Particles[i].radius + pipe.node_radius;
      if (norm2(branch) <= sum * sum) { // it means that i and j are in contact

        vec2r n = branch;
        double len = n.normalize();
        vec2r t(-n.y, n.x);

        // real relative velocities
        vec2r vel = pipe.vel[inode] - Particles[i].vel;
        vec2r realVel = Cell.h * vel + Cell.vh * sij;
        double dn = len - Particles[i].radius - pipe.node_radius;
        double Bi = Particles[i].radius + 0.5 * dn;

        // Normal force (elastic + viscuous)
        double vn = realVel * n;
        double fne = -kn * dn;
        double fnv = -InteractionsPipe[k].damp * vn;
        InteractionsPipe[k].fn = fne + fnv;

        // Tangential force (friction)
        double vijt = realVel * t - Particles[i].vrot * Bi;

        double ft = InteractionsPipe[k].ft - kt * dt * vijt;
        double ftest = mu * fne;
        if (fabs(ft) > ftest) {
          ft = (ft > 0.0) ? ftest : -ftest;
        }
        InteractionsPipe[k].ft = ft;

        // Resultant force and moment
        vec2r f = InteractionsPipe[k].fn * n + InteractionsPipe[k].ft * t;
        Particles[i].force -= f;
        pipe.force[inode] += f;
        Particles[i].moment -= ft * Bi;

        // Internal stress
        Sig.xx += f.x * branch.x;
        Sig.xy += f.x * branch.y;
        Sig.yx += f.y * branch.x;
        Sig.yy += f.y * branch.y;
      }

    } else {
      InteractionsPipe[k].fn = 0.0;
      InteractionsPipe[k].ft = 0.0;
    }
  } // end loop over particle-pipe interactions
}

void BuriedPipe::computeForces_internal_pipe() {
  // beam-like forces and moments
  for (size_t i = 0; i < pipe.pos.size(); i++) {
    size_t inext = (i + 1) % pipe.pos.size();
    size_t iprev = (i - 1 + pipe.pos.size()) % pipe.pos.size();

    // linear spring
    double dl = pipe.L[i] - pipe.L0;
    double fn = -pipe.k_stretch * dl; // positive is repulsive

    vec2r f = fn * pipe.u[i]; // from i to inext
    pipe.force[i] -= f;
    pipe.force[inext] += f;

    // angular spring
    pipe.moment[i] = pipe.k_bend * (pipe.angle[i] - pipe.angle0);

    /*
        // TODO: PLASTICITY in rotations
        if (pipe.moment[i] > pipe.yieldMoment) {
          pipe.moment[i] = pipe.yieldMoment;
        } else if (pipe.moment[i] < -pipe.yieldMoment) {
          pipe.moment[i] = -pipe.yieldMoment;
        }
    */

    // moment/2 on next
    vec2r finc_next = (0.5 * pipe.moment[i] / pipe.L[i]) * pipe.u[i].quarterLeftTurned();
    pipe.force[inext] += finc_next;
    pipe.force[i] -= finc_next;

    // -moment/2 on prev
    vec2r finc_prev = (0.5 * pipe.moment[i] / pipe.L[iprev]) * pipe.u[iprev].quarterLeftTurned();
    pipe.force[iprev] += finc_prev;
    pipe.force[i] -= finc_prev;
  }

  // internal pressure (TODO)
}

void BuriedPipe::accelerations() {
  // Set all forces and moments to zero
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].force.reset();
    Particles[i].moment = 0.0;
  }

  for (size_t i = 0; i < pipe.force.size(); i++) {
    pipe.force[i].reset();
    pipe.moment[i] = 0.0;
  }

  Sig.reset();

  computeForces_particle_particle();
  computeForces_particle_pipe();
  computeForces_internal_pipe();

  // ======================================================
  // Damping and finishing the computation of accelerations
  // ======================================================
  double invV = 1.0 / Cell.h.det();
  Sig *= invV;

  // Cundall damping for particles
  double factor = 0.0;
  double factorMinus = 0.0;
  double factorPlus = 0.0;
  if (numericalDampingCoeff > 0.0) {
    factor = 0.0;
    factorMinus = 1.0 - numericalDampingCoeff;
    factorPlus = 1.0 + numericalDampingCoeff;

    for (size_t i = 0; i < Particles.size(); i++) {
      factor = (Particles[i].force.x * Particles[i].vel.x > 0.0) ? factorMinus : factorPlus;
      Particles[i].force.x *= factor;
      factor = (Particles[i].force.y * Particles[i].vel.y > 0.0) ? factorMinus : factorPlus;
      Particles[i].force.y *= factor;
    }
  }

  //  TODO ?? Cundall damping for the pipe

  // Finally compute the accelerations (translation and rotation)
  // Compute inverse of the matrix h
  mat4r hinv = Cell.h.get_inverse();
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].acc = hinv * (Particles[i].force / Particles[i].mass);
    Particles[i].arot = Particles[i].moment / Particles[i].inertia;
  }

  double inv_pipe_node_mass = 1.0 / pipe.node_mass;
  for (size_t i = 0; i < pipe.force.size(); i++) {
    pipe.acc[i] = hinv * (pipe.force[i] * inv_pipe_node_mass);
  }

  if (Load.Drive.xx == ForceDriven) {
    Cell.ah.xx = ((Sig.xx - Load.Sig.xx) * Cell.h.yy - (Sig.yx - Load.Sig.yx) * Cell.h.xy) / Cell.mass;
    if (numericalDampingCoeff > 0.0) {
      Cell.ah.xx *= (Cell.ah.xx * Cell.vh.xx > 0.0) ? (factorMinus) : (factorPlus);
    }
  }
  if (Load.Drive.xy == ForceDriven) {
    Cell.ah.xy = ((Sig.xy - Load.Sig.xy) * Cell.h.yy - (Sig.yy - Load.Sig.yy) * Cell.h.xy) / Cell.mass;
    if (numericalDampingCoeff > 0.0) {
      Cell.ah.xy *= (Cell.ah.xy * Cell.vh.xy > 0.0) ? (factorMinus) : (factorPlus);
    }
  }
  if (Load.Drive.yx == ForceDriven) {
    Cell.ah.yx = ((Sig.yx - Load.Sig.yx) * Cell.h.xx - (Sig.xx - Load.Sig.xx) * Cell.h.yx) / Cell.mass;
    if (numericalDampingCoeff > 0.0) {
      Cell.ah.yx *= (Cell.ah.yx * Cell.vh.yx > 0.0) ? (factorMinus) : (factorPlus);
    }
  }
  if (Load.Drive.yy == ForceDriven) {
    Cell.ah.yy = ((Sig.yy - Load.Sig.yy) * Cell.h.xx - (Sig.xy - Load.Sig.xy) * Cell.h.yx) / Cell.mass;
    if (numericalDampingCoeff > 0.0) {
      Cell.ah.yy *= (Cell.ah.yy * Cell.vh.yy > 0.0) ? (factorMinus) : (factorPlus);
    }
  }
}
