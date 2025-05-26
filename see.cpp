#include "see.hpp"

void printHelp() {
  std::cout << std::endl;
  std::cout << "+         load next configuration file" << std::endl;
  std::cout << "-         load previous configuration file" << std::endl;
  std::cout << "=         fit the view" << std::endl;
  std::cout << "a/A       particle transparency" << std::endl;
  std::cout << "b/B       ghost particle transparency" << std::endl;
  std::cout << "c         show/hide periodic cell" << std::endl;
  std::cout << "C         show/hide contacts" << std::endl;
  std::cout << "f         show/hide normal forces" << std::endl;
  std::cout << "g         show/hide ghost particles" << std::endl;
  std::cout << "h         print this help" << std::endl;
  std::cout << "i         print information" << std::endl;
  std::cout << "m         show/hide polar plot" << std::endl;
  std::cout << "n         go to file (see terminal to enter the file number)" << std::endl;
  std::cout << "o         show/hide particle orientations" << std::endl;
  std::cout << "p         show/hide particles" << std::endl;
  std::cout << "P         switch pipe display" << std::endl;
  std::cout << "q         quit" << std::endl;
  std::cout << "s/S       tune vector sizes" << std::endl;
  std::cout << "w/W       tune displayed ghost width" << std::endl;
  std::cout << "x         show/hide surrouding material 'spider-map'" << std::endl;
  // std::cout << "x         xxxx" << std::endl;
  std::cout << std::endl;
  std::cout << "0         particles colored with light gray" << std::endl;
  std::cout << "1         particles colored with pressure" << std::endl;
  std::cout << std::endl;
}

void printInfo() {
  int V = glutGet(GLUT_VERSION);
  int V1 = V / 10000;
  int V2 = (V - V1 * 10000) / 100;
  int V3 = V - V1 * 10000 - V2 * 100;
  std::cout << "glut version " << V1 << "." << V2 << "." << V3 << "\n";
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

  case '0': {
    color_option = 0;
  } break;

  case '1': { // particle pressures
    double pmin = pdata[0].p;
    double pmax = pdata[0].p;
    for (size_t i = 1; i < pdata.size(); i++) {
      if (pdata[i].p < pmin) {
        pmin = pdata[i].p;
      }
      if (pdata[i].p > pmax) {
        pmax = pdata[i].p;
      }
    }
    colorTable.setMinMax(pmin, pmax);
    colorTable.setTableID(2);
    colorTable.Rebuild();
    color_option = 1;
  } break;

  case '2': {
    // colorTable.setMinMax(pmin, pmax);
    colorTable.setTableID(2);
    colorTable.Rebuild();
    color_option = 2;
  } break;

  case 'a': {
    alpha_particles = std::max(0.0f, alpha_particles - 0.05f);
  } break;
  case 'A': {
    alpha_particles = std::min(1.0f, alpha_particles + 0.05f);
  } break;

  case 'b': {
    alpha_ghosts = std::max(0.0f, alpha_ghosts - 0.05f);
  } break;
  case 'B': {
    alpha_ghosts = std::min(1.0f, alpha_ghosts + 0.05f);
  } break;

  case 'c': {
    show_cell = 1 - show_cell;
  } break;

  case 'C': {
    show_contacts = 1 - show_contacts;
  } break;

  case 'f': {
    show_forces = 1 - show_forces;
  } break;

  case 'g': {
    show_ghosts = 1 - show_ghosts;
  } break;

  case 'h': {
    printHelp();
  } break;

  case 'i': {
    printInfo();
  } break;

  case 'n': {
    std::cout << "Go to file number ";
    int conNumTry;
    std::cin >> conNumTry;
    try_to_readConf(conNumTry, Conf, confNum);
  } break;

  case 'm': {
    show_material_around = 1 - show_material_around;
  } break;

  case 'o': {
    showOrientations = 1 - showOrientations;
  } break;

  case 'p': {
    show_particles = 1 - show_particles;
  } break;

  case 'P': {
    show_pipe = 1 - show_pipe;
  } break;

  case 'q': {
    exit(0);
  } break;

  case 'S': {
    vScale *= 1.05;
  } break;

  case 's': {
    vScale *= 0.95;
    if (vScale < 0.0)
      vScale = 1.0;
  } break;

  case 'w': {
    ghost_width = std::max(0.0, ghost_width - 0.05);
  } break;
  case 'W': {
    ghost_width = std::min(0.5, ghost_width + 0.05);
  } break;

  case 'x': {
    show_pipe_plot = 1 - show_pipe_plot;
  } break;

  case '-': {
    if (confNum > 0) {
      try_to_readConf(confNum - 1, Conf, confNum);
    }
    updateTextLine();
  } break;

  case '+': {
    try_to_readConf(confNum + 1, Conf, confNum);
    updateTextLine();
  } break;

  case '=': {
    fit_view();
    reshape(width, height);
  } break;
  };

  glutPostRedisplay();
}

void updateTextLine() { textZone.addLine("conf %d,  t %0.4g s", confNum, Conf.t); }

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
    case GLUT_LEFT_BUTTON: {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
        mouse_mode = PAN;
      } else {
        mouse_mode = ROTATION;
      }
    } break;
    case GLUT_MIDDLE_BUTTON: {
      mouse_mode = ZOOM;
    } break;
    }
  }
}

void motion(int x, int y) {

  if (mouse_mode == NOTHING) {
    return;
  }

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

  case ZOOM: {
    double ddy = (worldBox.max.y - worldBox.min.y) * dy;
    double ddx = (worldBox.max.x - worldBox.min.x) * dy;
    worldBox.min.x -= ddx;
    worldBox.max.x += ddx;
    worldBox.min.y -= ddy;
    worldBox.max.y += ddy;
  } break;

  case PAN: {
    double ddx = (worldBox.max.x - worldBox.min.x) * dx;
    double ddy = (worldBox.max.y - worldBox.min.y) * dy;
    worldBox.min.x -= ddx;
    worldBox.max.x -= ddx;
    worldBox.min.y += ddy;
    worldBox.max.y += ddy;
  } break;

  default:
    break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  // display();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // glShadeModel(GL_SMOOTH);
  // glEnable(GL_DEPTH_TEST);

  if (show_cell) {
    drawBox();
  }

  if (show_material_around == 1) {
    renderMaterialAround();
  }

  if (show_particles == 1) {
    drawParticles();
  }

  if (show_pipe == 1) {
    drawPipe();
  }

  if (show_contacts == 1) {
    drawContacts();
  }

  if (show_forces == 1) {
    drawForces();
  }

  if (show_ghosts == 1) {
    drawGhosts();
  }

  if (show_pipe_plot == 1) {
    plotPipe();
  }

  textZone.draw();

  glFlush();
  glutSwapBuffers();
}

void fit_view() {
  //
  // 3 x ------- x 2
  //   |         |
  // 0 x ------- x 1
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = Conf.Cell.h.xy;
  double y1 = Conf.Cell.h.yy;
  double x2 = Conf.Cell.h.xy + Conf.Cell.h.xx;
  double y2 = Conf.Cell.h.yy + Conf.Cell.h.yx;
  double x3 = Conf.Cell.h.xx;
  double y3 = Conf.Cell.h.yx;
  worldBox.min.x = std::min(std::min(std::min(x0, x1), x2), x3);
  worldBox.min.y = std::min(std::min(std::min(y0, y1), y2), y3);
  worldBox.max.x = std::max(std::max(std::max(x0, x1), x2), x3);
  worldBox.max.y = std::max(std::max(std::max(y0, y1), y2), y3);
  reshape(width, height);
}

void reshape(int w, int h) {
  width = w;
  height = h;

  double left = worldBox.min.x;
  double right = worldBox.max.x;
  double bottom = worldBox.min.y;
  double top = worldBox.max.y;
  double worldW = right - left;
  double worldH = top - bottom;
  double dW = 0.1 * worldW;
  double dH = 0.1 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;

  if (worldW > worldH) {
    worldH = worldW * ((GLfloat)height / (GLfloat)width);
    top = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * ((GLfloat)width / (GLfloat)height);
    right = 0.5 * (left + right + worldW);
    left = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, bottom, top);

  glutPostRedisplay();
}

// Draw the periodic cell
void drawBox() {
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
  glLineWidth(1.0f);

  glBegin(GL_LINE_LOOP);
  glVertex2f(0.0f, 0.0f);
  glVertex2f(Conf.Cell.h.xy, Conf.Cell.h.yy);
  glVertex2f(Conf.Cell.h.xy + Conf.Cell.h.xx, Conf.Cell.h.yy + Conf.Cell.h.yx);
  glVertex2f(Conf.Cell.h.xx, Conf.Cell.h.yx);
  glEnd();
}

void setColor(int i, GLfloat alpha) {
  switch (color_option) {

  case 0: {
    glColor4f(0.8f, 0.8f, 0.9f, alpha);
  } break;

  case 1: { // pressure
    colorRGBA col;
    colorTable.getRGB(pdata[i].p, &col);
    glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
  } break;

  case 2: {
    /*
    colorRGBA col;
    colorTable.getRGB(Conf.grains[i].zncc2, &col);
    glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
    */
  } break;

  default: {
    glColor4f(0.8f, 0.8f, 0.9f, alpha);
  } break;
  }
}

void add_ghost_pos(int i, double mn, double mx, std::vector<vec2r> &lst) {
  lst.clear();
  vec2r pos = Conf.Particles[i].pos;
  if (pos.x > mn && pos.x < mx && pos.y > mn && pos.y < mx) {
    return;
  }

  vec2r ghostPos;

  if (pos.x <= mn) {
    ghostPos.set(pos.x + 1.0, pos.y);
    lst.push_back(ghostPos);

    if (pos.y <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.y >= mx) {
      ghostPos.set(pos.x + 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.x >= mx) {
    ghostPos.set(pos.x - 1.0, pos.y);
    lst.push_back(ghostPos);

    if (pos.y <= mn) {
      ghostPos.set(pos.x - 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.y >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.y <= mn) {
    ghostPos.set(pos.x, pos.y + 1.0);
    lst.push_back(ghostPos);

    if (pos.x <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.x >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.y >= mx) {
    ghostPos.set(pos.x, pos.y - 1.0);
    lst.push_back(ghostPos);

    if (pos.x <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.x >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }
}

void renderMaterialAround() {
  size_t nbNodes = Conf.pipe.pos.size();

  double zoneWidth = 2.0 * Conf.pipe.node_radius; // this will be changeable after

  std::vector<vec2r> radialDirection(nbNodes);
  std::vector<double> radiusPipe(nbNodes);
  vec2r center;
  for (size_t i = 0; i < nbNodes; i++) {
    radialDirection[i] = Conf.Cell.h * Conf.pipe.pos[i];
    center += radialDirection[i];
  }
  center /= (double)nbNodes;
  for (size_t i = 0; i < nbNodes; i++) {
    radialDirection[i] -= center;
    radiusPipe[i] = radialDirection[i].normalize() + Conf.pipe.node_radius;
  }

  std::vector<ZoneData> zones;

  for (size_t layer = 0; layer < 6; layer++) {
    for (size_t i = 0; i < nbNodes; i++) {
      size_t inext = i + 1;
      if (inext >= nbNodes) {
        inext = 0;
      }
      ZoneData Z;
      Z.corner[0] = center + (radiusPipe[i] + layer * zoneWidth) * radialDirection[i];
      Z.corner[1] = center + (radiusPipe[i] + (layer + 1) * zoneWidth) * radialDirection[i];
      Z.corner[2] = center + (radiusPipe[inext] + (layer + 1) * zoneWidth) * radialDirection[inext];
      Z.corner[3] = center + (radiusPipe[inext] + layer * zoneWidth) * radialDirection[inext];

      if (Z.insideCell(Conf.Cell.h)) {
        zones.push_back(Z);
      }
    }
  }

  // compute things
  for (size_t k = 0; k < Conf.Interactions.size(); k++) {
    size_t i = Conf.Interactions[k].i;
    size_t j = Conf.Interactions[k].j;
    double fn = Conf.Interactions[k].fn;
    double ft = Conf.Interactions[k].ft;

    vec2r sij = Conf.Particles[j].pos - Conf.Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r branch = Conf.Cell.h * sij;

    vec2r n = branch;
    double len = n.normalize();
    vec2r t(-n.y, n.x);
    vec2r f = fn * n + ft * t;

    vec2r mid = Conf.Cell.h * Conf.Particles[i].pos + 0.5 * len * n;

    for (size_t iz = 0; iz < zones.size(); iz++) {
      if (zones[iz].isInside(mid)) {
        zones[iz].Sigma.xx += f.x * branch.x;
        zones[iz].Sigma.xy += f.x * branch.y;
        zones[iz].Sigma.yx += f.y * branch.x;
        zones[iz].Sigma.yy += f.y * branch.y;
      }
    }
  }

  for (size_t k = 0; k < Conf.InteractionsPipe.size(); k++) {
    size_t i = Conf.InteractionsPipe[k].i;
    size_t inode = Conf.InteractionsPipe[k].inode;
    double fn = Conf.InteractionsPipe[k].fn;
    double ft = Conf.InteractionsPipe[k].ft;

    vec2r n = Conf.pipe.u[inode].quarterLeftTurned();
    vec2r branch = Conf.Particles[i].radius * n;
    vec2r t(-n.y, n.x);
    vec2r f = fn * n + ft * t;

    vec2r C = Conf.Cell.h * Conf.Particles[i].pos;

    for (size_t iz = 0; iz < zones.size(); iz++) {
      if (zones[iz].isInside(C)) {
        zones[iz].Sigma.xx += f.x * branch.x;
        zones[iz].Sigma.xy += f.x * branch.y;
        zones[iz].Sigma.yx += f.y * branch.x;
        zones[iz].Sigma.yy += f.y * branch.y;
      }
    }
  }

  for (size_t iz = 0; iz < zones.size(); iz++) {
    zones[iz].Sigma /= zones[iz].volume();
    zones[iz].pressure = 0.5 * (zones[iz].Sigma.xx + zones[iz].Sigma.yy);
    if (fabs(zones[iz].Sigma.yy) > 1e-8) {
      zones[iz].K = (zones[iz].Sigma.xx / zones[iz].Sigma.yy);
    }
  }

  double minval = 1e12;
  double maxval = -1e12;
  if (material_around_option == MAT_AROUND_K) {
    for (size_t iz = 0; iz < zones.size(); iz++) {
      minval = (zones[iz].K < minval) ? zones[iz].K : minval;
      maxval = (zones[iz].K > maxval) ? zones[iz].K : maxval;
    }
  } else if (material_around_option == MAT_AROUND_PRESSURE) {
    for (size_t iz = 0; iz < zones.size(); iz++) {
      minval = (zones[iz].pressure < minval) ? zones[iz].pressure : minval;
      maxval = (zones[iz].pressure > maxval) ? zones[iz].pressure : maxval;
    }
  }

  // display the zones
  ColorTable CT;
  CT.setMinMax(minval, maxval);
  colorRGBA col;
  glLineWidth(1.0f);

  for (size_t iz = 0; iz < zones.size(); iz++) {

    if (material_around_option == MAT_AROUND_K) {
      CT.getRGB(zones[iz].K, &col);
    } else if (material_around_option == MAT_AROUND_PRESSURE) {
      CT.getRGB(zones[iz].pressure, &col);
    }

    glColor3f(col.r / 255.0f, col.g / 255.0f, col.b / 255.0f);

    glBegin(GL_POLYGON);
    glVertex2f(zones[iz].corner[0].x, zones[iz].corner[0].y);
    glVertex2f(zones[iz].corner[1].x, zones[iz].corner[1].y);
    glVertex2f(zones[iz].corner[2].x, zones[iz].corner[2].y);
    glVertex2f(zones[iz].corner[3].x, zones[iz].corner[3].y);
    glEnd();
  }
}

void plotPipe() {
  if (mouse_mode != NOTHING) {
    return;
  }

  // Compute N, M, hoop-stress ...
  size_t nbNodes = Conf.pipe.pos.size();

  vec2r center;
  for (size_t i = 0; i < nbNodes; i++) {
    vec2r barPos = Conf.Cell.h * Conf.pipe.pos[i] + 0.5 * Conf.pipe.L[i] * Conf.pipe.u[i];
    center += barPos;
  }
  center /= (double)nbNodes;

  std::vector<PipeNodeData> pnd(nbNodes);
  double meanPipeRadius = 0.0;
  for (size_t i = 0; i < nbNodes; i++) {
    vec2r barPos = Conf.Cell.h * Conf.pipe.pos[i] + 0.5 * Conf.pipe.L[i] * Conf.pipe.u[i];
    barPos -= center;
    meanPipeRadius += norm(barPos);

    double A = atan2(barPos.y, barPos.x);
    if (A < 0.0) {
      A += 2 * M_PI;
    }
    pnd[i].angleRad = A;
    size_t inext = i + 1;
    if (inext >= nbNodes) {
      inext = 0;
    }
    vec2r df = Conf.pipe.force[inext] - Conf.pipe.force[i];
    pnd[i].N = df * Conf.pipe.u[i];

    pnd[i].M = cross(0.5 * Conf.pipe.L[i] * Conf.pipe.u[i], Conf.pipe.force[inext] + Conf.pipe.force[i]);
    double e = 2.0 * Conf.pipe.node_radius;
    pnd[i].externalHoopStress = (pnd[i].N / e - 6 * pnd[i].M / (e * e));
  }
  meanPipeRadius /= (double)(nbNodes);

  // plot

  glLineWidth(3.0f);
  glColor4f(1.0f, 0.1f, 0.9f, 1.0f);

  // rescale
  double Nmax = -1e12;
  double Nmin = 1e12;
  double Mmax = -1e12;
  double Mmin = 1e12;
  double Smax = -1e12;
  double Smin = 1e12;
  for (size_t i = 0; i < nbNodes; i++) {
    Nmin = (Nmin > pnd[i].N) ? pnd[i].N : Nmin;
    Nmax = (Nmax < pnd[i].N) ? pnd[i].N : Nmax;
    Mmin = (Mmin > pnd[i].M) ? pnd[i].M : Mmin;
    Mmax = (Mmax < pnd[i].M) ? pnd[i].M : Mmax;
    Smin = (Smin > pnd[i].externalHoopStress) ? pnd[i].externalHoopStress : Smin;
    Smax = (Smax < pnd[i].externalHoopStress) ? pnd[i].externalHoopStress : Smax;
  }
  double rescale = 1.0;
  if (pipe_plot_option == PLOT_AXIAL_FORCE) {
    rescale = 0.25*fabs(meanPipeRadius / (Nmax-Nmin));
  } else if (pipe_plot_option == PLOT_BENDING_MOMENT) {
    rescale = 0.25*fabs(meanPipeRadius / (Mmax-Mmin));
  } else if (pipe_plot_option == PLOT_HOOPSTRESS) {
    rescale = 0.25 * fabs(meanPipeRadius / (Smax-Smin));
  }

  glBegin(GL_LINES);
  vec2r pos, nextPos, nextPos2, prevPos;
  for (size_t i = 0; i < nbNodes; i++) {
    size_t inext = (i + 1) % nbNodes;
    size_t inext2 = (i + 2) % nbNodes;
    size_t iprev = (i - 1 + nbNodes) % nbNodes;
    pos = Conf.Cell.h * Conf.pipe.pos[i] + 0.5 * Conf.pipe.L[i] * Conf.pipe.u[i];
    nextPos = Conf.Cell.h * Conf.pipe.pos[inext] + 0.5 * Conf.pipe.L[inext] * Conf.pipe.u[inext];
    nextPos2 = Conf.Cell.h * Conf.pipe.pos[inext2] + 0.5 * Conf.pipe.L[inext2] * Conf.pipe.u[inext2];
    prevPos = Conf.Cell.h * Conf.pipe.pos[iprev] + 0.5 * Conf.pipe.L[iprev] * Conf.pipe.u[iprev];
    vec2r AxeNext = nextPos2 - nextPos;
    vec2r Axe = nextPos - pos;
    vec2r AxePrev = pos - prevPos;
    vec2r T1 = vec2r(Axe.y, -Axe.x) + vec2r(AxePrev.y, -AxePrev.x);
    T1.normalize();
    vec2r T2 = vec2r(AxeNext.y, -AxeNext.x) + vec2r(Axe.y, -Axe.x);
    T2.normalize();

    // set sign so that 'plus' is toward the center
    T1 *= -1.0;
    T2 *= -1.0;

    double f1 = 0.0, f2 = 0.0;
    if (pipe_plot_option == PLOT_AXIAL_FORCE) {
      f1 = pnd[i].N * rescale;
      f2 = pnd[inext].N * rescale;
    } else if (pipe_plot_option == PLOT_BENDING_MOMENT) {
      f1 = pnd[i].M * rescale;
      f2 = pnd[inext].M * rescale;
    } else if (pipe_plot_option == PLOT_HOOPSTRESS) {
      f1 = pnd[i].externalHoopStress * rescale;
      f2 = pnd[inext].externalHoopStress * rescale;
    }

    if (f1 > 1e-8) {
      glColor3f(1.0f, 0.0f, 0.0f);
    } else if (f1 < -1e-8) {
      glColor3f(0.0f, 0.0f, 1.0f);
    } else {
      glColor3f(0.0f, 1.0f, 0.0f);
    }
    glVertex2f(pos.x + f1 * T1.x, pos.y + f1 * T1.y);
    
    if (f2 > 1e-8) {
      glColor3f(1.0f, 0.0f, 0.0f);
    } else if (f2 < -1e-8) {
      glColor3f(0.0f, 0.0f, 1.0f);
    } else {
      glColor3f(0.0f, 1.0f, 0.0f);
    }
    glVertex2f(nextPos.x + f2 * T2.x, nextPos.y + f2 * T2.y);
  }
  glEnd();
}

void drawPipe() {
  if (mouse_mode != NOTHING) {
    return;
  }

  glLineWidth(1.0f);
  glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);

  // ...
  glBegin(GL_LINES);
  vec2r pos, nextPos;

  size_t N = Conf.pipe.pos.size();
  for (size_t i = 0; i < N; i++) {
    size_t inext = (i + 1) % N;
    // size_t prev = (i - 1 + N) % N;
    pos = Conf.Cell.h * Conf.pipe.pos[i];
    nextPos = Conf.Cell.h * Conf.pipe.pos[inext];
    vec2r T(nextPos.y - pos.y, -(nextPos.x - pos.x));
    T.normalize();
    T *= Conf.pipe.node_radius;

    glVertex2f(pos.x + T.x, pos.y + T.y);
    glVertex2f(nextPos.x + T.x, nextPos.y + T.y);

    glVertex2f(pos.x - T.x, pos.y - T.y);
    glVertex2f(nextPos.x - T.x, nextPos.y - T.y);
  }
  glEnd();

  // this draw the circles at each node
  if (show_pipe_nodes > 0) {
    for (size_t i = 0; i < Conf.pipe.pos.size(); i++) {
      pos = Conf.Cell.h * Conf.pipe.pos[i];
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + Conf.pipe.node_radius * cos(angle), pos.y + Conf.pipe.node_radius * sin(angle));
      }
      glEnd();
    }
  }
}

void drawParticles() {
  if (mouse_mode != NOTHING) {
    return;
  }

  glLineWidth(1.0f);

  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos;
    double R = Conf.Particles[i].radius;

    setColor(i, alpha_particles);
    glBegin(GL_POLYGON);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
    }
    glEnd();

    glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);
    glBegin(GL_LINE_LOOP);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
    }
    glEnd();

    if (showOrientations) {
      double rot = Conf.Particles[i].rot;
      glBegin(GL_LINES);
      glVertex2f(pos.x, pos.y);
      glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
      glEnd();
    }
  }
}

void drawGhosts() {
  if (mouse_mode != NOTHING) {
    return;
  }

  std::vector<vec2r> lst_pos; // list of reduced positions of ghost particles
  double mn = ghost_width;
  double mx = 1.0 - ghost_width;
  // GLColorRGBA color;
  glLineWidth(1.0f);
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    add_ghost_pos(i, mn, mx, lst_pos);
    double R = Conf.Particles[i].radius;
    for (size_t ig = 0; ig < lst_pos.size(); ig++) {

      vec2r pos = Conf.Cell.h * lst_pos[ig];

      setColor(i, alpha_ghosts);
      glBegin(GL_POLYGON);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();

      glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();

      if (showOrientations) {
        double rot = Conf.Particles[i].rot;
        glBegin(GL_LINES);
        glVertex2f(pos.x, pos.y);
        glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
        glEnd();
      }
    }
  }
}

void drawContacts() {
  if (mouse_mode != NOTHING) {
    return;
  }

  glLineWidth(1.5f);

  // grain-grain
  glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
  glBegin(GL_LINES);
  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i = Conf.Interactions[k].i;
    size_t j = Conf.Interactions[k].j;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.Particles[j].pos - Conf.Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r posj = posi + Conf.Cell.h * sij;
    glVertex2f(posi.x, posi.y);
    glVertex2f(posj.x, posj.y);
  }
  glEnd();

  // with pipe
  glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
  glBegin(GL_LINES);
  for (size_t k = 0; k < Conf.InteractionsPipe.size(); ++k) {
    size_t i = Conf.InteractionsPipe[k].i;
    size_t inode = Conf.InteractionsPipe[k].inode;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.pipe.pos[inode] - Conf.Particles[i].pos;
    double proj = -sij * Conf.pipe.u[inode];
    if (proj > 0.0) {
      vec2r n = Conf.pipe.u[inode].quarterLeftTurned();
      sij = (Conf.Particles[i].radius + Conf.pipe.node_radius) * n;
    } else {
      sij = Conf.Cell.h * sij;
    }
    vec2r posj = posi + sij;
    glVertex2f(posi.x, posi.y);
    glVertex2f(posj.x, posj.y);
  }
  glEnd();
}

void drawForces() {
  if (mouse_mode != NOTHING) {
    return;
  }

  // grain-grain
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i = Conf.Interactions[k].i;
    size_t j = Conf.Interactions[k].j;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.Particles[j].pos - Conf.Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r posj = posi + Conf.Cell.h * sij;

    // Calculate the width of the rectangle
    GLfloat width = Conf.radiusMean * (Conf.Interactions[k].fn / Conf.fnMax);

    // Calculate the direction vector and the perpendicular vector
    vec2r dir = posj - posi;
    vec2r perp(-dir.y, dir.x);
    perp.normalize();
    perp *= 0.5 * width;

    // Calculate the four corners of the rectangle
    vec2r p1 = posi + perp;
    vec2r p2 = posi - perp;
    vec2r p3 = posj - perp;
    vec2r p4 = posj + perp;

    // Draw the filled rectangle
    glBegin(GL_QUADS);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
    glVertex2f(p4.x, p4.y);
    glEnd();
  }

  // with pipe
  // glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
  glBegin(GL_LINES);
  for (size_t k = 0; k < Conf.InteractionsPipe.size(); ++k) {
    size_t i = Conf.InteractionsPipe[k].i;
    size_t inode = Conf.InteractionsPipe[k].inode;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.pipe.pos[inode] - Conf.Particles[i].pos;
    double proj = -sij * Conf.pipe.u[inode];
    if (proj > 0.0) {
      vec2r n = Conf.pipe.u[inode].quarterLeftTurned();
      sij = (Conf.Particles[i].radius + Conf.pipe.node_radius) * n;
    } else {
      sij = Conf.Cell.h * sij;
    }
    vec2r posj = posi + sij;

    // Calculate the width of the rectangle
    GLfloat width = Conf.radiusMean * (Conf.InteractionsPipe[k].fn / Conf.fnMax);

    // Calculate the direction vector and the perpendicular vector
    vec2r dir = posj - posi;
    vec2r perp(-dir.y, dir.x);
    perp.normalize();
    perp *= 0.5 * width;

    // Calculate the four corners of the rectangle
    vec2r p1 = posi + perp;
    vec2r p2 = posi - perp;
    vec2r p3 = posj - perp;
    vec2r p4 = posj + perp;

    // Draw the filled rectangle
    glBegin(GL_QUADS);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
    glVertex2f(p4.x, p4.y);
    glEnd();
  }
  glEnd();
}

void computeParticleData() {

  pdata.resize(Conf.Particles.size());
  for (size_t i = 0; i < pdata.size(); i++) {
    pdata[i].Sigma.reset();
    pdata[i].volume = M_PI * Conf.Particles[i].radius * Conf.Particles[i].radius;
  }

  for (size_t k = 0; k < Conf.Interactions.size(); k++) {
    size_t i = Conf.Interactions[k].i;
    size_t j = Conf.Interactions[k].j;
    vec2r sij = Conf.Particles[j].pos - Conf.Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r branch = Conf.Cell.h * sij;
    vec2r n = branch;
    double len = n.normalize();
    double dn = len - Conf.Particles[i].radius - Conf.Particles[j].radius;
    vec2r Bi = (Conf.Particles[i].radius + 0.5 * dn) * n;
    vec2r Bj = -(Conf.Particles[j].radius + 0.5 * dn) * n;

    vec2r t(-n.y, n.x);
    vec2r f = Conf.Interactions[k].fn * n + Conf.Interactions[k].ft * t;

    pdata[i].Sigma.xx += f.x * Bi.x;
    pdata[i].Sigma.xy += f.x * Bi.y;
    pdata[i].Sigma.yx += f.y * Bi.x;
    pdata[i].Sigma.yy += f.y * Bi.y;

    pdata[j].Sigma.xx -= f.x * Bj.x;
    pdata[j].Sigma.xy -= f.x * Bj.y;
    pdata[j].Sigma.yx -= f.y * Bj.x;
    pdata[j].Sigma.yy -= f.y * Bj.y;
  }

  for (size_t i = 0; i < pdata.size(); i++) {
    pdata[i].Sigma /= pdata[i].volume;
    pdata[i].p = 0.5 * (pdata[i].Sigma.xx + pdata[i].Sigma.yy);
  }
}

bool try_to_readConf(int num, BuriedPipe &CF, int &OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << file_name << std::endl;
    OKNum = num;
    CF.loadConf(file_name);
    CF.pipe.updateShape(CF.Cell.h);
    CF.accelerations();
    computeParticleData();
    return true;
  } else {
    std::cout << file_name << " does not exist" << std::endl;
  }
  return false;
}

void menu(int num) {
  switch (num) {

  case 0: {
    exit(0);
  } break;

  case 10: { // None
    show_pipe = 0;
  } break;
  case 11: { // Pipe alone
    show_pipe = 1;
    show_pipe_nodes = 0;
  } break;
  case 12: { // Show nodes
    show_pipe = 1;
    show_pipe_nodes = 1;
  } break;
  case 13: { // Show loading
    show_pipe = 1;
    show_pipe_nodes = 0;
  } break;
  case 14: { // Show internal stress
    show_pipe = 1;
    show_pipe_nodes = 0;
  } break;

  case 20: {
    show_pipe_plot = 1;
    pipe_plot_option = PLOT_AXIAL_FORCE;
  } break;
  case 21: {
    show_pipe_plot = 1;
    pipe_plot_option = PLOT_BENDING_MOMENT;
  } break;
  case 22: {
    show_pipe_plot = 1;
    pipe_plot_option = PLOT_HOOPSTRESS;
  } break;

  case 30: {
    show_material_around = 1;
    material_around_option = MAT_AROUND_K;
  } break;
  case 31: {
    show_material_around = 1;
    material_around_option = MAT_AROUND_PRESSURE;
  } break;

    // case xx: {} break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  int submenu1 = glutCreateMenu(menu);
  glutAddMenuEntry("None", 10);
  glutAddMenuEntry("Pipe alone", 11);
  glutAddMenuEntry("Show nodes", 12);
  glutAddMenuEntry("Show loading", 13);
  glutAddMenuEntry("Show internal stress", 14);

  int submenu2 = glutCreateMenu(menu);
  glutAddMenuEntry("Axial", 20);
  glutAddMenuEntry("Bending moment", 21);
  glutAddMenuEntry("External Hoop Stress", 22);

  int submenu3 = glutCreateMenu(menu);
  glutAddMenuEntry("K", 30);
  glutAddMenuEntry("pressure", 31);

  glutCreateMenu(menu); // Main menu
  glutAddSubMenu("Pipe Display Options", submenu1);
  glutAddSubMenu("Polar Plot Options", submenu2);
  glutAddSubMenu("Spider Map Options", submenu3);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char *argv[]) {

  if (argc == 1) {
    confNum = 0;
  } else if (argc == 2) {
    confNum = atoi(argv[1]);
  }

  std::cout << "Current Configuration: ";
  try_to_readConf(confNum, Conf, confNum);

  mouse_mode = NOTHING;

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  int X0 = (glutGet(GLUT_SCREEN_WIDTH) - width) / 2;
  int Y0 = (glutGet(GLUT_SCREEN_HEIGHT) - height) / 2;
  glutInitWindowPosition(X0, Y0);
  glutInitWindowSize(width, height);

  main_window = glutCreateWindow("CONF VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  // ==== Other initialisations
  glText::init();
  updateTextLine();

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
