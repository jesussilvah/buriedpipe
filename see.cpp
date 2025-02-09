#include "see.hpp"

void printHelp() {
  using namespace std;
  cout << endl;
  cout << "+         load next configuration file" << endl;
  cout << "-         load previous configuration file" << endl;
  cout << "=         fit the view" << endl;
  cout << "q         quit" << endl;
  // cout << "" << endl;
  cout << endl;
}

void printInfo() {
  using namespace std;

  cout << "Reference Conf = " << refConfNum << "\n";
  cout << "Current Conf = " << confNum << "\n";
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

    case '0': {
      color_option = 0;
    } break;

    case '1': {
      colorTable.setMinMax(0.5, 1.0);
      colorTable.setTableID(2);
      colorTable.Rebuild();
      color_option = 1;
    } break;

    case '2': {
      colorTable.setMinMax(0.5, 1.0);
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

    case 'd': {
      show_displacements = 1 - show_displacements;
      if (show_displacements == 1 && show_fluctuations == 1) {
        show_fluctuations = 0;
      }
    } break;

    case 'f': {
      show_fluctuations = 1 - show_fluctuations;
      if (show_fluctuations == 1 && show_displacements == 1) {
        show_displacements = 0;
      }
    } break;

    case 'g': {
      show_ghosts = 1 - show_ghosts;
    } break;

    case 'i': {
      printInfo();
    } break;

    case 'n': {
      std::cout << "Go to file number ";
      int conNumTry;
      std::cin >> conNumTry;
      bool fileRead = try_to_readConf(conNumTry, Conf, confNum);
      if (fileRead == true) {
        triangulate();
      }

    } break;

    case 'o': {
      showOrientations = 1 - showOrientations;
    } break;

    case 'p': {
      show_particles = 1 - show_particles;
    } break;

    case 'q': {
      exit(0);
    } break;

    case 'r': {
      ref_fixed = 1 - ref_fixed;
      updateTextLine();
    } break;

    case 'S': {
      vScale *= 1.05;
    } break;

    case 's': {
      vScale *= 0.95;
      if (vScale < 0.0) vScale = 1.0;
    } break;

    case 'w': {
      ghost_width = std::max(0.0, ghost_width - 0.05);
    } break;
    case 'W': {
      ghost_width = std::min(0.5, ghost_width + 0.05);
    } break;

    case '-': {
      int deltaNum = confNum - refConfNum;

      if (ref_fixed == 0 && refConfNum > 0) {
        std::cout << "Reference Configuration: ";
        try_to_readConf(refConfNum - 1, RefConf, refConfNum);
      }
      if (confNum - refConfNum > deltaNum || (ref_fixed == 1 && refConfNum == 0)) {
        std::cout << "Current Configuration: ";
        bool fileRead = try_to_readConf(confNum - 1, Conf, confNum);
        if (fileRead == true) {
          triangulate();
        }
      }

      updateTextLine();
    } break;

    case '+': {
      std::cout << "Current Configuration: ";
      bool fileRead = try_to_readConf(confNum + 1, Conf, confNum);
      if (fileRead == true) {
        triangulate();
      }
      if (fileRead == true && ref_fixed == 0) {
        std::cout << "Reference Configuration: ";
        try_to_readConf(refConfNum + 1, RefConf, refConfNum);
      }

      updateTextLine();
    } break;

    case '=': {
      fit_view();
      reshape(width, height);
    } break;
  };

  glutPostRedisplay();
}

void updateTextLine() {
  textZone.addLine("%s[%d](%0.4g s) -> [%d](%0.4g s), delta = [%d](%0.4g s)", (ref_fixed == 1) ? "*" : " ", refConfNum,
                   RefConf.t, confNum, Conf.t, confNum - refConfNum, Conf.t - RefConf.t);
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
          mouse_mode = PAN;
        else
          mouse_mode = ROTATION;
        break;
      case GLUT_MIDDLE_BUTTON:
        mouse_mode = ZOOM;
        break;
    }
  }
}

void motion(int x, int y) {

  if (mouse_mode == NOTHING) return;

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
  display();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  if (show_cell) {
    drawBox();
  }

  if (show_particles == 1) {
    drawParticles();
  }

  drawPipe();
  // drawContacts();
  drawForces();

  if (show_ghosts == 1) {
    drawGhosts();
  }

  // drawMesh();

  if (show_displacements == 1) {
    drawDisplacements();
  }

  if (show_fluctuations == 1) {
    drawFluctuations();
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

void setColor(int /*i*/, GLfloat alpha) {
  switch (color_option) {

    case 0: {
      glColor4f(0.8f, 0.8f, 0.9f, alpha);
    } break;

    case 1: {
      /*
      colorRGBA col;
      colorTable.getRGB(Conf.grains[i].zncc1, &col);
      glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
      */
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

void add_ghost_pos(int i, double mn, double mx, std::vector<vec2r>& lst) {
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

void drawPipe() {
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
  /*
  for (size_t i = 0; i < Conf.pipe.pos.size(); i++) {
    pos = Conf.Cell.h * Conf.pipe.pos[i];
    glBegin(GL_LINE_LOOP);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
      glVertex2f(pos.x + Conf.pipe.node_radius * cos(angle), pos.y + Conf.pipe.node_radius * sin(angle));
    }
    glEnd();
  }
  */
}

void drawParticles() {
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
  // if (mouse_mode != NOTHING && box.Particles.size() > 2000) return;

  std::vector<vec2r> lst_pos;  // list of reduced positions of ghost particles
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
    perp *= 0.5*width;

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
    perp *= 0.5*width;

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

void drawDisplacements() {
  glLineWidth(1.5f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  glBegin(GL_LINES);
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r displ = pos - RefConf.Cell.h * RefConf.Particles[i].pos;

    glVertex2f(pos.x - displ.x, pos.y - displ.y);
    glVertex2f(pos.x, pos.y);
  }
  glEnd();
}

void drawFluctuations() {
  glLineWidth(1.5f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  vec2r vecFluctMean;
  std::vector<vec2r> allFluct(Conf.Particles.size());
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r delta = Conf.Particles[i].pos - RefConf.Particles[i].pos;
    delta.x -= floor(delta.x + 0.5);
    delta.y -= floor(delta.y + 0.5);
    allFluct[i] = Conf.Cell.h * delta;
    vecFluctMean.x += fabs(allFluct[i].x);
    vecFluctMean.y += fabs(allFluct[i].y);
  }
  vecFluctMean /= (double)(Conf.Particles.size());
  double normalisationFactor = 1.0 / norm(vecFluctMean);

  glBegin(GL_LINES);
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r fluct = allFluct[i] * normalisationFactor;
    fluct *= vScale;

    glVertex2f(pos.x - fluct.x, pos.y - fluct.y);
    glVertex2f(pos.x, pos.y);
  }
  glEnd();
}

void triangulate() {
  std::vector<Point2D> points;
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos;
    points.push_back(Point2D(pos.x, pos.y));
  }
  triangles.clear();
  Delaunay TRIANGULATION(points);
  triangle_element T;
  double L2mini = 4 * Conf.radiusMin * Conf.radiusMin;
  for (int i = 0; i < TRIANGULATION.ntree; i++) {
    if (TRIANGULATION.thelist[i].stat > 0) {
      T.p0 = TRIANGULATION.thelist[i].p[0];
      T.p1 = TRIANGULATION.thelist[i].p[1];
      T.p2 = TRIANGULATION.thelist[i].p[2];
      // il faudrait tester si le triangle n'est pas trop plat ou d'autres critères
      // avant le l'ajouter
      vec2r x0 = Conf.Cell.h * Conf.Particles[T.p0].pos;
      vec2r x1 = Conf.Cell.h * Conf.Particles[T.p1].pos;
      vec2r x2 = Conf.Cell.h * Conf.Particles[T.p2].pos;

      if (norm2(x1 - x0) < L2mini || norm2(x2 - x0) < L2mini || norm2(x2 - x1) < L2mini) {
        continue;
      }
      triangles.push_back(T);
    }
  }
  computeTriangleStrain();
}

void computeTriangleStrain() {
  if (triangles.empty()) {
    return;
  }

  I2min = 1e20;
  I2max = -1e20;

  for (size_t i = 0; i < triangles.size(); ++i) {
    vec2r x0 = Conf.Cell.h * Conf.Particles[triangles[i].p0].pos;
    vec2r x1 = Conf.Cell.h * Conf.Particles[triangles[i].p1].pos;
    vec2r x2 = Conf.Cell.h * Conf.Particles[triangles[i].p2].pos;

    vec2r X0 = RefConf.Cell.h * RefConf.Particles[triangles[i].p0].pos;
    vec2r X1 = RefConf.Cell.h * RefConf.Particles[triangles[i].p1].pos;
    vec2r X2 = RefConf.Cell.h * RefConf.Particles[triangles[i].p2].pos;

    mat9r X(X0.x, X0.y, 1.0, X1.x, X1.y, 1.0, X2.x, X2.y, 1.0);
    mat9r x(x0.x, x0.y, 1.0, x1.x, x1.y, 1.0, x2.x, x2.y, 1.0);
    mat9r X_inv = X.get_inverse();
    mat9r F = X_inv * x;

    triangles[i].F.xx = F.xx;
    triangles[i].F.xy = F.xy;
    triangles[i].F.yx = F.yx;
    triangles[i].F.yy = F.yy;
    triangles[i].translation.x = F.zx;
    triangles[i].translation.y = F.zy;

    mat4r C = triangles[i].F.transposed() * triangles[i].F;  // Cauchy-Green droit
    triangles[i].E = 0.5 * (C - mat4r::unit());
    triangles[i].I1 = C.trace();
    triangles[i].I2 = 0.5 * (triangles[i].I1 * triangles[i].I1 - (C * C).trace());
    if (triangles[i].I2 < I2min) {
      I2min = triangles[i].I2;
    }
    if (triangles[i].I2 > I2max) {
      I2max = triangles[i].I2;
    }
    triangles[i].I3 = C.det();
  }
  // std::cout << I2min << " " << I2max << std::endl;
}

void drawMesh() {
  if (triangles.empty()) {
    return;
  }

  glLineWidth(1.0f);
  triangleColorTable.setMinMax(I2min, 1. /*I2max*/);
  // triangleColorTable.setMinMax(0.0, 0.2);

  for (size_t i = 0; i < triangles.size(); ++i) {
    vec2r pos0 = Conf.Cell.h * Conf.Particles[triangles[i].p0].pos;
    vec2r pos1 = Conf.Cell.h * Conf.Particles[triangles[i].p1].pos;
    vec2r pos2 = Conf.Cell.h * Conf.Particles[triangles[i].p2].pos;

    color4f col;
    triangleColorTable.getColor4f(fabs(triangles[i].I2), &col);
    // triangleColorTable.getColor4f(fabs(triangles[i].E.xy), &col);
    glColor4f(col.r, col.g, col.r, 0.2f);

    glBegin(GL_POLYGON);
    // glBegin(GL_LINE_LOOP);
    glVertex2f(pos0.x, pos0.y);
    glVertex2f(pos1.x, pos1.y);
    glVertex2f(pos2.x, pos2.y);
    glEnd();
  }
}

bool try_to_readConf(int num, BuriedPipe& CF, int& OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << file_name << std::endl;
    OKNum = num;
    CF.loadConf(file_name);
    return true;
  } else {
    std::cout << file_name << " does not exist" << std::endl;
  }
  return false;
}

void menu(int num) {
  switch (num) {

    case 0:
      exit(0);
      break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  int submenu1 = glutCreateMenu(menu);  // Particle Colors
  glutAddMenuEntry("None", 100);
  glutAddMenuEntry("Velocity Magnitude", 101);
  glutAddMenuEntry("Sum of Normal Contact Forces", 102);

  int submenu2 = glutCreateMenu(menu);  // Force Colors
  glutAddMenuEntry("None", 200);
  glutAddMenuEntry("Magnitude", 201);

  glutCreateMenu(menu);  // Main menu
  glutAddSubMenu("Particle Colors", submenu1);
  glutAddSubMenu("Force Colors", submenu2);
  glutAddSubMenu("Velocity Colors", submenu2);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {

  if (argc == 1) {
    refConfNum = confNum = 0;
    ref_fixed = 1;
  } else if (argc == 2) {
    refConfNum = confNum = atoi(argv[1]);
    ref_fixed = 1;
  } else if (argc == 3) {
    refConfNum = atoi(argv[1]);
    confNum = atoi(argv[2]);
    ref_fixed = 0;
  }

  std::cout << "Reference Configuration: ";
  try_to_readConf(refConfNum, RefConf, refConfNum);
  std::cout << "Current Configuration: ";
  bool fileRead = try_to_readConf(confNum, Conf, confNum);
  if (fileRead == true) {
    triangulate();
  }

  // triangleColorTable.setTableID(MATLAB_HOT);
  // triangleColorTable.setSwap(true);  // setInvert ?
  // triangleColorTable.Rebuild();
  // triangleColorTable.setMinMax(0.0, 0.5);

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowPosition(50, 50);
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

  mouse_mode = NOTHING;

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Other initialisations
  glText::init();
  // textZone.addLine("%s[%d] -> [%d]", (ref_fixed == 1) ? "*" : " ", refConfNum, confNum);
  updateTextLine();

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
