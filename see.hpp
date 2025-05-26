#pragma once

#ifndef GL_SILENCE_DEPRECATION
#define GL_SILENCE_DEPRECATION
#endif

#include <GL/freeglut.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "BuriedPipe.hpp"

// toofus
#include "AABB.hpp"
#include "ColorTable.hpp"
#include "delaunay2D.hpp"
#include "fileTool.hpp"
#include "glTools.hpp"
#include "mat4.hpp"

struct ParticleData {
  mat4r Sigma;
  double volume;
  double q;
  double p;
};

struct PipeNodeData {
  double angleRad; // position angle of the node with respect to x horizontal axis
  double N;        // axial force
  double M;        // bending moment
  double externalHoopStress;
};

struct ZoneData {
  vec2r corner[4];
  mat4r Sigma;
  double K{0.0};
  double pressure{0.0};
  size_t nbParticles{0};

  double volume() {
    double h1 = 0.5 * fabs(cross(corner[1] - corner[0], corner[3] - corner[0]));
    double h2 = 0.5 * fabs(cross(corner[1] - corner[2], corner[3] - corner[2]));
    return h1 + h2;
  }

  bool insideCell(const mat4r &H) {
    ZoneData theCell;
    theCell.corner[0].set(0.0, 0.0);
    theCell.corner[1].set(H.xx, H.yx);
    theCell.corner[2].set(H.xx + H.xy, H.yx + H.yy);
    theCell.corner[3].set(H.xy, H.yy);
    bool res = theCell.isInside(corner[0]) && theCell.isInside(corner[1]) && theCell.isInside(corner[2]) &&
               theCell.isInside(corner[3]);
    return res;
  }

  bool isInside(const vec2r &point) {
    // Compute cross products for each edge
    double cp1 = cross(corner[1] - corner[0], point - corner[0]); // crossProduct(corner0, corner1, point);
    double cp2 = cross(corner[2] - corner[1], point - corner[1]); // crossProduct(corner1, corner2, point);
    double cp3 = cross(corner[3] - corner[2], point - corner[2]); // crossProduct(corner2, corner3, point);
    double cp4 = cross(corner[0] - corner[3], point - corner[3]); // crossProduct(corner3, corner0, point);

    // Determine if the point is on the same side of all edges
    bool side1 = (cp1 >= 0) && (cp2 >= 0) && (cp3 >= 0) && (cp4 >= 0);
    bool side2 = (cp1 <= 0) && (cp2 <= 0) && (cp3 <= 0) && (cp4 <= 0);

    return side1 || side2;
  }
};

std::vector<ParticleData> pdata;

// the conf files
BuriedPipe Conf;
int confNum = 1;

AABB worldBox;

int main_window;

#define PLOT_AXIAL_FORCE 0
#define PLOT_BENDING_MOMENT 1
#define PLOT_HOOPSTRESS 2

#define MAT_AROUND_PRESSURE 0
#define MAT_AROUND_K 1

// flags
int show_background = 1; // not used
int show_particles = 1;
int show_ghosts = 0;
int show_cell = 1;
int show_forces = 0;
int show_contacts = 0;
int showOrientations = 0;
int show_pipe = 1;
int show_pipe_nodes = 0;

int show_pipe_plot = 0;
int pipe_plot_option = PLOT_HOOPSTRESS;

int show_material_around = 0;
int material_around_option = MAT_AROUND_K;

int color_option = 0;
ColorTable colorTable;

GLfloat alpha_particles = 1.0f;
GLfloat alpha_ghosts = 0.15f;

double ghost_width = 0.15;
double arrowSize = 0.0005;
double arrowAngle = 0.7;
double vScale = 0.01;

double forceTubeFactor = 1.0;

int width = 800;
int height = 800;
float wh_ratio = (float)width / (float)height;
glTextZone textZone(1, &width, &height);

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int display_mode = 0; // sample or slice rotation
int mouse_start[2];

// Drawing functions
void setColor(int i, GLfloat alpha); // this will set a color depending on the selected option
void drawForces();
void drawContacts();
void drawBox();
void drawParticles();
void drawPipe();
void drawGhosts();

// Graphics
void renderMaterialAround(); // this function draw stress information in material around the pipe
void plotPipe(); // this function makes a polaar plot where 0 is the pipe mean axis and negative values are outgoing

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);
void menu(int num);

// Helper functions
void buildMenu();
void printHelp();
void fit_view();
bool try_to_readConf(int num, BuriedPipe &CF, int &OKNum);
void updateTextLine();
void add_ghost_pos(int i, double mn, double mx, std::vector<vec2r> &lst);
void computeParticleData();
