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
#include "fileTool.hpp"
#include "glTools.hpp"
#include "delaunay2D.hpp"
#include "mat4.hpp"
#include "mat9.hpp"

// triangulation
struct triangle_element {
  size_t p0{0};
  size_t p1{0};
  size_t p2{0};
  mat4r F;
  vec2r translation;
  mat4r E;
  double I1{0.0};
  double I2{0.0};
  double I3{0.0};
};

std::vector<triangle_element> triangles;
double I2min, I2max;
ColorTable triangleColorTable;

// the conf files
BuriedPipe RefConf;
BuriedPipe Conf;
int refConfNum = 1;
int confNum = 1;

AABB worldBox;

int main_window;

// flags
int ref_fixed = 1;
int show_background = 1;
int show_particles = 1;
int show_ghosts = 0;
int show_displacements = 0;
int show_fluctuations = 0;
int show_cell = 1;
int show_forces = 0;
int showOrientations = 0;
int showMesh = 1;

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
int display_mode = 0;  // sample or slice rotation
int mouse_start[2];

// Drawing functions
void setColor(int i, GLfloat alpha);
void drawForces();
void drawContacts();
void drawFluctuations();
void drawDisplacements();
void drawBox();
void drawParticles();
void drawPipe();
void drawGhosts();
void drawMesh();

// triangulation
void triangulate();
void computeTriangleStrain();

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
bool try_to_readConf(int num, BuriedPipe& CF, int& OKNum);
void updateTextLine();
void add_ghost_pos(int i, double mn, double mx, std::vector<vec2r> & lst);


