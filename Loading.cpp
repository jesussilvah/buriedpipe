#include <cstdio>

#include "Loading.hpp"

Loading::Loading() {}

void Loading::BiaxialCompression(double pressure, double velocity) {
  snprintf(StoredCommand, 256, "BiaxialCompression %g %g", pressure, velocity);
  Drive.xx = true;
  Drive.yy = false;
  Drive.xy = Drive.yx = false;
  Sig.xx = pressure;
  Sig.xy = Sig.yx = Sig.yy = 0.0;
  vh.yy = -velocity;
  vh.xy = vh.yx = 0.0;
  vh.xx = 0.0;  // free in fact
  ServoFunction = nullptr;
}

void Loading::IsostaticCompression(double pressure) {
  snprintf(StoredCommand, 256, "IsostaticCompression %g", pressure);
  Drive.xx = Drive.yy = ForceDriven;
  Drive.xy = Drive.yx = VelocityDriven;
  Sig.xx = Sig.yy = pressure;
  Sig.xy = Sig.yx = 0.0;
  vh.xx = vh.yy = 0.0;  // free in fact
  vh.xy = vh.yx = 0.0;
  ServoFunction = nullptr;
}

void Loading::SimpleShear(double pressure, double gammaDot) {
  snprintf(StoredCommand, 256, "SimpleShear %g %g", pressure, gammaDot);
  Drive.xx = Drive.xy = Drive.yx = VelocityDriven;
  Drive.yy = ForceDriven;
  Sig.xx = Sig.xy = Sig.yx = 0.0;
  Sig.yy = pressure;
  vh.xx = vh.yx = vh.yy = 0.0;
  vh.xy = 0.0;  // will be driven by the servoFunction
  ServoFunction = [gammaDot](Loading& load, PeriodicCell& cell) -> void { load.vh.xy = gammaDot * cell.h.yy; };
}

void Loading::VelocityControl(double Vxx, double Vxy, double Vyx, double Vyy) {
  snprintf(StoredCommand, 256, "VelocityControl %g %g %g %g", Vxx, Vxy, Vyx, Vyy);
  Drive.xx = Drive.xy = Drive.yx = Drive.yy = VelocityDriven;
  Sig.xx = Sig.xy = Sig.yx = Sig.yy = 0.0;
  vh.xx = Vxx;
  vh.xy = Vxy;
  vh.yx = Vyx;
  vh.yy = Vyy;
  ServoFunction = nullptr;
}
