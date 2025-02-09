#pragma once

#include <cstddef>

struct Interaction {
  size_t i, j;
  double fn, ft;
  double damp;

  Interaction();
  Interaction(size_t I, size_t J, double Damp);
};

