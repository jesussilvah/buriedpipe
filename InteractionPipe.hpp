#pragma once

#include <cstddef>

struct InteractionPipe {
  size_t i, inode;
  double fn, ft;
  double damp;

  InteractionPipe();
  InteractionPipe(size_t I, size_t Inode, double Damp);
};
