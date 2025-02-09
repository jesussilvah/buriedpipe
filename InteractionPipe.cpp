#include "InteractionPipe.hpp"

InteractionPipe::InteractionPipe() : i(0), inode(0), fn(0.0), ft(0.0), damp(0.0) {}
InteractionPipe::InteractionPipe(size_t I, size_t Inode, double Damp)
    : i(I), inode(Inode), fn(0.0), ft(0.0), damp(Damp) {}
