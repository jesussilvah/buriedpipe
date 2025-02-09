#include "Interaction.hpp"

Interaction::Interaction() : i(0), j(0), fn(0.0), ft(0.0), damp(0.0) {}
Interaction::Interaction(size_t I, size_t J, double Damp) : i(I), j(J), fn(0.0), ft(0.0), damp(Damp) {}
