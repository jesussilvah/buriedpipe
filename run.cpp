#include "BuriedPipe.hpp"

int main(int argc, char const* argv[]) {
  BuriedPipe simu;

  if (argc < 2) {
    simu.setSample();
    simu.saveConf("input.txt");
    return 0;
  } else {
    simu.loadConf(argv[1]);
  }

  std::cout << "Beginning iterations." << std::endl;
  simu.integrate();
  std::cout << "End of iterations." << std::endl;
  return 0;
}
