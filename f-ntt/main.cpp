/**
 * Systolic array emulation
 * Topologies: 1-D linear, 2-D planar
 */

#define USE_MAC_VERSION_3
#define DEBUG_NTT

#include "ntt.hpp"
#include <algorithm>
#include <iostream>

/// Compile requirements - have cmake installed.
// run the following commands to generate build files
// cmake -B build -S .
// and then compile and link
// cmake --build build
// to run, cd in to the folder where the cmake file is. then
// ./build/ntt_implementation

int main(int argc, char **argv) {

  auto input = random_vector_gen(1024);
  auto T1 = ntt_2d(input, (uint16_t)32, (uint)32, (uint16_t)4);
  //auto result = inv_ntt(ntt(input));
  return EXIT_SUCCESS;
}