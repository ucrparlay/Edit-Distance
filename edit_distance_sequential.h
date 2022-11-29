#ifndef EDIT_DISTANCE_SEQUENTIAL_H_
#define EDIT_DISTANCE_SEQUENTIAL_H_

#include <string>
#include "utils.h"

class EditDistanceSequential {
 public:
  size_t Solve(const parlay::sequence<uint32_t>& a, const parlay::sequence<uint32_t>& b);
};

#endif  // namespace EDIT_DISTANCE_SEQUENTIAL_H_
