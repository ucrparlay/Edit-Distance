#ifndef EDIT_DISTANCE_PARALLAL_H_
#define EDIT_DISTANCE_PARALLAL_H_

#include <string>
#include "parlay/sequence.h"

class EditDistanceParallel {
 public:
  size_t Solve(const parlay::sequence<uint32_t>& a, const parlay::sequence<uint32_t>& b);
};

#endif // EDIT_DISTANCE_PARALLAL_H_
