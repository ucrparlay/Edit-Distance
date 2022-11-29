#ifndef EDIT_DISTANCE_DP_
#define EDIT_DISTANCE_DP_

#include <string>
#include "utils.h"

class EditDistanceDP {
 public:
  size_t Solve(const parlay::sequence<uint32_t>& a, const parlay::sequence<uint32_t>& b);
};

#endif  // namespace EDIT_DISTANCE_DP_
