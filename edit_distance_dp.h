#ifndef EDIT_DISTANCE_DP_
#define EDIT_DISTANCE_DP_

#include <string>

#include "utils.h"

template <typename T>
class EditDistanceDP {
 public:
  size_t Solve(const parlay::sequence<T>& a, const parlay::sequence<T>& b);
};

#endif  // namespace EDIT_DISTANCE_DP_
