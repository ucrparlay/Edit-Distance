#include "edit_distance_dp.h"

size_t EditDistanceDP::Solve(const std::vector<int>& a,
                             const std::vector<int>& b) {
  return a.size() + b.size();
}