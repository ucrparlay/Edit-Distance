#include "edit_distance_sequential.h"

size_t EditDistanceSequential::Solve(const std::vector<int>& a,
                                     const std::vector<int>& b) {
  return a.size() + b.size();
}
