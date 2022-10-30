#ifndef EDIT_DISTANCE_SEQUENTIAL_H_
#define EDIT_DISTANCE_SEQUENTIAL_H_

#include <vector>

class EditDistanceSequential {
 public:
  size_t Solve(const std::vector<int>& a, const std::vector<int>& b);
};

#endif  // EDIT_DISTANCE_SEQUENTIAL_H_
