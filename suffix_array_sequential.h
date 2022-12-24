#ifndef SUFFIX_ARRAY_SEQUENTIAL_H_
#define SUFFIX_ARRAY_SEQUENTIAL_H_

#include <array>
#include <string>
#include <vector>

struct SuffixArraySequential {
  int n;
  std::vector<int> a, sa, rank, height;

  void Resort(std::vector<std::array<int, 2>>& key);
  void Build(const std::vector<int>& a_);
};

#endif  // namespace SUFFIX_ARRAY_SEQUENTIAL_H_
