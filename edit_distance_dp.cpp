#include "edit_distance_dp.h"

#include <algorithm>
#include <vector>

size_t EditDistanceDP::Solve(const std::string a, const std::string b) {
  int n = a.size(), m = b.size();
  std::vector<std::vector<int>> f(n + 1, std::vector<int>(m + 1, n + m));
  f[n][m] = 0;
  for (int i = n; i >= 0; i--) {
    for (int j = m; j >= 0; j--) {
      if (i < n) {
        f[i][j] = std::min(f[i][j], f[i + 1][j] + 1);
      }
      if (j < m) {
        f[i][j] = std::min(f[i][j], f[i][j + 1] + 1);
      }
      if (i < n && j < m) {
        if (a[i] == b[j]) {
          f[i][j] = std::min(f[i][j], f[i + 1][j + 1]);
        } else {
          f[i][j] = std::min(f[i][j], f[i + 1][j + 1] + 1);
        }
      }
    }
  }
  return f[0][0];
}