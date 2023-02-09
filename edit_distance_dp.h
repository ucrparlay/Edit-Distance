#ifndef EDIT_DISTANCE_DP_
#define EDIT_DISTANCE_DP_

#include <string>

#include "utils.h"

template <typename T>
size_t EditDistanceDP(const parlay::sequence<T>& a,
                      const parlay::sequence<T>& b) {
  size_t n = a.size(), m = b.size();
  constexpr uint32_t MAX_VAL = std::numeric_limits<uint32_t>::max() / 2;
  std::vector<std::vector<uint32_t>> dp(n + 1,
                                        std::vector<uint32_t>(m + 1, MAX_VAL));
  for (int i = 0; i <= n; i++) {
    dp[i][0] = i;
  }
  for (int j = 0; j <= m; j++) {
    dp[0][j] = j;
  }
  for (size_t i = 1; i <= n; i++) {
    for (size_t j = 1; j <= m; j++) {
      if (a[i - 1] == b[j - 1]) {
        dp[i][j] = dp[i - 1][j - 1];
      } else {
        dp[i][j] = std::min({dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1]}) + 1;
      }
    }
  }
  return dp[n][m];
}

#endif  // namespace EDIT_DISTANCE_DP_
