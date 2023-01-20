#include "hash_lcp_parallel.h"

// Returns an 32-bit integer of edit distance for two sequences A and B
// using bfs and hashing table query for lcp, in parallel
template <typename Seq>
int EditDistanceHashParallel(const Seq &a, const Seq &b) {
  // build sparse table

  /**
   * Record of building time
   */
  Timer tmr;
  int n = a.size();
  int m = b.size();
  if (n == 0) {
    return m;
  }
  if (m == 0) {
    return n;
  }
  parlay::sequence<int> logN1;
  parlay::sequence<int> powerN1;
  parlay::sequence<parlay::sequence<int>> table_s1;
  parlay::sequence<parlay::sequence<int>> table_s2;

  build_hash_table(a, b, table_s1, table_s2, powerN1, logN1);
  auto Diag = [&](int i, int j) { return i - j + m; };
  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);

  double building_tm = tmr.elapsed();
  std::cout << " building time of hashing: " << building_tm << std::endl;
  max_row[Diag(0, 0)] = query_lcp(a, b, table_s1, table_s2, logN1, 0, 0);
  // bfs for path
  int k = 0;
  for (;;) {
    if (max_row[Diag(n, m)] == n) break;  // find path
    k++;
    int l = Diag(0, std::min(k, int(m)));
    int r = Diag(std::min(k, int(n)), 0);
    parlay::parallel_for(l, r + 1, [&](int id) {
      int t = -1;
      if (max_row[id] != -1) {
        int i = max_row[id];
        int j = i + m - id;
        if (i == n || j == m) {
          t = i;
        } else {
          int get_lcp =
              query_lcp(a, b, table_s1, table_s2, logN1, i + 1, j + 1);
          t = i + 1 + get_lcp;
        }
      }
      if (id > 0 && max_row[id - 1] != -1) {
        int i = max_row[id - 1];
        int j = i + m - id + 1;
        if (i == n) {
          t = n;
        } else {
          t = std::max(
              t, i + 1 + query_lcp(a, b, table_s1, table_s2, logN1, i + 1, j));
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          t = std::max(
              t, i + query_lcp(a, b, table_s1, table_s2, logN1, i, j + 1));
        }
      }
      // assert(t <= n);
      temp[id] = t;
    });
    parlay::parallel_for(l, r + 1,
                         [&](int id) { max_row[id] = std::min(temp[id], id); });
  }
  return k;
}
