#include "edit_distance_parallel.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include "parlay/sequence.h"
#include "range_min.h"
#include "sparse_table_sequential.h"
#include "suffix_array_parallel.h"
#include "utils.h"

size_t EditDistanceSA(const parlay::sequence<uint32_t>& a,
                      const parlay::sequence<uint32_t>& b,
                      double* building_tm) {
  Timer tmr;
  int n = a.size(), m = b.size();
  auto c = parlay::sequence<uint32_t>(n + m + 1);
  parlay::parallel_for(0, n, [&](int i) { c[i] = a[i]; });
  parlay::parallel_for(0, m, [&](int i) { c[i + n + 1] = b[i]; });
  c[n] = std::numeric_limits<uint32_t>::max();
  auto rank = parlay::sequence<unsigned int>();
  auto sa = parlay::sequence<unsigned int>();
  auto lcp = parlay::sequence<unsigned int>();
  std::tie(rank, sa, lcp) = suffix_array_large_alphabet(c);
  auto rmq = range_min(lcp);
  auto GetLcp = [&](int i, int j) -> int {
    // std::cout << "GetLcp " << i << ' ' << j << '\n';
    if (i == n || j == m) return 0;
    assert(0 <= i && i < n && 0 <= j && j < m);
    for (int k = 0; k < 8; k++) {
      if (i + k >= n || j + k >= m || a[i + k] != b[j + k]) {
        return k;
      }
    }
    int l = rank[i], r = rank[j + n + 1];
    if (l > r) std::swap(l, r);
    assert(l < r);
    int id = rmq.query(l + 1, r);
    return lcp[id];
  };
  *building_tm = tmr.elapsed();
  // std::cout << " building time of SA: " << building_tm << std::endl;

  auto Diag = [&](int i, int j) { return i - j + m; };

  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);
  max_row[Diag(0, 0)] = GetLcp(0, 0);

  int k = 0;
  for (;;) {
    if (max_row[Diag(n, m)] == n) break;  // find path
    k++;
    int l = Diag(0, std::min(k, m)), r = Diag(std::min(k, n), 0);
    parlay::parallel_for(l, r + 1, [&](int id) {
      int t = -1;
      if (max_row[id] != -1) {
        int i = max_row[id];
        int j = i + m - id;
        if (i == n || j == m) {
          t = i;
        } else {
          t = i + 1 + GetLcp(i + 1, j + 1);
        }
      }
      if (id > 0 && max_row[id - 1] != -1) {
        int i = max_row[id - 1];
        int j = i + m - id + 1;
        if (i == n) {
          t = n;
        } else {
          t = std::max(t, i + 1 + GetLcp(i + 1, j));
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          t = std::max(t, i + GetLcp(i, j + 1));
        }
      }
      assert(t <= n);
      temp[id] = t;
    });
    parlay::parallel_for(l, r + 1,
                         [&](int id) { max_row[id] = std::min(temp[id], id); });
    // std::cout << "k: " << k << ", ";
    // for (int i = 0; i < n + m + 1; i++) std::cout << max_row[i] << ' ';
    // std::cout << std::endl;
  }
  return k;
}
