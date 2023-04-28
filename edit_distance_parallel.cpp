#include "edit_distance_parallel.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include "dc3.h"
#include "parlay/sequence.h"
#include "range_min.h"
#include "sparse_table_sequential.h"
#include "suffix_array_parallel.h"
#include "utils.h"

// #define COMPUTE_AVERAGE_LCP

size_t EditDistanceSA(const parlay::sequence<uint32_t>& a,
                      const parlay::sequence<uint32_t>& b, double* building_tm,
                      bool use_DC3) {
  Timer tmr;
#ifdef COMPUTE_AVERAGE_LCP
  std::atomic_uint64_t lcp_total = 0;
  std::atomic_uint64_t lcp_cnt = 0;
#endif
  int n = a.size(), m = b.size();
  auto c = parlay::sequence<uint32_t>(n + m);
  parlay::parallel_for(0, n, [&](int i) { c[i] = a[i]; });
  parlay::parallel_for(0, m, [&](int i) { c[i + n] = b[i]; });
  auto rank = parlay::sequence<unsigned int>();
  auto sa = parlay::sequence<unsigned int>();
  auto lcp = parlay::sequence<unsigned int>();
  if (use_DC3) {
    std::tie(rank, sa, lcp) = DC3(c);
  } else {
    std::tie(rank, sa, lcp) = suffix_array_large_alphabet(c);
  }
  auto rmq = range_min(lcp);
  auto GetLcp = [&](int i, int j) -> int {
    if (i == n || j == m) return 0;
    assert(0 <= i && i < n && 0 <= j && j < m);
    for (int k = 0; k < 8; k++) {
      if (i + k >= n || j + k >= m || a[i + k] != b[j + k]) {
        return k;
      }
    }
    int l = rank[i], r = rank[j + n];
    if (l > r) std::swap(l, r);
    int id = rmq.query(l + 1, r);
    return std::min(lcp[id], (unsigned int)n - i);
  };
  *building_tm = tmr.elapsed();
  // std::cout << " building time of SA: " << building_tm << std::endl;

  auto Diag = [&](int i, int j) { return i - j + m; };

  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);
  max_row[Diag(0, 0)] = GetLcp(0, 0);
#ifdef COMPUTE_AVERAGE_LCP
  lcp_total += GetLcp(0, 0);
  lcp_cnt++;
#endif
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
          int _lcp = GetLcp(i + 1, j + 1);
          t = i + 1 + _lcp;
#ifdef COMPUTE_AVERAGE_LCP
          lcp_total += _lcp;
          lcp_cnt++;
#endif
        }
      }
      if (id > 0 && max_row[id - 1] != -1) {
        int i = max_row[id - 1];
        int j = i + m - id + 1;
        if (i == n) {
          t = n;
        } else {
          int _lcp_2 = GetLcp(i + 1, j);
          t = std::max(t, i + 1 + _lcp_2);
#ifdef COMPUTE_AVERAGE_LCP
          lcp_total += _lcp_2;
          lcp_cnt++;
#endif
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          int _lcp_3 = GetLcp(i, j + 1);
          t = std::max(t, i + _lcp_3);
#ifdef COMPUTE_AVERAGE_LCP
          lcp_total += _lcp_3;
          lcp_cnt++;
#endif
        }
      }
      assert(t <= n);
      temp[id] = t;
    });
    parlay::parallel_for(l, r + 1,
                         [&](int id) { max_row[id] = std::min(temp[id], id); });
  }
#ifdef COMPUTE_AVERAGE_LCP
  std::cout << "Lcp total: " << lcp_total << std::endl;
  std::cout << "lcp_cnt: " << lcp_cnt << std::endl;
  std::cout << "average: " << (double)(lcp_total) / lcp_cnt << std::endl;
#endif
  return k;
}
