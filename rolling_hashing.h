#ifndef ROLLING_HASHING_H
#define ROLLING_HASHING_H
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "utils.h"

using namespace parlay;
// static constexpr int PRIME = 479;
using hash_r_T = uint64_t;

hash_r_T q_power(hash_r_T base, size_t n) {
  hash_r_T ret = 1;
  hash_r_T a = base;
  while (n) {
    if (n & 1) {
      ret = ret * a;
    }
    a = a * a;
    n >>= 1;
  }
  return ret;
}

/**
 * Inplace inclusive scan
 */
template <typename Seq>
void s_inplace_scan_inclusive(Seq &A, size_t n) {
  auto block_size = std::max((size_t)(5000), (size_t)std::sqrt(n));
  if (n <= 100000) {
    for (size_t i = 1; i < n; i++) {
      A[i] += A[i - 1] * PRIME;
    }
  } else {
    hash_r_T power_table[block_size];
    power_table[0] = PRIME;
    for (int i = 1; i < block_size; i++) {
      power_table[i] = (hash_r_T)(PRIME * power_table[i - 1]);
    }
    size_t num_blocks = (n - 1) / block_size + 1;
    hash_r_T p_sum[num_blocks];

    parlay::parallel_for(0, num_blocks, [&](size_t i) {
      for (size_t j = i * block_size + 1; j < std::min((i + 1) * block_size, n);
           j++) {
        A[j] += A[j - 1] * PRIME;
      }
    });
    p_sum[0] = A[block_size - 1];
    for (size_t k = 1; k < num_blocks; k++) {
      p_sum[k] = A[(k + 1) * block_size - 1] +
                 p_sum[k - 1] * power_table[block_size - 1];
    }

    parlay::parallel_for(1, num_blocks, [&](size_t i) {
      for (size_t j = i * block_size;
           j < std::min(i * block_size + block_size, n); j++) {
        A[j] += p_sum[j / block_size - 1] *
                power_table[j - block_size * (j / block_size)];
      }
    });
  }
}

template <typename Seq>
void build_rolling(const Seq &s1, const Seq &s2,
                   parlay::sequence<hash_r_T> &table_s1,
                   parlay::sequence<hash_r_T> &table_s2) {
  size_t table1_size = s1.size();
  size_t table2_size = s2.size();
  table_s1.resize(table1_size);
  table_s2.resize(table2_size);

  // prefix sum
  // for table_1
  parlay::parallel_for(0, table1_size,
                       [&](uint32_t i) { table_s1[i] = (hash_r_T)(s1[i]); });
  s_inplace_scan_inclusive(table_s1, table1_size);

  // for table 2
  parlay::parallel_for(0, table2_size,
                       [&](uint32_t i) { table_s2[i] = (hash_r_T)(s2[i]); });
  s_inplace_scan_inclusive(table_s2, table2_size);
}

// query hash value from i to j
hash_r_T get_hash(const parlay::sequence<hash_r_T> &hash_table, size_t i,
                  size_t j) {
  hash_r_T value_s;
  if (i == 0) {
    value_s = 0;
  } else {
    value_s = hash_table[i - 1];
  }
  hash_r_T value_t = hash_table[j];
  hash_r_T pw_diff = q_power(PRIME, j - i + 1);
  hash_r_T res = value_t - value_s * pw_diff;
  return res;
}

template <typename T>
int query_rolling(const parlay::sequence<T> &s1, const parlay::sequence<T> &s2,
                  const parlay::sequence<hash_r_T> &table1,
                  const parlay::sequence<hash_r_T> &table2, size_t i,
                  size_t j) {
  if ((uint32_t)i >= s1.size() || (uint32_t)j >= s2.size()) return 0;
  if ((hash_r_T)s1[i] != (hash_r_T)s2[j]) return 0;
  uint32_t try_r = 1;
  uint32_t r = std::min(s1.size() - i, s2.size() - j);
  uint32_t l = 0;
  while (try_r <= r &&
         get_hash(table1, i, i + try_r) == get_hash(table2, j, j + try_r)) {
    l = try_r;
    try_r *= 2;
  }

  r = std::min(try_r, r);

  uint32_t res = 0;
  while (l <= r) {
    uint32_t m = l + (r - l) / 2;
    if (get_hash(table1, i, i + m - 1) == get_hash(table2, j, j + m - 1)) {
      res = m;
      l = m + 1;
    } else {
      r = m - 1;
    }
  }
  return res;
}

#endif
