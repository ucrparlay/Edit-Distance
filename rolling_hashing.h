#ifndef ROLLING_HASHING_H
#define ROLLING_HASHING_H
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "utils.h"

using namespace parlay;
// static constexpr int PRIME = 479;

uint32_t q_power(uint32_t base, size_t n) {
  uint32_t ret = 1;
  uint32_t a = base;
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
  auto block_size = std::max((size_t)(1000), (size_t)std::sqrt(n));
  if (n <= 1000000) {
    for (size_t i = 1; i < n; i++) {
      A[i] += A[i - 1] * PRIME;
    }
  } else {
    uint32_t power_table[block_size];
    power_table[0] = PRIME;
    for (int i = 1; i < block_size; i++) {
      power_table[i] = (uint32_t)(PRIME * power_table[i - 1]);
    }
    size_t num_blocks = (n - 1) / block_size + 1;
    uint32_t p_sum[num_blocks];

    parlay::parallel_for(0, num_blocks, [&](size_t i) {
      // aux[i] = inplace_seq_scan_exclusive_direct(
      //     A, l + i * block_size, std::min(l + i * block_size + block_size,
      //     r));
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
                   parlay::sequence<uint32_t> &table_s1,
                   parlay::sequence<uint32_t> &table_s2) {
  size_t table1_size = s1.size();
  size_t table2_size = s2.size();
  table_s1.resize(table1_size);
  table_s2.resize(table2_size);

  // prefix sum
  // for table_1
  parlay::parallel_for(0, table1_size,
                       [&](uint32_t i) { table_s1[i] = (uint32_t)(s1[i]); });
  s_inplace_scan_inclusive(table_s1, table1_size);

  // for table 2
  parlay::parallel_for(0, table2_size,
                       [&](int i) { table_s2[i] = (uint32_t)(s2[i]); });
  s_inplace_scan_inclusive(table_s2, table2_size);
}

// query hash value from i to j
uint32_t get_hash(const parlay::sequence<uint32_t> &hash_table, size_t i,
                  size_t j) {
  uint32_t value_s = hash_table[i - 1];
  uint32_t value_t = hash_table[j];
  uint32_t pw_diff = q_power(PRIME, j - i + 1);
  uint32_t res = value_t - value_s * pw_diff;
  return res;
}

template <typename T>
int query_rolling(const parlay::sequence<T> &s1, const parlay::sequence<T> &s2,
                  const parlay::sequence<uint32_t> &table1,
                  const parlay::sequence<uint32_t> &table2, size_t i,
                  size_t j) {
  if ((uint32_t)s1[i] != (uint32_t)s2[j]) return 0;
  int try_r = 1;
  int r = std::min(s1.size() - i, s2.size() - j) - 1;
  while (get_hash(table1, i, i + try_r) == get_hash(table2, j, j + try_r) &&
         try_r <= r) {
    try_r *= 2;
  }

  r = try_r;
  int l = try_r / 2;
  int res = 0;
  if (r == 0) {
    return 1;
  }
  while (l <= r) {
    int m = l + (r - l) / 2;
    if (get_hash(table1, i, i + m) == get_hash(table2, j, j + m)) {
      res = m + 1;
      l = m + 1;
    } else {
      r = m - 1;
    }
  }
  return res;
}

#endif
