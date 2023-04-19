#ifndef ROLLING_HASHING_BLK_H
#define ROLLING_HASHING_BLK_H
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "utils.h"

using namespace parlay;
using hash_r_b_T = uint64_t;
size_t compression_ratio = 4;

hash_r_b_T q_power_blk(hash_r_b_T base, size_t n) {
  hash_r_b_T ret = 1;
  hash_r_b_T a = base;
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
void s_inplace_scan_inclusive_blk(Seq &A, size_t n) {
  auto block_size = std::max((size_t)(5000), (size_t)std::sqrt(n));
  // size_t compression_ratio = sizeof(hash_r_b_T) / sizeof(char);
  hash_r_b_T real_prime = q_power_blk(PRIME, compression_ratio);
  if (n <= 100000) {
    for (size_t i = 1; i < n; i++) {
      A[i] += A[i - 1] * real_prime;
    }
  } else {
    hash_r_b_T power_table[block_size];
    power_table[0] = real_prime;
    for (int i = 1; i < block_size; i++) {
      power_table[i] = (hash_r_b_T)(real_prime * power_table[i - 1]);
    }
    size_t num_blocks = (n - 1) / block_size + 1;
    hash_r_b_T p_sum[num_blocks];

    parlay::parallel_for(0, num_blocks, [&](size_t i) {
      for (size_t j = i * block_size + 1; j < std::min((i + 1) * block_size, n);
           j++) {
        A[j] += A[j - 1] * real_prime;
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
void build_rolling_blk(const Seq &s1, const Seq &s2,
                       parlay::sequence<hash_r_b_T> &table_s1,
                       parlay::sequence<hash_r_b_T> &table_s2) {
  // size_t compression_ratio = sizeof(hash_r_b_T) / sizeof(char);
  size_t table1_size = s1.size() / compression_ratio;
  size_t table2_size = s2.size() / compression_ratio;
  table_s1.resize(table1_size);
  table_s2.resize(table2_size);
  // for table_1
  // reload index

  // cout << "table 1 size: " << table1_size << endl;
  // cout << "table 2 size: " << table2_size << endl;

  parlay::parallel_for(0, table1_size, [&](uint32_t i) {
    hash_r_b_T res = (hash_r_b_T)(s1[i * compression_ratio]);
    // assert(res <= 256 && res >= 0);
    for (int j = i * compression_ratio + 1; j < (i + 1) * compression_ratio;
         j++) {
      res = res * PRIME + (hash_r_b_T)(s1[j]);
    }
    table_s1[i] = res;
  });
  s_inplace_scan_inclusive_blk(table_s1, table1_size);

  // for table 2
  parlay::parallel_for(0, table2_size, [&](uint32_t i) {
    hash_r_b_T res = (hash_r_b_T)(s2[i * compression_ratio]);
    // assert(res <= 256 && res >= 0);
    for (int j = i * compression_ratio + 1; j < (i + 1) * compression_ratio;
         j++) {
      res = res * PRIME + (hash_r_b_T)(s2[j]);
    }
    table_s2[i] = res;
  });
  s_inplace_scan_inclusive_blk(table_s2, table2_size);
}

// query hash value at point i (include)
template <typename Seq>
hash_r_b_T point_hash_blk(const Seq &s,
                          const parlay::sequence<hash_r_b_T> &b_hash_table,
                          size_t i) {
  // size_t compression_ratio = sizeof(hash_r_b_T) / sizeof(char);

  // if ((i + 1) % compression_ratio == 0) {
  //   return b_hash_table[i / compression_ratio];
  // } else {
  if (i < compression_ratio) {
    size_t num_remain = i + 1;
    hash_r_b_T value_blk = 0;
    for (size_t r = 0; r < num_remain; r++) {
      // assert(r < s.size());
      value_blk = value_blk * PRIME + (hash_r_b_T)(s[r]);
    }
    return value_blk;
  } else {
    size_t pos = (size_t)(i / compression_ratio - 1);
    // cout << "i: " << i << " ratio: " << compression_ratio << endl;
    size_t num_remain = i - (pos + 1) * compression_ratio + 1;
    // cout << "pos: " << pos << " size: " << b_hash_table.size() << endl;
    // assert(pos < b_hash_table.size());
    hash_r_b_T value_blk = b_hash_table[pos];
    for (size_t r = 0; r < num_remain; r++) {
      assert((pos + 1) * compression_ratio + r < s.size());
      value_blk = value_blk * PRIME +
                  (hash_r_b_T)(s[(pos + 1) * compression_ratio + r]);
    }
    return value_blk;
  }
  // }
}

// query hash value from i to j
template <typename Seq>
hash_r_b_T get_inter_hash_blk(const Seq &s,
                              const parlay::sequence<hash_r_b_T> &b_hash_table,
                              size_t i, size_t j) {
  hash_r_b_T value_s;

  if (i == 0) {
    value_s = 0;
  } else {
    value_s = point_hash_blk(s, b_hash_table, i - 1);
  }
  hash_r_b_T value_t = point_hash_blk(s, b_hash_table, j);
  hash_r_b_T pw_diff = q_power_blk(PRIME, j - i + 1);
  hash_r_b_T res = value_t - value_s * pw_diff;
  return res;
}

template <typename T>
int query_rolling_blk(const parlay::sequence<T> &s1,
                      const parlay::sequence<T> &s2,
                      const parlay::sequence<hash_r_b_T> &table1,
                      const parlay::sequence<hash_r_b_T> &table2, size_t i,
                      size_t j) {
  if ((size_t)i >= s1.size() || (size_t)j >= s2.size()) return 0;
  if ((hash_r_b_T)s1[i] != (hash_r_b_T)s2[j]) return 0;
  size_t try_r = 1;
  size_t r = std::min(s1.size() - i, s2.size() - j);
  size_t l = 0;
  while (try_r <= r && get_inter_hash_blk(s1, table1, i, i + try_r) ==
                           get_inter_hash_blk(s2, table2, j, j + try_r)) {
    l = try_r;
    try_r *= 2;
  }

  r = std::min(try_r, r);
  // std::cout << "l: " << l << " "
  //           << "r: " << r << std::endl;
  size_t res = 0;
  while (l <= r) {
    size_t m = l + (r - l) / 2;
    if (get_inter_hash_blk(s1, table1, i, i + m - 1) ==
        get_inter_hash_blk(s2, table2, j, j + m - 1)) {
      res = m;
      l = m + 1;
    } else {
      r = m - 1;
    }
  }
  // std::cout << "res: " << res << std::endl;
  return res;
}

#endif
