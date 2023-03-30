#ifndef ROLLING_HASHING_H
#define ROLLING_HASHING_H
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "utils.h"

using namespace parlay;
static constexpr int PRIME = 479;

// a binary associative operator
auto f = [](const uint32_t &a, const uint32_t &b) -> uint32_t {
  return a * PRIME + b;
};
auto mul_add = parlay::binary_op(f, 0);

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
  // size_t block_size_1 = std::sqrt(table1_size);

  // auto ress = parlay::scan_inclusive(table_s1, mul_add);
  for (int j = 1; j < table1_size; j++) {
    table_s1[j] = (uint32_t)(table_s1[j - 1] * PRIME + table_s1[j]);
  }

  // for table 2
  parlay::parallel_for(0, table2_size,
                       [&](int i) { table_s2[i] = int(s2[i]); });
  // size_t block_size_1 = std::sqrt(table1_size);
  for (int j = 1; j < table2_size; j++) {
    table_s2[j] += table_s2[j - 1] * PRIME;
  }
  // parlay::scan_inclusive_inplace(parlay::make_slice(table_s2), mul_add);
}

// query hash value from i to j
int get_hash(const parlay::sequence<uint32_t> &hash_table, size_t i, size_t j) {
  int value_s = hash_table[i - 1];
  int value_t = hash_table[j];
  int pw_diff = quick_power(PRIME, j - i + 1);
  int res = value_t - value_s * pw_diff;
  return res;
}

template <typename T>
int query_rolling(const parlay::sequence<T> &s1, const parlay::sequence<T> &s2,
                  const parlay::sequence<uint32_t> &table1,
                  const parlay::sequence<uint32_t> &table2, size_t i,
                  size_t j) {
  if (s1[i] != s2[j]) return 0;
  int r = std::min(s1.size() - i, s2.size() - j);
  // std::cout << "r value: " << r << std::endl;

  int l = 0;
  int res = 0;
  if (r == 1) {
    return 1;
  }
  while (l <= r) {
    int m = l + (r - l) / 2;
    if (get_hash(table1, i, i + m) == get_hash(table2, j, j + m)) {
      res = m + 1;
      l = m + 1;
      // std::cout << "l value: " << l << std::endl;
    } else {
      r = m - 1;
      // std::cout << "r value: " << r << std::endl;
    }
  }
  return res;
}

#endif
