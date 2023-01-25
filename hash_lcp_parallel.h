#ifndef lcp_parallel_h
#define lcp_parallel_h

#include "hash_block_parallel.h"
#include "utils.h"
#define ull unsigned long long
using namespace std;

// Function for building two sparse hash tables for sequences A and B
// Input: two seuences `s1` and `s2`, two 2-d vectors as tables.
// Pre-computed Power tables [p^0, p^1, ...]. logn scale table
// from 1 to logn.
template <typename T>
void build_hash_table(const parlay::sequence<T> &s1,
                      const parlay::sequence<T> &s2,
                      parlay::sequence<parlay::sequence<int>> &table_s1,
                      parlay::sequence<parlay::sequence<int>> &table_s2,
                      parlay::sequence<int> &powerN1,
                      parlay::sequence<int> &logN1) {
  // build the first leaves layer
  size_t table1_d2 = s1.size();
  size_t table2_d2 = s2.size();
  size_t table1_d1 = 0;
  size_t table2_d1 = 0;
  table1_d1 = FASTLOG2(table1_d2) + 1;
  table2_d1 = FASTLOG2(table2_d2) + 1;
  // initialize size
  table_s1.resize(table1_d1);
  table_s2.resize(table2_d1);
  table_s1[0].resize(table1_d2);
  table_s2[0].resize(table2_d2);
  // build the first layer of ST in parallel
  parlay::parallel_for(0, table1_d2,
                       [&](int i) { table_s1[0][i] = int(s1[i]); });
  parlay::parallel_for(0, table2_d2,
                       [&](int j) { table_s2[0][j] = int(s2[j]); });

  int aux_log_size = std::max(table1_d1, table2_d1);
  logN1.resize(aux_log_size);
  powerN1.resize(aux_log_size + 1);
  int len = 1;
  powerN1[0] = 1;
  powerN1[1] = PRIME_BASE;
  for (size_t i = 0; i < aux_log_size; i++) {
    logN1[i] = len;
    len *= 2;
  }
  for (size_t i = 1; i < aux_log_size; i++) {
    powerN1[i + 1] = powerN1[i] * powerN1[i];
  }

  // build the second to k-th layer
  if (table1_d1 > 1) {
    for (size_t i = 1; i < table1_d1; i++) {
      table_s1[i].resize(table1_d2 - logN1[i] + 1);
      parlay::parallel_for(0, table1_d2 - logN1[i] + 1, [&](int j) {
        table_s1[i][j] =
            table_s1[i - 1][j + logN1[i - 1]] + table_s1[i - 1][j] * powerN1[i];
      });
    }
  }
  if (table2_d1 > 1) {
    for (size_t i = 1; i < table2_d1; i++) {
      table_s2[i].resize(table2_d2 - logN1[i] + 1);
      parlay::parallel_for(0, table2_d2 - logN1[i] + 1, [&](int j) {
        table_s2[i][j] =
            table_s2[i - 1][j + logN1[i - 1]] + table_s2[i - 1][j] * powerN1[i];
      });
    }
  }
}

// Returns an 32-bit integer for each query range, from `i` to `j`,
// for the longest common prefix sequence of global sequences A and B
// Pass the two hash value tables for A and B, and the logn scale table
// from 1 to logn.
//
// This method is equivalent to function:
//    auto lcp(Seq1 const &s, Seq2 const &SA);
template <typename T>
int query_lcp(const parlay::sequence<T> &s1, const parlay::sequence<T> &s2,
              parlay::sequence<parlay::sequence<int>> &table1,
              parlay::sequence<parlay::sequence<int>> &table2,
              parlay::sequence<int> &logN1, int i, int j) {
  if (i >= (int)(s1.size()) || j >= (int)(s2.size())) {
    return 0;
  }

  if (s1[i] != s2[j]) {
    return 0;
  }
  int ii = i;
  int jj = j;
  int longest_try_power = 0;
  for (int p = 0; p < logN1.size() + 1; p++) {
    if (p >= logN1.size() || ii + logN1[p] - 1 >= table1[0].size() ||
        jj + logN1[p] - 1 >= table2[0].size() ||
        table1[p][ii] != table2[p][jj]) {
      longest_try_power = p;
      break;
    }
  }

  longest_try_power--;
  if (longest_try_power == 0) {
    return 1;
  }
  int res = 0;
  ii += logN1[longest_try_power];
  jj += logN1[longest_try_power];
  res += logN1[longest_try_power];
  longest_try_power--;

  while (longest_try_power >= 0) {
    if ((ii + (logN1[longest_try_power]) - 1) < table1[0].size() &&
        (jj + (logN1[longest_try_power]) - 1 < table2[0].size())) {
      if (table1[longest_try_power][ii] == table2[longest_try_power][jj]) {
        ii += (logN1[longest_try_power]);
        jj += (logN1[longest_try_power]);
        res += logN1[longest_try_power];
      }
    }
    longest_try_power--;
  }
  if (res != 0) {
    return res;
  } else {
    return 1;
  }
}
#endif
