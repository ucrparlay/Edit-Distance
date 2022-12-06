#ifndef lcp_parallel_h
#define lcp_parallel_h

#include "utils.h"
#define ull unsigned long long
using namespace std;
const int PRIME_BASE = 139;

// Function for building two sparse hash tables for sequences A and B
// Input: two seuences `s1` and `s2`, two 2-d vectors as tables.
// Pre-computed Power tables [p^0, p^1, ...]. logn scale table
// from 1 to logn.
template <typename T>
void build_hash_table(const parlay::sequence<T>& s1,
                      const parlay::sequence<T>& s2,
                      vector<vector<int>>& table_s1,
                      vector<vector<int>>& table_s2, vector<int>& powerN1,
                      vector<int>& powerN2, vector<int>& logN1,
                      vector<int>& logN2) {
  // build the first leaves layer
  ull table1_d2 = s1.size();
  ull table2_d2 = s2.size();
  ull table1_d1 = 0;
  int upp1 = 1;
  ull table2_d1 = 0;
  int upp2 = 1;
  while (upp1 <= table1_d2) {
    table1_d1++;
    upp1 *= 2;
  }
  while (upp2 < table2_d2) {
    table2_d1++;
    upp2 *= 2;
  }
  // initialize size
  table_s1.resize(table1_d1);
  table_s2.resize(table2_d1);
  for (int i = 0; i < table1_d1; i++) table_s1[i].resize(table1_d2);
  for (int i = 0; i < table2_d1; i++) table_s2[i].resize(table2_d2);
  powerN1.resize(table1_d1);
  powerN2.resize(table2_d1);
  logN1.resize(table1_d1);
  logN2.resize(table2_d1);
  // build the first layer of ST in parallel
  parlay::parallel_for(0, table1_d2,
                       [&](int i) { table_s1[0][i] = int(s1[i]); });
  parlay::parallel_for(0, table2_d2,
                       [&](int i) { table_s2[0][i] = int(s2[i]); });
  int len = 1;
  powerN1[0] = 1;
  for (int i = 0; i < table1_d1; i++) {
    logN1[i] = len;
    len *= 2;
    powerN1[i + 1] = (PRIME_BASE << i);
  }
  len = 1;
  powerN2[0] = 1;
  for (int i = 0; i < table2_d1; i++) {
    logN2[i] = len;
    len *= 2;
    powerN2[i + 1] = (PRIME_BASE << i);
  }

  // build the second to k-th layer
  for (int i = 1; i < table1_d1; i++) {
    parlay::parallel_for(0, table1_d2 - logN1[i] + 1, [&](int j) {
      for (int pw = 0; pw < logN1[i - 1]; pw++) {
        table_s1[i][j] = table_s1[i - 1][j] * PRIME_BASE;
      }
      table_s1[i][j] += table_s1[i - 1][j + logN1[i - 1]];
    });
  }
  for (int i = 1; i < table2_d1; i++) {
    parlay::parallel_for(0, table2_d2 - logN2[i] + 1, [&](int j) {
      for (int pw = 0; pw < logN2[i - 1]; pw++) {
        table_s2[i][j] = table_s2[i - 1][j] * PRIME_BASE;
      }
      table_s2[i][j] += table_s2[i - 1][j + logN2[i - 1]];
    });
  }
}

// Returns an 32-bit integer for each query range, from `i` to `j`,
// for the longest common prefix sequence of global sequences A and B
// Pass the two hash value tables for A and B, and the logn scale table
// from 1 to logn.
//
// This method is equivalent to function:
//    auto lcp(Seq1 const &s, Seq2 const &SA);
int query_lcp(vector<vector<int>>& table1, vector<vector<int>>& table2,
              vector<int>& logN1, vector<int>& logN2, int i, int j) {
  int v1;
  int v2;
  // possible value is from 0 to the smaller one of the remaining sequences
  int pos_i = i;
  int pos_j = j;
  int l = 0;
  int r = min(table1[0].size() - i, table2[0].size() - j);
  // search from the power first
  int mid_try = 1;
  int power_id = 0;
  int power_try;
  while (l != r) {
    for (int i = 0; i < logN1.size(); i++) {
      if (l + logN1[i] <= r) {
        power_id = i;
        power_try = logN1[i];
        mid_try = l + power_try;
      }
    }
    // test the hash value (whether equal)
    v1 = table1[power_id][pos_i];
    v2 = table2[power_id][pos_j];
    if (v1 == v2) {
      pos_i += power_try;
      pos_j += power_try;
      l = mid_try;
    } else {
      r = mid_try - 1;
    }
  }
  return l;
}

#endif
