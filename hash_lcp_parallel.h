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
void build_hash_table(const parlay::sequence<T> &s1,
                      const parlay::sequence<T> &s2,
                      vector<vector<int>> &table_s1,
                      vector<vector<int>> &table_s2, vector<int> &powerN1,
                      vector<int> &powerN2, vector<int> &logN1,
                      vector<int> &logN2)
{
  // build the first leaves layer
  ull table1_d2 = s1.size();
  ull table2_d2 = s2.size();
  ull table1_d1 = 0;
  int upp1 = 1;
  ull table2_d1 = 0;
  int upp2 = 1;
  while (upp1 <= table1_d2)
  {
    table1_d1++;
    upp1 *= 2;
  }
  while (upp2 < table2_d2)
  {
    table2_d1++;
    upp2 *= 2;
  }
  // initialize size
  table_s1.resize(table1_d1);
  table_s2.resize(table2_d1);
  table_s1[0].resize(table1_d2);
  table_s2[0].resize(table2_d2);
  powerN1.resize(table1_d1 + 1);
  powerN2.resize(table2_d1 + 1);
  logN1.resize(table1_d1);
  logN2.resize(table2_d1);
  // build the first layer of ST in parallel
  parlay::parallel_for(0, table1_d2,
                       [&](int i)
                       { table_s1[0][i] = int(s1[i]); });
  parlay::parallel_for(0, table2_d2,
                       [&](int i)
                       { table_s2[0][i] = int(s2[i]); });
  int len = 1;
  powerN1[0] = 1;
  for (int i = 0; i < table1_d1; i++)
  {
    logN1[i] = len;
    len *= 2;
    powerN1[i + 1] = (PRIME_BASE << i);
  }
  len = 1;
  powerN2[0] = 1;
  for (int i = 0; i < table2_d1; i++)
  {
    logN2[i] = len;
    len *= 2;
    powerN2[i + 1] = (PRIME_BASE << i);
  }

  // build the second to k-th layer
  for (int i = 1; i < table1_d1; i++)
  {
    table_s1[i].resize(table1_d2 - logN1[i] + 1);
    parlay::parallel_for(0, table1_d2 - logN1[i] + 1, [&](int j)
                         {
      for (int pw = 0; pw < logN1[i - 1]; pw++) {
        table_s1[i][j] = table_s1[i - 1][j] * PRIME_BASE;
      }
      table_s1[i][j] += table_s1[i - 1][j + logN1[i - 1]]; });
  }
  for (int i = 1; i < table2_d1; i++)
  {
    table_s2[i].resize(table2_d2 - logN2[i] + 1);
    parlay::parallel_for(0, table2_d2 - logN2[i] + 1, [&](int j)
                         {
      for (int pw = 0; pw < logN2[i - 1]; pw++) {
        table_s2[i][j] = table_s2[i - 1][j] * PRIME_BASE;
      }
      table_s2[i][j] += table_s2[i - 1][j + logN2[i - 1]]; });
  }
}

// Returns an 32-bit integer for each query range, from `i` to `j`,
// for the longest common prefix sequence of global sequences A and B
// Pass the two hash value tables for A and B, and the logn scale table
// from 1 to logn.
//
// This method is equivalent to function:
//    auto lcp(Seq1 const &s, Seq2 const &SA);
int query_lcp(vector<vector<int>> &table1, vector<vector<int>> &table2,
              vector<int> &logN1, vector<int> &logN2, int i, int j)
{
  if (table1[0][i] != table2[0][j])
    return 0;
  if (i == table1[0].size() || i == table1[0].size())
    return 1;
  int v1;
  int v2;
  // possible value is from 0 to the smaller one of the remaining sequences
  int max_range = std::min(table1[0].size() - i, table2[0].size() - j);
  // search from the power first
  for (int log_i = logN1.size() - 1; log_i != 0; log_i--)
  {
    if (logN1[log_i] <= max_range)
    {
      v1 = table1[log_i][i];
      v2 = table2[log_i][j];
      if (v1 == v2)
      {
        return (logN1[log_i] + query_lcp(table1, table2, logN1, logN2,
                                         i + logN1[log_i], j + logN1[log_i]));
      }
    }
  }
  return 1;
}
#endif
