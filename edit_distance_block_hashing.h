#include "hash_block_parallel.h"

// Returns an 32-bit integer of edit distance for two sequences A and B
// using bfs and hashing table query for lcp, in parallel
template <typename Seq>
int EditDistanceBlockHashParallel(const Seq &a, const Seq &b)
{
  // build sparse table
  size_t n = a.size();
  size_t m = b.size();
  vector<vector<int>> table_A;
  vector<vector<int>> table_B;
  vector<int> S_A;
  vector<int> S_B;
  vector<int> P_A;
  vector<int> P_B;
  vector<int> aux_table;
  construct_table(a, b, table_A, table_B, P_A, P_B, S_A, S_B, aux_table, 34);
  auto Diag = [&](int i, int j)
  { return i - j + m; };
  parlay::sequence<int> max_row(n + m + 1, -1), temp(n + m + 1);
  max_row[Diag(0, 0)] = query_lcp(0, 0, a, b, table_A, table_B, S_A, P_A, S_B, P_B, aux_table);

  // bfs for path
  int k = 0;
  for (;;)
  {
    if (max_row[Diag(n, m)] == n)
      break; // find path
    k++;
    int l = Diag(0, std::min(k, int(m)));
    int r = Diag(std::min(k, int(n)), 0);
    parlay::parallel_for(l, r + 1, [&](int id)
                         {
      int t = -1;
      if (max_row[id] != -1) {
        int i = max_row[id];
        int j = i + m - id;
        if (i == n || j == m) {
          t = i;
        } else {
          int get_lcp =
              query_lcp(i + 1, j + 1, a, b, table_A, table_B, S_A, P_A, S_B, P_B, aux_table);
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
              t, i + 1 + query_lcp(i + 1, j, a, b, table_A, table_B, S_A, P_A, S_B, P_B, aux_table));
        }
      }
      if (id < n + m && max_row[id + 1] != -1) {
        int i = max_row[id + 1];
        int j = i + m - id - 1;
        if (j == m) {
          t = std::max(t, i);
        } else {
          t = std::max(
              t, i + query_lcp(i, j + 1, a, b, table_A, table_B, S_A, P_A, S_B, P_B, aux_table));
        }
      }
      // assert(t <= n);
      temp[id] = t; });
    parlay::parallel_for(l, r + 1,
                         [&](int id)
                         { max_row[id] = std::min(temp[id], id); });
  }
  return k;
}
